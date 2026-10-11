# =============================================================================
# 02b_slc_cox_screen.R
# -----------------------------------------------------------------------------
# PURPOSE
#   Multivariable Cox screen of the solute carrier (SLC) family in TCGA-LGG and
#   TCGA-GBM. The screen locates SLC7A7 within its own gene family; it is not a
#   test of SLC7A7. Inference is in 04_multicohort_validation.R.
#
# INPUT   cohorts loaded through load_cohort() (R/cohorts.R), which harmonises
#         clinical tables to: sample, time, event, age, grade, idh, codel, mgmt
# OUTPUT  data/processed/LGG_all_cox_zscore.csv   -> Additional file 1, sheet TCGA_LGG
#         data/processed/GBM_all_cox_zscore.csv   -> Additional file 1, sheet TCGA_GBM
#
# DESIGN
#   * Genes are selected by symbol match on GENE_FAMILY = "SLC" (symbol contains
#     the string), not from a curated panel. Pseudogenes, antisense and
#     read-through loci carrying the string are therefore included and are
#     visible in the output; on the TCGA annotation used here all 457 matching
#     symbols also begin with SLC.
#   * Expression enters as a cohort-internal z-score of the log expression
#     value, so the hazard ratio is per standard deviation.
#   * Covariate sets differ by cohort and by stratum, because a covariate that
#     does not vary inside a cohort or stratum cannot be fitted there. The sets
#     below are the ones used for the published screen (SCREEN_SPEC). The
#     reasons they differ:
#         - WHO grade is constant in TCGA-GBM (all grade IV), so it is absent
#           from every TCGA-GBM model.
#         - IDH status is constant inside each IDH stratum, so it appears only
#           in the all-patient models.
#         - 1p/19q codeletion is essentially restricted to IDH-mutant glioma,
#           so it is not fitted inside the IDH-wildtype stratum.
#   * A model is fitted only when the stratum reaches the event threshold in
#     SCREEN_SPEC. Those thresholds correspond to 10 events per variable
#     throughout, except the TCGA-GBM all-patient model, whose threshold of 20
#     events corresponds to 6.7 per variable; that cohort carries far more
#     events than either threshold, so the gate is not active there.
#   * Under monotone likelihood (separation) the fit falls back to Firth
#     penalised Cox regression (coxphf); if that also fails the gene carries NA.
#     Genes that carry NA in the published tables are pseudogenes and antisense
#     loci whose expression does not support a fit.
#   * Benjamini-Hochberg correction is applied across the screened genes within
#     each cohort and stratum.
#
# EXPECTED OUTPUT (matches the published Additional file 1)
#   TCGA-LGG  457 screened; fitted 455 / 452 / 450 (all patients / mutant / wildtype)
#   TCGA-GBM  457 screened; fitted 456 /   0 / 456 - the IDH-mutant stratum of
#             TCGA-GBM never reaches 20 events and is not fitted for any gene.
#   SLC7A7    TCGA-LGG ranks 13 / 455, 29 / 452, 17 / 450
#             TCGA-GBM ranks  8 / 456, not fitted, 17 / 456
# =============================================================================

source(here::here("R", "setup.R"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(survival); library(coxphf)
})

# ---- covariate sets and event thresholds actually used -----------------------
# rhs        right-hand side of the Cox model; x is the z-scored gene
# min_events minimum number of deaths required in the stratum before fitting
SCREEN_SPEC <- list(
  TCGA_LGG = list(
    overall  = list(rhs = "x + age + idh + codel + grade", min_events = 50),
    mutant   = list(rhs = "x + age + codel + grade",       min_events = 40),
    wildtype = list(rhs = "x + age + grade",               min_events = 30)
  ),
  TCGA_GBM = list(
    overall  = list(rhs = "x + age + mgmt",                min_events = 20),
    mutant   = list(rhs = "x + age",                       min_events = 20),
    wildtype = list(rhs = "x + age",                       min_events = 20)
  )
)

# ---- Cox with Firth fallback -------------------------------------------------
safe_cox_screen <- function(formula, data) {
  warn_msg <- NULL
  fit <- withCallingHandlers(
    tryCatch(coxph(formula, data = data), error = function(e) NULL),
    warning = function(w) { warn_msg <<- conditionMessage(w); invokeRestart("muffleWarning") }
  )

  if (!is.null(warn_msg) && grepl("infinite", warn_msg)) {      # separation
    fit <- tryCatch(coxphf(formula, data = data), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    return(data.frame(HR = exp(fit$coefficients[1]), lower95 = exp(fit$ci.lower[1]),
                      upper95 = exp(fit$ci.upper[1]), p = fit$prob[1], model = "firth"))
  }
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  data.frame(HR = s$conf.int[1, "exp(coef)"], lower95 = s$conf.int[1, "lower .95"],
             upper95 = s$conf.int[1, "upper .95"], p = s$coefficients[1, "Pr(>|z|)"],
             model = "cox")
}

# ---- one gene, three strata --------------------------------------------------
analyze_gene <- function(gene, expr, clin, spec) {

  d <- clin
  d$x <- scale(as.numeric(expr[gene, ]))[, 1]        # cohort-internal z-score

  fit_stratum <- function(dd, s) {
    ev <- sum(dd$event == 1, na.rm = TRUE)
    if (!is.finite(ev) || ev < s$min_events) return(NULL)
    safe_cox_screen(stats::as.formula(paste("Surv(time, event) ~", s$rhs)), dd)
  }

  overall <- fit_stratum(d,                                    spec$overall)
  mut_res <- fit_stratum(dplyr::filter(d, idh == "Mutant"),    spec$mutant)
  wt_res  <- fit_stratum(dplyr::filter(d, idh == "WT"),        spec$wildtype)

  pick <- function(r, f) if (is.null(r)) NA_real_ else as.numeric(r[[f]])
  data.frame(
    gene = gene,
    overall_HR = pick(overall, "HR"), overall_lower = pick(overall, "lower95"),
    overall_upper = pick(overall, "upper95"), overall_p = pick(overall, "p"),
    mut_HR = pick(mut_res, "HR"), mut_lower = pick(mut_res, "lower95"),
    mut_upper = pick(mut_res, "upper95"), mut_p = pick(mut_res, "p"),
    wt_HR = pick(wt_res, "HR"), wt_lower = pick(wt_res, "lower95"),
    wt_upper = pick(wt_res, "upper95"), wt_p = pick(wt_res, "p")
  )
}

# ---- run one cohort ----------------------------------------------------------
run_screen <- function(cohort_name, out_path) {

  spec <- SCREEN_SPEC[[cohort_name]]
  if (is.null(spec)) stop("no screen specification for cohort ", cohort_name)

  co    <- load_cohort(cohort_name)
  expr  <- co$expr
  clin  <- co$clin

  genes <- rownames(expr)
  genes <- genes[str_detect(genes, GENE_FAMILY)]
  message(sprintf("[%s] %d gene symbols matching '%s', %d samples, %d events",
                  cohort_name, length(genes), GENE_FAMILY, ncol(expr),
                  sum(clin$event == 1, na.rm = TRUE)))
  for (nm in names(spec))
    message(sprintf("[%s]   %-8s ~ %s   (>= %d events)",
                    cohort_name, nm, spec[[nm]]$rhs, spec[[nm]]$min_events))

  res <- map_dfr(genes, analyze_gene, expr = expr, clin = clin, spec = spec) %>%
    mutate(FDR_overall = p.adjust(overall_p, "BH"),
           FDR_mutant  = p.adjust(mut_p,     "BH"),
           FDR_wt      = p.adjust(wt_p,      "BH"))

  write.csv(res, out_path, row.names = FALSE)
  message(sprintf("[%s] written: %s", cohort_name, out_path))

  # ---- reproducibility report: fitted counts and the gene's rank per stratum
  for (pre in c("overall", "mut", "wt")) {
    ok <- res[is.finite(res[[paste0(pre, "_p")]]) &
              is.finite(res[[paste0(pre, "_HR")]]), ]
    ok <- ok[order(ok[[paste0(pre, "_p")]]), ]
    if (!nrow(ok)) {
      message(sprintf("[%s] %-8s no gene fitted (event threshold not reached)",
                      cohort_name, pre))
      next
    }
    i <- which(ok$gene == GENE_OF_INTEREST)
    if (length(i))
      message(sprintf("[%s] %-8s %s ranks %d of %d fitted: HR %.3f (%.3f-%.3f), p %.3g",
                      cohort_name, pre, GENE_OF_INTEREST, i, nrow(ok),
                      ok[[paste0(pre, "_HR")]][i], ok[[paste0(pre, "_lower")]][i],
                      ok[[paste0(pre, "_upper")]][i], ok[[paste0(pre, "_p")]][i]))
    else
      message(sprintf("[%s] %-8s %d genes fitted; %s not among them",
                      cohort_name, pre, nrow(ok), GENE_OF_INTEREST))
  }
  invisible(res)
}

run_screen("TCGA_LGG", PATH_LGG_COX_ALL)
run_screen("TCGA_GBM", PATH_GBM_COX_ALL)
