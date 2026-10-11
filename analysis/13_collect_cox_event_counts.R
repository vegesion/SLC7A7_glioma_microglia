# ============================================================================
# 13_collect_cox_event_counts.R
# ----------------------------------------------------------------------------
# PURPOSE  Recover n and the event count for every Cox model without redrawing
#          figures.
#
# INPUT    fitted cohorts
# OUTPUT   event-count table, feeding Additional file 4
#
# Sources R/setup.R, which loads config/config.R. See docs/PIPELINE.md for
# the execution order and docs/OUTPUT_MAP.md for the mapping to additional
# files and figures. Comments below this block are the working notes kept
# from development and are in Korean.
# ============================================================================

# =============================================================================
# 13_collect_cox_event_counts.R
# 원고에 보고해야 하는 모든 Cox 모형의 n / event 수를 한 번에 수집
# -----------------------------------------------------------------------------
# 목적: 결과 CSV 들이 nevent 을 저장하지 않아 §7(통계 보고 규약)의
#       "Cox 결과에는 n 과 event 수를 함께 쓴다" 를 만족시킬 수 없음.
#       피규어를 다시 그리지 않고 적합만 반복해 숫자만 회수한다.
#
# 출력: results/cox_event_counts.csv
#       열: block, cohort, stratum, model, term, n, nevent, covars
#
# 주의: 이 스크립트는 HR/CI/p 를 다시 보고하지 않는다. 기존 CSV 의 값을
#       덮어쓰지 않기 위해 n/nevent 만 뽑는다. HR 이 기존 CSV 와 다르게 나오면
#       그 자체가 보고해야 할 불일치이므로, 아래 CHECK 블록이 경고를 띄운다.
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({ library(survival); library(dplyr) })

GENE <- GENE_OF_INTEREST
out  <- list()

add <- function(block, cohort, stratum, model, term, fit, covars) {
  data.frame(block = block, cohort = cohort, stratum = stratum, model = model,
             term = term, n = fit$n, nevent = fit$nevent,
             covars = paste(covars, collapse = " + "), stringsAsFactors = FALSE)
}

cohorts <- load_all_cohorts(); stopifnot(length(cohorts) > 0)

## =============================================================================
## A. Fig 3 본문 다변량 Cox  (04_multicohort_validation.R:31 과 동일 호출)
##    covars 는 multivariable_cox 의 기본값. 코호트에 없는 변수는
##    signature.R:63 의 keep 필터가 자동으로 떨어뜨리므로 실제 사용된
##    공변량을 fit 에서 다시 읽어 기록한다.
## =============================================================================
DEFAULT_COVARS <- c("age", "grade", "idh", "mgmt")

for (nm in names(cohorts)) {
  co <- cohorts[[nm]]
  if (!GENE %in% rownames(co$expr)) next
  mv <- multivariable_cox(co$expr, co$clin, GENE, covars = DEFAULT_COVARS)
  if (is.null(mv)) { message("[A] 표본 부족으로 건너뜀: ", nm); next }
  used <- setdiff(all.vars(formula(mv$fit))[-(1:2)], c("time", "event"))
  out[[length(out) + 1]] <- add("Fig3_multivariable", nm, "all", "cox", "gene",
                                mv$fit, used)
  # CHECK: 기존 CSV 와 HR 이 일치하는지
  r <- mv$table[mv$table$term == "gene", ]
  message(sprintf("[A] %-10s n=%d events=%d  HR=%.4f (%.4f-%.4f) p=%.3g",
                  nm, mv$fit$n, mv$fit$nevent, r$HR, r$lower, r$upper, r$p))
}

## =============================================================================
## B. Fig 3A / Fig 3D  IDH 층화
##    stratified_cox(by="idh") 는 covars = setdiff(c("age","grade","mgmt"), "idh")
##    를 넘긴다 (signature.R:81). 층별 fit 을 직접 얻기 위해 같은 인자로
##    multivariable_cox 를 호출한다.
## =============================================================================
STRAT_COVARS <- c("age", "grade", "mgmt")

for (nm in names(cohorts)) {
  co <- cohorts[[nm]]
  if (!GENE %in% rownames(co$expr)) next
  for (l in stats::na.omit(unique(co$clin$idh))) {
    idx <- which(co$clin$idh == l)
    mv <- multivariable_cox(co$expr[, idx, drop = FALSE],
                            co$clin[idx, , drop = FALSE], GENE,
                            covars = STRAT_COVARS)
    if (is.null(mv)) { message("[B] 표본 부족: ", nm, " idh=", l); next }
    used <- setdiff(all.vars(formula(mv$fit))[-(1:2)], c("time", "event"))
    out[[length(out) + 1]] <- add("Fig3_stratified", nm, paste0("idh=", l),
                                  "cox", "gene", mv$fit, used)
    r <- mv$table[mv$table$term == "gene", ]
    message(sprintf("[B] %-10s idh=%-7s n=%d events=%d  HR=%.4f p=%.3g",
                    nm, l, mv$fit$n, mv$fit$nevent, r$HR, r$p))
  }
}

## =============================================================================
## C. Fig 4C head-to-head  (Fig4CD_specificity_panels.R:42-57 과 동일 모형)
##    노출변수만 바뀌고 공변량 결측 패턴은 같으므로 코호트당 n/events 는
##    마커 전체에서 동일해야 한다. 그렇지 않으면 그 마커에 결측이 있다는
##    뜻이므로 마커별로 모두 기록한다.
## =============================================================================
H2H_COHORTS <- c("TCGA_LGG", "CGGA_693", "CGGA_325")
MARKERS     <- c("SLC7A7", "PTPRC", "C1QB", "C1QA", "TYROBP",
                 "ITGAM", "AIF1", "CSF1R", "CD68")

for (nm in intersect(H2H_COHORTS, names(cohorts))) {
  co <- cohorts[[nm]]; cl <- co$clin
  for (g in intersect(MARKERS, rownames(co$expr))) {
    d <- data.frame(time = cl$time, event = cl$event,
                    x = as.numeric(scale(as.numeric(co$expr[g, ]))),
                    age = cl$age, idh = cl$idh)
    d <- d[stats::complete.cases(d) & d$time > 0, , drop = FALSE]
    if (nrow(d) < 30 || length(unique(d$idh)) < 2 || sd(d$x) == 0) next
    fit <- coxph(Surv(time, event) ~ x + age + idh, data = d)
    out[[length(out) + 1]] <- add("Fig4C_head_to_head", nm, "all", "cox", g,
                                  fit, c("age", "idh"))
  }
  rng <- do.call(rbind, out)
  rng <- rng[rng$block == "Fig4C_head_to_head" & rng$cohort == nm, ]
  message(sprintf("[C] %-10s n=%s events=%s (마커 %d종)",
                  nm, paste(unique(rng$n), collapse = "/"),
                  paste(unique(rng$nevent), collapse = "/"), nrow(rng)))
}

## =============================================================================
## D. Fig 4D 침윤 보정 민감도  (R/deconvolution.R:45 adjust_for_fraction 과 동일)
##    ESTIMATE / xCell 재실행을 피하기 위해 05 가 저장한
##    Fig4_deconvolution_scores.csv 를 읽는다.
## =============================================================================
DECON_FILE <- file.path(DIR_RESULTS, "Fig4_deconvolution_scores.csv")
MYE_COLS   <- c("Macrophages", "Macrophages M1", "Macrophages M2", "Monocytes",
                "ImmuneScore.x", "MicroenvironmentScore")
D_COVARS   <- c("age", "grade", "idh")

if (!file.exists(DECON_FILE)) {
  message("[D] 건너뜀: ", DECON_FILE, " 없음 (05 의 (c) 블록을 먼저 실행)")
} else {
  frac <- read.csv(DECON_FILE, check.names = FALSE, stringsAsFactors = FALSE)
  co   <- cohorts[[TRAIN_COHORT]]
  cl   <- co$clin
  v    <- as.numeric(scale(as.numeric(co$expr[GENE, ])))
  for (cc in intersect(MYE_COLS, names(frac))) {
    fr <- frac[match(colnames(co$expr), frac$sample), cc]
    d  <- data.frame(time = cl$time, event = cl$event, gene = v,
                     myeloid = as.numeric(scale(suppressWarnings(as.numeric(fr)))),
                     cl[, intersect(D_COVARS, names(cl)), drop = FALSE])
    keep <- vapply(d, function(x) length(unique(stats::na.omit(x))) > 1, logical(1))
    d <- d[, keep, drop = FALSE]; d <- d[stats::complete.cases(d) & d$time > 0, ]
    if (nrow(d) < 30) { message("[D] 표본 부족: ", cc); next }
    rhs <- setdiff(names(d), c("time", "event"))
    fit <- coxph(stats::as.formula(paste("Surv(time, event) ~",
                                         paste(rhs, collapse = " + "))), data = d)
    out[[length(out) + 1]] <- add("Fig4D_infiltration", TRAIN_COHORT, "all",
                                  "cox", paste0("gene | adj=", cc), fit,
                                  c(D_COVARS, cc))
    message(sprintf("[D] adj=%-22s n=%d events=%d", cc, fit$n, fit$nevent))
  }
}

## =============================================================================
## E. Fig 2A / 2B  pan-SLC 457개 스크린
##    ★ 457회 반복이 불필요하다. analyze_gene_cox 는 유전자마다 노출변수만
##      바꾸고 공변량 집합과 발현 완전성은 동일하므로, 한 subgroup 안의
##      분석 n 과 event 수는 457개 유전자 전부 같다. subgroup 당 1회 적합으로
##      충분하다 (LGG 3군 + GBM 3군 = 6회).
##    ★ analyze_gene_cox:74 의 events 는 coxph 가 공변량 결측 행을 떨어뜨리기
##      전의 수이므로 그 값을 보고하면 안 된다 (상한값).
##    ★ 이 블록은 screen 이 쓴 원본 임상표(paper_* 열) 를 요구한다.
##      load_cohort 의 harmonise_clinical 을 거친 표가 아니다.
## =============================================================================
SCREEN_SPECS <- list(LGG = list(tag = "LGG", spec = cohort_spec_LGG()),
                     GBM = list(tag = "GBM", spec = cohort_spec_GBM()))

for (sp in SCREEN_SPECS) {
  rds <- file.path(DIR_DATA_PROC, paste0(sp$tag, "_counts_clinical.rds"))
  nrm <- if (sp$tag == "LGG") PATH_LGG_NORM else PATH_GBM_NORM
  if (!file.exists(rds) || !file.exists(nrm)) {
    message("[E] 건너뜀 (", sp$tag, "): 원본 입력 없음 — ", rds, " / ", nrm)
    next
  }
  expr_df  <- load_expr_matrix(nrm)
  clin_raw <- readRDS(rds)$clinical
  merged   <- build_merged_df(expr_df[expr_df$Gene == GENE, , drop = FALSE], clin_raw)
  if (isTRUE(COX_SCALE_EXPR)) merged$logCPM_scaled <- scale(merged$logCPM)[, 1]
  else                        merged$logCPM_scaled <- merged$logCPM

  for (sgn in names(sp$spec)) {
    sg <- sp$spec[[sgn]]
    d  <- if (is.null(sg$filter)) merged else merged[which(eval(sg$filter, envir = merged)), ]
    f  <- stats::as.formula(paste0(
      "Surv(paper_Survival..months., paper_Vital.status..1.dead.) ~ logCPM_scaled + ",
      sg$covars))
    fit <- tryCatch(coxph(f, data = d), error = function(e) NULL)
    if (is.null(fit)) { message("[E] 적합 실패: ", sp$tag, "/", sgn); next }
    out[[length(out) + 1]] <- add("Fig2_panSLC_screen", paste0("TCGA_", sp$tag),
                                  sgn, "cox", "logCPM_scaled", fit,
                                  strsplit(sg$covars, " \\+ ")[[1]])
    gate <- sum(d$paper_Vital.status..1.dead. == 1, na.rm = TRUE)
    message(sprintf("[E] %s/%-7s n=%d events=%d   (min_events 관문값 %d — 보고용 아님)",
                    sp$tag, sgn, fit$n, fit$nevent, gate))
  }
}

## =============================================================================
res <- dplyr::bind_rows(out)
write_result(res, "cox_event_counts.csv")
print(as.data.frame(res))
message("\n완료: results/cox_event_counts.csv")
message("Firth 적합(coxphf)으로 넘어간 행은 nevent 가 없으므로 이 표에 나타나지 않는다.")
