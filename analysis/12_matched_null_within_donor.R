# ============================================================================
# 12_matched_null_within_donor.R
# ----------------------------------------------------------------------------
# PURPOSE  The matched-control gene null built for the within-donor contrast.
#
# INPUT    cell-level matrix and the donor index
# OUTPUT   Supp_wdnull_*, feeding Additional file 10
#
# Sources R/setup.R, which loads config/config.R. See docs/PIPELINE.md for
# the execution order and docs/OUTPUT_MAP.md for the mapping to additional
# files and figures. Comments below this block are the working notes kept
# from development and are in Korean.
# ============================================================================

# =============================================================================
# 12_matched_null_within_donor.R
# within-donor 대비에 대한 유전자 매칭 empirical null
# 전제: 08b 의 §1-§2 실행 후 (meta, expr_t, genes 존재)
# =============================================================================
suppressPackageStartupMessages({ library(Matrix); library(fgsea); library(dplyr) })

CALIPER_Q       <- 0.05
CALIPER_Q_TIGHT <- 0.02
EXCLUDE_COR_ABS <- 0.50
N_SUB <- 450; N_DRAW <- 20
TARGET_PATHWAYS <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_GLYCOLYSIS","HALLMARK_G2M_CHECKPOINT","HALLMARK_HYPOXIA",
  "HALLMARK_MITOTIC_SPINDLE","HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE")

## ── 사전 계산 (1회) ─────────────────────────────────────────────────────────
donor <- factor(meta$donor)
D  <- sparse.model.matrix(~ 0 + donor)        # cells x donors
nd <- Matrix::colSums(D)
Y  <- expr_t                                   # cells x genes (dgCMatrix)
N  <- nrow(Y); G <- ncol(Y); K <- length(nd)
DF <- N - K - 2

demean <- function(v) v - as.numeric(D %*% (Matrix::crossprod(D, v) / nd))

yy  <- Matrix::colSums(Y * Y)
S   <- as.matrix(Matrix::crossprod(D, Y))      # donors x genes
yyw <- yy - Matrix::colSums(S * S / nd)        # ỹᵀỹ, 유전자별
rm(S); gc(verbose = FALSE)

b  <- demean(meta$depth)                       # d_within (원척도)
bb <- sum(b * b)
bY <- as.numeric(Matrix::crossprod(Y, b))

#' 예측변수 x(세포 벡터)에 대한 전 유전자 t_within
t_within_fast <- function(x) {
  a  <- demean(x)
  aa <- sum(a * a); if (!is.finite(aa) || aa == 0) return(NULL)
  ab <- sum(a * b)
  det <- aa * bb - ab * ab
  aY <- as.numeric(Matrix::crossprod(Y, a))
  beta_a <- ( bb * aY - ab * bY) / det
  beta_b <- (-ab * aY + aa * bY) / det
  sse <- yyw - (beta_a * aY + beta_b * bY)
  sse[sse < .Machine$double.eps] <- NA_real_
  se_a <- sqrt((sse / DF) * (bb / det))
  setNames(beta_a / se_a, colnames(Y))
}

## 08b 의 run_gsea_pair 와 동일한 설정. HALLMARK 는 08b 가 만든 것을 그대로 쓴다.
N_PERM_SIMPLE <- 100000          # 08b 는 100000. null 전체에는 과하므로 낮추되 §0에서 검증한다.

nes_of <- function(tv, seed) {
  if (is.null(tv)) return(setNames(rep(NA_real_, length(TARGET_PATHWAYS)), TARGET_PATHWAYS))
  tv <- sort(tv[is.finite(tv)], decreasing = TRUE)     # 08b 와 동일하게 정렬
  set.seed(seed)
  r <- suppressWarnings(fgsea::fgsea(
    pathways = HALLMARK, stats = tv,
    minSize = 10, maxSize = 500,                  # ← 15 가 아니라 10
    eps = 0,                                      # ← 핵심. 없으면 극단값이 NA 가 된다
    nPermSimple = N_PERM_SIMPLE, nproc = 1))
  out <- setNames(rep(NA_real_, length(TARGET_PATHWAYS)), TARGET_PATHWAYS)
  h <- match(TARGET_PATHWAYS, r$pathway); out[!is.na(h)] <- r$NES[h[!is.na(h)]]
  out
}

## ── 단계 0: 근사 검증 ───────────────────────────────────────────────────────
t_fast <- t_within_fast(meta$g)
ref <- read.csv(file.path(DIR_RESULTS, "Fig6_wb_gene_stats.csv"), stringsAsFactors = FALSE)
cm  <- intersect(names(t_fast), ref$gene)
rho <- cor(t_fast[cm], ref$t_within[match(cm, ref$gene)], method = "spearman")
message(sprintf("[검증] 공통 %d개 유전자, Spearman rho = %.4f", length(cm), rho))

## 관측 NES 가 원고 보고값(08b 산출물)을 재현하는지 — 이게 진짜 검증이다
obs_nes <- nes_of(t_fast, GSEA_SEED)
pub <- read.csv(file.path(DIR_RESULTS, "Fig6_wb_gsea_within.csv"), stringsAsFactors = FALSE)
chk <- data.frame(pathway = TARGET_PATHWAYS,
                  NES_null_framework = round(obs_nes, 3),
                  NES_published      = round(pub$NES[match(TARGET_PATHWAYS, pub$pathway)], 3))
chk$diff <- round(chk$NES_null_framework - chk$NES_published, 3)
print(chk)
write.csv(chk, file.path(DIR_RESULTS, "Supp_wdnull_observed_check.csv"), row.names = FALSE)

if (any(is.na(obs_nes)))
  stop("아직 NA 가 있습니다. eps = 0 과 minSize = 10 이 적용됐는지 확인하세요.")
if (max(abs(chk$diff), na.rm = TRUE) > 0.05)
  warning("관측 NES 가 원고값과 0.05 이상 차이납니다. N_PERM_SIMPLE 을 100000 으로 올려 다시 확인하세요.")

## ── 단계 1: 세포 수준 매칭 ──────────────────────────────────────────────────
gw <- demean(meta$g); gww <- sum(gw * gw)
gstat_c <- data.frame(
  gene      = colnames(Y),
  det_cell  = Matrix::colMeans(Y > 0),
  mean_cell = Matrix::colMeans(Y),
  sd_within = sqrt(yyw / (N - K)),
  rho_w     = as.numeric(Matrix::crossprod(Y, gw)) / sqrt(gww * yyw),
  stringsAsFactors = FALSE)
gstat_c <- gstat_c[is.finite(gstat_c$sd_within) & gstat_c$sd_within > 0, ]

ref_c <- c(det  = mean(meta$g > 0),
           mean = mean(meta$g),
           sdw  = sd(gw))
z  <- scale(gstat_c[, c("mean_cell", "sd_within", "det_cell")])
z0 <- c((ref_c["mean"] - attr(z,"scaled:center")[1]) / attr(z,"scaled:scale")[1],
        (ref_c["sdw"]  - attr(z,"scaled:center")[2]) / attr(z,"scaled:scale")[2],
        (ref_c["det"]  - attr(z,"scaled:center")[3]) / attr(z,"scaled:scale")[3])
gstat_c$dist <- sqrt(rowSums((z - matrix(z0, nrow(z), 3, byrow = TRUE))^2))

base_all <- gstat_c %>% filter(gene != GENE_OF_INTEREST, gene != "SLC3A2")

## ── A. caliper 분위수별 풀 크기 (세포 수준 매칭 기준) ───────────────────────
## §1 의 gstat_c, base_dec, base_all 생성까지 실행한 뒤 이 블록을 돌린다.
cat(sprintf("eligible: dec %d / all %d (총 universe %d)\n",
            nrow(base_dec), nrow(base_all), nrow(gstat_c)))
for (q in c(0.04, 0.06, 0.08, 0.10, 0.15, 0.20, 0.30)) {
  s <- base_dec[base_dec$dist <= quantile(base_dec$dist, q), ]
  cat(sprintf("  q=%.2f: n=%4d | mean %.2f~%.2f | det %.3f~%.3f | sd_w %.3f~%.3f\n",
              q, nrow(s), min(s$mean_cell), max(s$mean_cell),
              min(s$det_cell), max(s$det_cell),
              min(s$sd_within), max(s$sd_within)))
}
cat(sprintf("SLC7A7: mean %.3f | det %.3f | sd_w %.3f\n",
            mean(meta$g), mean(meta$g > 0), sd(demean(meta$g))))

## ── B. 대조 유전자 1개 처리 시간 ────────────────────────────────────────────
g1 <- base_dec$gene[which.min(base_dec$dist)]
print(system.time({
  tv <- t_within_fast(as.numeric(Y[, g1]))
  nes_of(tv, GSEA_SEED + 1)
}))


base_dec <- base_all %>% filter(abs(rho_w) < EXCLUDE_COR_ABS)
cal <- function(d, q) d %>% filter(dist <= quantile(dist, q, na.rm = TRUE))
sets <- list(dec_main  = cal(base_dec, CALIPER_Q)$gene,
             dec_tight = cal(base_dec, CALIPER_Q_TIGHT)$gene,
             all_main  = cal(base_all, CALIPER_Q)$gene,
             all_tight = cal(base_all, CALIPER_Q_TIGHT)$gene)
fit_genes <- union(sets$dec_main, sets$all_main)
message(sprintf("[매칭] dec_main %d | all_main %d | 적합 대상 %d",
                length(sets$dec_main), length(sets$all_main), length(fit_genes)))
write.csv(gstat_c[gstat_c$gene %in% fit_genes, ],
          file.path(DIR_RESULTS, "Supp_wdnull_matched_genes.csv"), row.names = FALSE)

## ── 단계 2: null 실행 ───────────────────────────────────────────────────────
null_mat <- matrix(NA_real_, length(fit_genes), length(TARGET_PATHWAYS),
                   dimnames = list(fit_genes, TARGET_PATHWAYS))
pb <- txtProgressBar(0, length(fit_genes), style = 3)
for (i in seq_along(fit_genes)) {
  xi <- as.numeric(Y[, fit_genes[i]])
  null_mat[i, ] <- nes_of(t_within_fast(xi), GSEA_SEED + i)
  setTxtProgressBar(pb, i)
}
close(pb)

emp_for <- function(gs, lab) {
  idx <- which(rownames(null_mat) %in% gs)
  bind_rows(lapply(TARGET_PATHWAYS, function(p) {
    nl <- null_mat[idx, p]; nl <- nl[is.finite(nl)]; ob <- obs_nes[[p]]
    data.frame(null_set = lab, pathway = p, observed_NES = ob, n_null = length(nl),
               null_mean = mean(nl), null_sd = sd(nl),
               z_vs_null = (ob - mean(nl)) / sd(nl),
               p_emp_two = (1 + sum(abs(nl) >= abs(ob))) / (1 + length(nl)),
               p_emp_one = (1 + sum(if (ob >= 0) nl >= ob else nl <= ob)) / (1 + length(nl)))
  }))
}
emp <- bind_rows(lapply(names(sets), function(s) emp_for(sets[[s]], s))) %>%
  group_by(null_set) %>% mutate(p_emp_two_BH = p.adjust(p_emp_two, "BH")) %>% ungroup()
write.csv(emp, file.path(DIR_RESULTS, "Supp_wdnull_empirical_p.csv"), row.names = FALSE)
print(as.data.frame(emp))

stab <- bind_rows(lapply(names(sets), function(s) {
  ia <- which(rownames(null_mat) %in% sets[[s]]); if (length(ia) <= N_SUB) return(NULL)
  set.seed(GSEA_SEED + 10000L)
  dr <- replicate(N_DRAW, sample(ia, N_SUB), simplify = FALSE)
  bind_rows(lapply(TARGET_PATHWAYS, function(p) {
    ob <- obs_nes[[p]]
    ps <- vapply(dr, function(ix) { nl <- null_mat[ix, p]; nl <- nl[is.finite(nl)]
    (1 + sum(abs(nl) >= abs(ob))) / (1 + length(nl)) }, numeric(1))
    data.frame(null_set = s, pathway = p, p_emp_min = min(ps), p_emp_median = median(ps),
               p_emp_max = max(ps), n_draw_sig_0.05 = sum(ps < 0.05))
  }))
}))
write.csv(stab, file.path(DIR_RESULTS, "Supp_wdnull_draw_stability.csv"), row.names = FALSE)

nd_long <- as.data.frame(null_mat) %>% tibble::rownames_to_column("control_gene") %>%
  tidyr::pivot_longer(-control_gene, names_to = "pathway", values_to = "NES")
write.csv(nd_long, file.path(DIR_RESULTS, "Supp_wdnull_NES_distribution.csv"), row.names = FALSE)


summ <- bind_rows(lapply(names(sets), function(s) {
  ia <- which(rownames(null_mat) %in% sets[[s]])
  bind_rows(lapply(TARGET_PATHWAYS, function(p) {
    ob <- obs_nes[[p]]; nl <- null_mat[ia, p]; nl <- nl[is.finite(nl)]
    data.frame(null_set = s, pathway = p, n_null = length(nl), observed_NES = ob,
               null_median = median(nl),
               pct_same_sign = mean(sign(nl) == sign(ob)) * 100,
               pct_abs_ge_obs = mean(abs(nl) >= abs(ob)) * 100,
               obs_percentile = mean(nl < ob) * 100)
  }))
}))
write.csv(summ, file.path(DIR_RESULTS, "Supp_wdnull_summary.csv"), row.names = FALSE)
all_nl <- null_mat[which(rownames(null_mat) %in% sets$dec_main), TARGET_PATHWAYS]
all_nl <- all_nl[is.finite(all_nl)]
message(sprintf("[전체] n=%d | |NES|>2 %.1f%% | |NES|>3 %.1f%% | NA %d",
                length(all_nl), mean(abs(all_nl) > 2) * 100, mean(abs(all_nl) > 3) * 100,
                length(sets$dec_main) * length(TARGET_PATHWAYS) - length(all_nl)))
message("done.")


