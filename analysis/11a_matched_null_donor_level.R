# ============================================================================
# 11a_matched_null_donor_level.R
# ----------------------------------------------------------------------------
# PURPOSE  Matched-control gene null for the donor-level contrast, at two
#          caliper widths and with and without excluding co-regulated genes.
#
# INPUT    donor pseudobulk
# OUTPUT   Supp_null_*, feeding Additional files 5 and 10
#
# Sources R/setup.R, which loads config/config.R. See docs/PIPELINE.md for
# the execution order and docs/OUTPUT_MAP.md for the mapping to additional
# files and figures. Comments below this block are the working notes kept
# from development and are in Korean.
# ============================================================================

# =============================================================================
# 11a_matched_null_donor_level.R  —  음성대조(negative control) 검증  [전면 개정판]
# -----------------------------------------------------------------------------
# 목적 / Purpose
#   채택한 "donor 수준 연속형" 분석이 sequencing depth 인공산물이 아님을,
#   SLC7A7과 발현수준·검출률·분산이 매칭된 대조 유전자 200개로 만든
#   empirical null 분포로 증명한다.
#
#   핵심 논리: SLC7A7을 예측변수로 넣었을 때 나온 OXPHOS NES(+2.053)가,
#   "SLC7A7과 통계적 성질이 같지만 생물학적으로 무관한 유전자 200개"를
#   같은 자리에 넣었을 때의 NES 분포에서 얼마나 극단적인가?
#
# 입력 / Input
#   data/processed/myeloid.rds   (Seurat v5, myeloid 세포. 06에서 생성)
#   config/config.R, R/setup.R
#
# 출력 / Output
#   results/Supp_null_matched_genes.csv        선정된 대조 유전자 200개 + 매칭 지표
#   results/Supp_null_NES_distribution.csv     200회 × pathway NES 원자료
#   results/Supp_null_empirical_p.csv          ★ 논문 Supplementary Table
#   results/Supp_null_SLC3A2_partialcor.csv    SLC3A2 공발현의 매칭 대조
#   results/Supp_depth_sensitivity_legacy.csv  (구) 세포단위 이분 대조 — 폐기근거용
#   figures/SuppX_null_NES_hist.pdf            ★ 논문 Supplementary Figure
#
# 원본 대체: 기존 11_depth_sensitivity.R(depth비 · glmer OR · 매칭 5개)를
#            §5 legacy 블록으로 축소 보존하고, §2–4 를 신규 추가.
# =============================================================================

source(here::here("R", "setup.R"))

suppressPackageStartupMessages({
  library(Matrix); library(edgeR); library(limma); library(fgsea)
  library(msigdbr); library(ggplot2); library(dplyr); library(tidyr); library(tibble)
})

# ── ★ 분석 구획 (compartment) ───────────────────────────────────────────────
#  음성대조는 "본문에서 실제로 쓰는 집단"에서 돌려야 합니다.
#  microglia 논문이면 COMPARTMENT <- "microglia" 로 두고,
#  MDM 은 특이성 대조로 한 번 더 돌리세요 (같은 스크립트, 이 줄만 변경).
#
#    "microglia" : TAM-MG 계열만            (n=77 donor, median 211 cells/donor)
#    "mac"       : TAM-BDM + Mono 계열만    (n=77 donor, median 287)  ← 특이성 대조
#    "myeloid"   : 전체                      (n=91 donor, median 469)
COMPARTMENT <- "microglia"

# 구획 정의 — obj$annotation_level_4 라벨 기준 (필요시 수정)
COMPARTMENT_LABELS <- list(
  microglia = c("TAM-MG aging sig","TAM-MG pro-infl I","TAM-MG pro-infl II",
                "TAM-MG prolif"),
  mac       = c("TAM-BDM INF","TAM-BDM MHC","TAM-BDM anti-infl",
                "TAM-BDM hypoxia/MES","Mono naive","Mono anti-infl","Mono hypoxia"),
  myeloid   = NULL   # NULL = 필터 없음
)
LABEL_COL <- "annotation_level_4"

# ── 파라미터 ────────────────────────────────────────────────────────────────
# ── matched-control null: caliper 매칭 ──────────────────────────────────────
# 표준화 (mean_pb, sd_don, det) 공간에서 SLC7A7 까지의 거리가 CALIPER_Q 분위수
# 이내인 유전자를 모두 대조로 쓴다. CALIPER_Q_TIGHT 는 그 부분집합이므로
# 추가 적합 없이 "매칭을 더 좁혔을 때"의 결과가 함께 나온다.
CALIPER_Q        <- 0.05     # 주 분석 caliper (진단: n≈893, det 0.235-0.576)
CALIPER_Q_TIGHT  <- 0.02     # 중첩 민감도    (진단: n≈357, det 0.309-0.500)
N_SUB            <- 450      # 안정성 확인용 부분추출 크기 (주 caliper 의 약 절반)
N_DRAW           <- 20       # 부분추출 반복 횟수 (저장된 NES 재사용 — 추가 적합 없음)
POOL_CAP         <- 1800     # 실제 적합할 유전자 수 상한 (런타임 보호)
RUN_BOTH_RHO     <- TRUE     # TRUE: 공조절 제외 null 과 미제외 null 을 함께 산출
MATCH_K_POOL     <- 3000     # (미사용 — 거리순 상위 선택을 쓰지 않는다)
EXCLUDE_COR_ABS  <- 0.50     # |Spearman(donor)| 이 이 값 이상이면 제외 (공조절 유전자 배제)
DET_SUBSAMPLE    <- 50000    # 검출률 계산용 세포 subsample
MIN_DONOR_CELLS  <- 50       # ★ microglia 는 donor당 세포수가 적음(median 211).
#   20(config 기본)은 pseudobulk 가 너무 시끄러움.
#   50 / 100 두 값으로 민감도 확인 권장:
#     20 → 77 donor / 50 → 69 / 100 → 54
MIN_DATASET_DONOR<- 2        # dataset 더미가 자유도를 다 먹는 것 방지
NPROC            <- 1        # fgsea 병렬 (Windows 는 1 권장)

TARGET_PATHWAYS <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",   # ← 주 검정 대상
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_HYPOXIA",
  "HALLMARK_MITOTIC_SPINDLE",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",      # ← null 이어야 하는 음성 기준
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE"
)

GENE <- GENE_OF_INTEREST          # "SLC7A7"
PARTNER <- GENE_PARTNER           # "SLC3A2"

# =============================================================================
# §1. donor 수준 행렬 구축 (관측 분석과 완전히 동일한 전처리)
# =============================================================================

message("[1] loading microglia object ...")
obj <- readRDS(file.path(DIR_DATA_PROC, "myeloid.rds"))

# ★ 구획 서브셋 — 본문 분석과 동일한 집단으로 좁힌다
labs <- COMPARTMENT_LABELS[[COMPARTMENT]]
if (!is.null(labs)) {
  stopifnot(LABEL_COL %in% colnames(obj@meta.data))
  missing_lab <- setdiff(labs, unique(obj@meta.data[[LABEL_COL]]))
  if (length(missing_lab))
    warning("라벨 불일치: ", paste(missing_lab, collapse = ", "))
  obj <- obj[, obj@meta.data[[LABEL_COL]] %in% labs]
}
message(sprintf("    COMPARTMENT = %s | cells = %d", COMPARTMENT, ncol(obj)))

# 출력 파일 접미사 (구획별로 덮어쓰지 않게)
SFX <- paste0("_", COMPARTMENT)

md <- obj@meta.data
stopifnot(any(c("donor_id", "dataset") %in% colnames(md)) ||
            any(c("donor", "dataset", "author") %in% colnames(md)))
donor_col <- if ("donor_id" %in% colnames(md)) "donor_id" else "donor"

# depth 공변량: config 의 DEPTH_COVARIATE (기본 nCount_RNA)
stopifnot(DEPTH_COVARIATE %in% colnames(md))

cnt <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
dat <- SeuratObject::LayerData(obj, assay = "RNA", layer = "data")   # LogNormalize

donor_f <- factor(as.character(md[[donor_col]]))
keep_donor <- names(which(table(donor_f) >= MIN_DONOR_CELLS))
cell_ok <- donor_f %in% keep_donor

cnt <- cnt[, cell_ok, drop = FALSE]
dat <- dat[, cell_ok, drop = FALSE]
md  <- md[cell_ok, , drop = FALSE]

donor_f <- droplevels(donor_f[cell_ok])

# 희소 지시행렬 (cells × donors) — AggregateExpression 버전 버그 회피
ind <- Matrix::sparse.model.matrix(~ 0 + donor_f)
colnames(ind) <- levels(donor_f)
n_cell <- Matrix::colSums(ind)

message(sprintf("    donors = %d, cells = %d", ncol(ind), nrow(md)))

## (a) pseudobulk counts (genes × donors) — DE 입력
pb <- as.matrix(cnt %*% ind)
colnames(pb) <- colnames(ind)

## (b) donor 평균 정규화 발현 (genes × donors) — 예측변수 후보
dm <- sweep(as.matrix(dat %*% ind), 2, n_cell, "/")

## (c) donor 메타 (depth · dataset · n)
donor_meta <- md %>%
  dplyr::mutate(donor_id = as.character(.data[[donor_col]])) %>%
  dplyr::group_by(donor_id ) %>%
  dplyr::summarise(depth   = mean(.data[[DEPTH_COVARIATE]]),
                   dataset = names(which.max(table(author))),
                   n_cell  = dplyr::n(), .groups = "drop") %>%
  as.data.frame()

rownames(donor_meta) <- donor_meta$donor_id 
donor_meta <- donor_meta[colnames(pb), ]

# dataset 이 donor 1명뿐이면 제외 (더미가 자유도를 소진)
ds_ok <- names(which(table(donor_meta$dataset) >= MIN_DATASET_DONOR))
keep2  <- donor_meta$dataset %in% ds_ok
pb <- pb[, keep2, drop = FALSE]; dm <- dm[, keep2, drop = FALSE]
donor_meta <- donor_meta[keep2, ]
donor_meta$dataset <- factor(donor_meta$dataset)

message(sprintf("    after dataset filter: donors = %d, datasets = %d",
                nrow(donor_meta), nlevels(donor_meta$dataset)))

## (d) 세포수준 검출률 (매칭용) — subsample 로 계산
set.seed(SUBSAMPLE_SEED)
sidx <- sample(ncol(cnt), min(DET_SUBSAMPLE, ncol(cnt)))
cnt_sub  <- as(cnt[, sidx, drop = FALSE], "dgCMatrix")
det_rate <- Matrix::rowMeans(cnt_sub > 0)
rm(cnt_sub); gc()

## (e) 유전자 필터 — 관측/null 전 반복에서 동일한 gene universe 를 쓴다 (중요)
dge0 <- DGEList(pb)
design_base <- model.matrix(~ scale(dm[GENE, ]) + donor_meta$depth + donor_meta$dataset)
keep_g <- filterByExpr(dge0, design = design_base, min.count = PB_FILTER_MIN_COUNT)
dge0 <- dge0[keep_g, , keep.lib.sizes = FALSE]
dge0 <- calcNormFactors(dge0, method = "TMM")
universe <- rownames(dge0)
message(sprintf("    genes after filterByExpr: %d", length(universe)))

# =============================================================================
# §2. 핵심 함수 — 관측/null 이 완전히 동일한 경로를 타도록 하나로 묶는다
# =============================================================================

#' donor 연속형 DE → t 통계량 랭킹 (self-gene 및 SLC7A7 제외)
#' @param g 예측변수로 쓸 유전자명
run_donor_continuous <- function(g) {
  x <- scale(as.numeric(dm[g, colnames(dge0$counts)]))[, 1]
  if (!is.finite(sd(x)) || sd(x) == 0) return(NULL)
  design <- model.matrix(~ x + donor_meta$depth + donor_meta$dataset)
  colnames(design) <- make.names(colnames(design))
  v   <- voom(dge0, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  tt  <- topTable(fit, coef = "x", number = Inf, sort.by = "none")
  st  <- setNames(tt$t, rownames(tt))
  # 순환논리 차단: 예측변수 자신 + 항상 SLC7A7 제거 (universe 동일성 유지)
  st[setdiff(names(st), unique(c(g, GENE)))]
}

#' 랭킹 → Hallmark fgsea → 관심 pathway NES 벡터
nes_of <- function(stats_vec, seed) {
  if (is.null(stats_vec)) return(setNames(rep(NA_real_, length(TARGET_PATHWAYS)),
                                          TARGET_PATHWAYS))
  set.seed(seed)
  res <- suppressWarnings(fgsea(pathways = HALLMARK,
                                stats    = stats_vec,
                                minSize  = GSEA_MIN_SIZE,
                                maxSize  = GSEA_MAX_SIZE,
                                nproc    = NPROC))
  out <- setNames(rep(NA_real_, length(TARGET_PATHWAYS)), TARGET_PATHWAYS)
  hit <- match(TARGET_PATHWAYS, res$pathway)
  out[!is.na(hit)] <- res$NES[hit[!is.na(hit)]]
  out
}

# Hallmark 집합 (관측 분석과 동일 버전 사용)
# msigdbr >= 10 은 collection=, 그 이전은 category= 인자를 씁니다.
HALLMARK <- local({
  m <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
                error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
  split(m$gene_symbol, m$gs_name)
})
stopifnot(length(HALLMARK) == 50)

# =============================================================================
# §3. 발현수준 매칭 대조 유전자 200개 선정
# =============================================================================

message("[2] selecting expression-matched control genes ...")

gstat <- data.frame(
  gene    = universe,
  mean_pb = rowMeans(cpm(dge0, log = TRUE, prior.count = 1)),
  sd_don  = apply(dm[universe, colnames(dge0$counts), drop = FALSE], 1, sd),
  det     = det_rate[universe],
  stringsAsFactors = FALSE
)
gstat <- gstat[is.finite(gstat$mean_pb) & is.finite(gstat$sd_don) &
                 is.finite(gstat$det), ]

# donor 수준 상관 (공조절 유전자 제외용)
x_obs <- as.numeric(dm[GENE, colnames(dge0$counts)])
rho_g <- apply(dm[gstat$gene, colnames(dge0$counts), drop = FALSE], 1,
               function(y) suppressWarnings(cor(x_obs, y, method = "spearman")))
gstat$rho_slc7a7 <- rho_g

# 표준화 3D 공간에서 최근접 이웃
z  <- scale(gstat[, c("mean_pb", "sd_don", "det")])
i0 <- which(gstat$gene == GENE)
if (length(i0) != 1) {
  stop(sprintf(paste0("%s 가 filterByExpr 이후 gene universe 에 없습니다.\n",
                      "  → PB_FILTER_MIN_COUNT 를 낮추거나, filterByExpr 에 ",
                      "'keep[GENE] <- TRUE' 를 강제하세요."), GENE))
}
d2 <- sqrt(rowSums((z - matrix(z[i0, ], nrow(z), 3, byrow = TRUE))^2))
gstat$dist <- d2

## ── caliper 매칭 ────────────────────────────────────────────────────────────
ref <- list(mean_pb = gstat$mean_pb[i0], det = gstat$det[i0], sd_don = gstat$sd_don[i0])

## 두 모집단: 공조절 유전자 제외본(dec) / 미제외본(all)
base_all <- gstat %>% filter(gene != GENE, gene != PARTNER)
base_dec <- base_all %>% filter(abs(rho_slc7a7) < EXCLUDE_COR_ABS)

## 각 모집단 안에서 거리 분위수로 caliper 를 잡는다
cal <- function(df, q) df %>% filter(dist <= quantile(dist, q, na.rm = TRUE))

set_dec_main  <- cal(base_dec, CALIPER_Q)
set_dec_tight <- cal(base_dec, CALIPER_Q_TIGHT)
set_all_main  <- cal(base_all, CALIPER_Q)
set_all_tight <- cal(base_all, CALIPER_Q_TIGHT)

## 실제로 NES 를 적합할 유전자 = 두 주 caliper 의 합집합
fit_genes <- union(set_dec_main$gene, if (RUN_BOTH_RHO) set_all_main$gene else character(0))
if (length(fit_genes) > POOL_CAP) {
  set.seed(GSEA_SEED)
  keep_extra <- setdiff(fit_genes, set_dec_main$gene)          # 주 분석 집합은 반드시 보존
  n_extra    <- max(0, POOL_CAP - length(set_dec_main$gene))
  fit_genes  <- c(set_dec_main$gene, sample(keep_extra, min(n_extra, length(keep_extra))))
  message("    ! POOL_CAP 에 걸려 미제외본 쪽을 일부만 적합합니다")
}
matched <- gstat %>% filter(gene %in% fit_genes) %>% arrange(dist)

rng <- function(d, col) sprintf("%.3f-%.3f", min(d[[col]]), max(d[[col]]))
message(sprintf("    SLC7A7: mean_pb %.3f | det %.4f | sd_don %.3f",
                ref$mean_pb, ref$det, ref$sd_don))
for (nm in c("set_dec_main","set_dec_tight","set_all_main","set_all_tight")) {
  d <- get(nm)
  message(sprintf("    %-14s n=%4d | d<=%.3f | mean %s | det %s",
                  nm, nrow(d), max(d$dist), rng(d,"mean_pb"), rng(d,"det")))
}
message(sprintf("    → 적합 대상 %d개 (두 주 caliper 의 합집합)", nrow(matched)))

## 재현을 위해 멤버십을 그대로 저장
memb <- gstat %>%
  filter(gene %in% matched$gene) %>%
  transmute(gene, mean_pb, sd_don, det, rho_slc7a7, dist,
            in_dec_main  = gene %in% set_dec_main$gene,
            in_dec_tight = gene %in% set_dec_tight$gene,
            in_all_main  = gene %in% set_all_main$gene,
            in_all_tight = gene %in% set_all_tight$gene)
write.csv(memb, file.path(DIR_RESULTS, paste0("Supp_null_matched_genes", SFX, ".csv")),
          row.names = FALSE)

write.csv(data.frame(
  parameter = c("caliper_q","caliper_q_tight","exclude_cor_abs","run_both_rho",
                "slc7a7_mean_pb","slc7a7_det","slc7a7_sd_don",
                "n_base_dec","n_base_all",
                "n_dec_main","n_dec_tight","n_all_main","n_all_tight","n_fitted",
                "dist_cut_dec_main","dist_cut_all_main","n_sub","n_draw"),
  value = c(CALIPER_Q, CALIPER_Q_TIGHT, EXCLUDE_COR_ABS, RUN_BOTH_RHO,
            ref$mean_pb, ref$det, ref$sd_don,
            nrow(base_dec), nrow(base_all),
            nrow(set_dec_main), nrow(set_dec_tight),
            nrow(set_all_main), nrow(set_all_tight), nrow(matched),
            max(set_dec_main$dist), max(set_all_main$dist), N_SUB, N_DRAW)),
  file.path(DIR_RESULTS, paste0("Supp_null_matching_params", SFX, ".csv")),
  row.names = FALSE)

# =============================================================================
# §4. 관측값 + empirical null 200회
# =============================================================================

message("[3] observed run ...")
obs_nes <- nes_of(run_donor_continuous(GENE), seed = GSEA_SEED)
print(round(obs_nes, 3))

N_FIT <- nrow(matched)
message("[4] null runs (적합 ", N_FIT, "개 — 30-45분 소요) ...")
null_mat <- matrix(NA_real_, nrow = N_FIT, ncol = length(TARGET_PATHWAYS),
                   dimnames = list(matched$gene, TARGET_PATHWAYS))
pb_bar <- txtProgressBar(min = 0, max = N_FIT, style = 3)
for (i in seq_len(N_FIT)) {
  null_mat[i, ] <- nes_of(run_donor_continuous(matched$gene[i]),
                          seed = GSEA_SEED + i)
  setTxtProgressBar(pb_bar, i)
}
close(pb_bar)

null_df <- as.data.frame(null_mat) %>%
  tibble::rownames_to_column("control_gene") %>%
  pivot_longer(-control_gene, names_to = "pathway", values_to = "NES")
write.csv(null_df, file.path(DIR_RESULTS, paste0("Supp_null_NES_distribution", SFX, ".csv")),
          row.names = FALSE)

## empirical p — (1 + #{|NES_null| >= |NES_obs|}) / (1 + n)  [Phipson & Smyth 2010]
## empirical p — (1 + #{|NES_null| >= |NES_obs|}) / (1 + n)  [Phipson & Smyth 2010]
emp_for <- function(genes, label) {
  idx <- which(rownames(null_mat) %in% genes)
  lapply(TARGET_PATHWAYS, function(p) {
    nl <- null_mat[idx, p]; nl <- nl[is.finite(nl)]
    ob <- obs_nes[[p]]
    data.frame(
      null_set       = label,
      pathway        = p,
      observed_NES   = ob,
      n_null         = length(nl),
      null_mean      = mean(nl),
      null_sd        = sd(nl),
      null_q025      = unname(quantile(nl, 0.025)),
      null_q975      = unname(quantile(nl, 0.975)),
      z_vs_null      = (ob - mean(nl)) / sd(nl),
      p_emp_two      = (1 + sum(abs(nl) >= abs(ob))) / (1 + length(nl)),
      p_emp_one      = (1 + sum(if (ob >= 0) nl >= ob else nl <= ob)) / (1 + length(nl))
    )
  }) %>% bind_rows()
}

sets <- list(
  "decorrelated_caliper05" = set_dec_main$gene,
  "decorrelated_caliper02" = set_dec_tight$gene)
if (RUN_BOTH_RHO) sets <- c(sets, list(
  "all_caliper05"          = set_all_main$gene,
  "all_caliper02"          = set_all_tight$gene))

emp <- bind_rows(lapply(names(sets), function(s) emp_for(sets[[s]], s)))
emp <- emp %>% group_by(null_set) %>%
  mutate(p_emp_two_BH = p.adjust(p_emp_two, "BH")) %>% ungroup()   # 집합별로 따로 보정

write.csv(emp, file.path(DIR_RESULTS, paste0("Supp_null_empirical_p", SFX, ".csv")),
          row.names = FALSE)
print(as.data.frame(emp))

## 집합 간 비교 요약 — 결론이 매칭 폭이나 공조절 제외에 의존하는지 한눈에
cmp <- emp %>%
  select(null_set, pathway, p_emp_two, p_emp_two_BH) %>%
  tidyr::pivot_wider(names_from = null_set, values_from = c(p_emp_two, p_emp_two_BH))
write.csv(cmp, file.path(DIR_RESULTS, paste0("Supp_null_set_comparison", SFX, ".csv")),
          row.names = FALSE)
print(as.data.frame(cmp))

## 그림 — pathway별 null 분포 + 관측값 위치
pdf(file.path(DIR_FIGURES, paste0("SuppX_null_NES_hist", SFX, ".pdf")), width = 9, height = 6)
print(
  ggplot(null_df, aes(NES)) +
    geom_histogram(bins = 40, fill = "grey75", colour = "white") +
    geom_vline(data = emp, aes(xintercept = observed_NES),
               colour = PALETTE_TWO[1], linewidth = 1) +
    geom_text(data = emp, aes(x = observed_NES, y = Inf,
                              label = sprintf("P[emp] == %.3g", p_emp_two)),
              parse = TRUE, vjust = 1.6, hjust = -0.05, size = 3) +
    facet_wrap(~ pathway, scales = "free", ncol = 3) +
    labs(x = sprintf("NES from %d matched control genes (caliper %g, co-regulated excluded)",
                     nrow(set_dec_main), CALIPER_Q),
         y = "count") +
    theme_bw(base_size = 9)
)
dev.off()


## ── 대조 유전자 선택에 대한 민감도 ──────────────────────────────────────────
## 주 caliper 안에서 N_SUB 개를 비복원 추출하는 것을 N_DRAW 회 반복한다.
## 저장된 NES 를 쓰므로 추가 적합이 없다.
stab <- bind_rows(lapply(names(sets), function(s) {
  idx_all <- which(rownames(null_mat) %in% sets[[s]])
  if (length(idx_all) <= N_SUB) return(NULL)
  set.seed(GSEA_SEED + 10000L)
  draws <- replicate(N_DRAW, sample(idx_all, N_SUB), simplify = FALSE)
  bind_rows(lapply(TARGET_PATHWAYS, function(p) {
    ob <- obs_nes[[p]]
    ps <- vapply(draws, function(ix) {
      nl <- null_mat[ix, p]; nl <- nl[is.finite(nl)]
      (1 + sum(abs(nl) >= abs(ob))) / (1 + length(nl))
    }, numeric(1))
    data.frame(null_set = s, pathway = p, n_sub = N_SUB, n_draw = N_DRAW,
               p_emp_min = min(ps), p_emp_median = median(ps),
               p_emp_max = max(ps), p_emp_sd = sd(ps),
               n_draw_sig_0.05 = sum(ps < 0.05))
  }))
}))
write.csv(stab, file.path(DIR_RESULTS, paste0("Supp_null_draw_stability", SFX, ".csv")),
          row.names = FALSE)
print(as.data.frame(stab))


# =============================================================================
# §5. SLC3A2 공발현의 매칭 대조 (§9.1 depth 보정 null 대응)
# =============================================================================

message("[5] SLC3A2 partial-correlation null ...")

# SLC3A2 에 매칭된 200개 (SLC7A7 매칭과 별개)
iP <- which(gstat$gene == PARTNER)
if (length(iP) == 1) {
  dP <- sqrt(rowSums((z - matrix(z[iP, ], nrow(z), 3, byrow = TRUE))^2))
  gstatP   <- gstat %>% filter(gene != PARTNER, gene != GENE) %>% mutate(distP = dP[gene])
  matchedP <- gstatP %>% filter(distP <= quantile(distP, CALIPER_Q, na.rm = TRUE))
  message(sprintf("    SLC3A2 null: caliper %g 안의 %d개 사용", CALIPER_Q, nrow(matchedP)))
  
  t_of <- function(g) {
    y <- as.numeric(dm[g, colnames(dge0$counts)])
    coef(summary(lm(y ~ x_obs + donor_meta$depth)))["x_obs", "t value"]
  }
  t_obs  <- t_of(PARTNER)
  t_null <- vapply(matchedP$gene, t_of, numeric(1))
  
  outP <- data.frame(
    gene       = PARTNER,
    t_observed = t_obs,
    n_null     = length(t_null),
    null_mean  = mean(t_null), null_sd = sd(t_null),
    percentile = mean(t_null < t_obs) * 100,
    p_emp_two  = (1 + sum(abs(t_null) >= abs(t_obs))) / (1 + length(t_null)),
    p_emp_one  = (1 + sum(t_null >= t_obs)) / (1 + length(t_null))
  )
  write.csv(outP, file.path(DIR_RESULTS, paste0("Supp_null_SLC3A2_partialcor", SFX, ".csv")),
            row.names = FALSE)
  print(outP)
}

# =============================================================================
# §6. [legacy] 세포단위 이분 분석의 음성대조 — 폐기 근거로 보존
# -----------------------------------------------------------------------------
#  이 블록은 "왜 세포단위 이분을 버렸는가"의 증거이므로 지우지 말고 남긴다.
#  구버전과 달리 매칭 유전자 200개를 사용하여 empirical p 를 재계산한다.
#  기대 결과: oxphos empirical p ≈ 0.5 (통과 실패) → 폐기 결정의 정당화
# =============================================================================

if (file.exists(file.path(DIR_RESULTS, "Supp_depth_sensitivity.csv"))) {
  file.rename(file.path(DIR_RESULTS, "Supp_depth_sensitivity.csv"),
              file.path(DIR_RESULTS, "Supp_depth_sensitivity_legacy.csv"))
  message("[6] 기존 Supp_depth_sensitivity.csv → *_legacy.csv 로 이동")
}
# (원본 세포단위 코드는 supplementary/celllevel_slc7a7_deg.R 로 이관 권장)

message("done. 논문에 넣을 표: results/Supp_null_empirical_p.csv")
