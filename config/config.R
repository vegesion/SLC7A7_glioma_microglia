# =============================================================================
# config.R  ---  프로젝트 전역 설정 / Project-wide configuration
# -----------------------------------------------------------------------------
# 모든 경로 · 파라미터 · 임계값을 이 파일 한 곳에서 관리합니다.
# 개별 스크립트에서 절대경로를 하드코딩하지 마세요. here::here() 기준 상대경로 사용.
# All paths / parameters / thresholds live here. Never hard-code absolute paths
# in analysis scripts; everything is relative to the project root via here().
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
  library(here)
})

# ── 연구 대상 / Study target ─────────────────────────────────────────────────
# 축(axis): 아르기닌 "대사효소"(선행연구, ASL)가 아니라 "수송체"(y+L/CAT 계열).
GENE_OF_INTEREST <- "SLC7A7"   # 사전지정 후보 (y+LAT1)
GENE_PARTNER     <- "SLC3A2"   # 이형이합체 파트너 (4F2hc)
GENE_FAMILY      <- "SLC"      # 스크리닝 대상 접두사

MSIGDB_TRANSPORT_SETS <- c("GOBP_L_ARGININE_IMPORT_ACROSS_PLASMA_MEMBRANE",
                           "GOBP_L_AMINO_ACID_TRANSPORT",
                           "GOBP_L_ARGININE_TRANSMEMBRANE_TRANSPORT")
ARG_ENZYME_SET <- c("GOBP_ARGININE_METABOLIC_PROCESS",
                    "GOBP_L_ARGININE_BIOSYNTHETIC_PROCESS",
                    "GOBP_L_ARGININE_CATABOLIC_PROCESS")

# ── 디렉터리 / Directories (all relative to project root) ────────────────────
DIR_DATA_RAW   <- here("data", "raw")        # 원본 데이터 (git 미추적)
DIR_DATA_PROC  <- here("data", "processed")  # 중간 산출 CSV/RData (git 미추적)
DIR_RESULTS    <- here("results")            # 결과 표/RDS
DIR_FIGURES    <- here("figures")            # 내보낸 그림 (tiff/pdf/emf)

for (d in c(DIR_DATA_RAW, DIR_DATA_PROC, DIR_RESULTS, DIR_FIGURES)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- External data location -------------------------------------------------
# The repository ships with no data. Set SLC7A7_DATA_ROOT in .Renviron if your
# copy of the primary files lives outside the project; leave it unset when the
# files are under data/raw or data/processed, which are searched first.
# See docs/DATA_ACCESS.md for accessions, versions and access dates.
# Do not point this at a cloud-synced drive: file.exists() then waits on the
# network and the session appears to hang.
EXTERNAL_DATA_ROOT <- Sys.getenv("SLC7A7_DATA_ROOT", unset = "")

DATA_SEARCH_ROOTS <- c(
  DIR_DATA_RAW,                                      # 1) repository data/raw
  DIR_DATA_PROC,                                     # 2) repository data/processed
  if (nzchar(EXTERNAL_DATA_ROOT)) EXTERNAL_DATA_ROOT else NULL
)

# 후보 경로 중 실제 존재하는 첫 파일 반환. 없으면 default(첫 후보) 반환(에러 없음).
first_existing <- function(filename, roots = DATA_SEARCH_ROOTS) {
  cand <- file.path(roots, filename)
  hit  <- cand[file.exists(cand)]
  if (length(hit)) hit[1] else cand[1]
}
# 실행 시점에 파일이 반드시 있어야 하는 경우(명확한 에러 메시지). 스크립트에서 사용 권장.
require_data <- function(filename, roots = DATA_SEARCH_ROOTS) {
  cand <- file.path(roots, filename)
  hit  <- cand[file.exists(cand)]
  if (length(hit)) return(hit[1])
  stop(sprintf(
    "데이터 파일을 찾지 못했습니다: %s\n  탐색한 위치:\n%s\n  → data/raw/에 넣거나 config의 DATA_SEARCH_ROOTS를 수정하세요.",
    filename, paste0("    - ", cand, collapse = "\n")), call. = FALSE)
}

# ── 주요 입력 파일 경로 / Key input file paths ───────────────────────────────
# config 로드 시점엔 에러 없이 "찾은 경로(없으면 첫 후보)"로 잡고,
# 실제 읽을 때는 스크립트에서 require_data() 로 존재를 강제할 수 있습니다.
PATH_LGG_NORM <- first_existing("LGG_TCGA_norm.csv")   # LGG 정규화 발현행렬(logCPM)
PATH_GBM_NORM <- first_existing("GBM_TCGA_norm.csv")   # GBM 정규화 발현행렬
# GBmap h5ad, kept under its original CELLxGENE version UUID as the file name.
PATH_H5AD     <- first_existing("861acfd8-25f0-418b-a445-aa96da232827.h5ad")

# GDC 원본 다운로드 폴더 (한 번 받으면 여기 영구 보관 → 재다운로드 방지). git 미추적.
DIR_GDC <- file.path(DIR_DATA_RAW, "GDCdata")

# TCGA 다운로드 방식/안정화 파라미터
#  - "api"(기본): TCGAbiolinks 내장 다운로드. 백신/프록시가 큰 tar.gz 를 끊으면 실패할 수 있음
#                 → files.per.chunk 를 작게(1~2) 하면 각 청크가 작아 통과 확률↑
#  - "client": GDC 공식 gdc-client 사용(재개 가능, 대용량에 가장 안정적). PATH 에 gdc-client 필요.
TCGA_DOWNLOAD_METHOD <- "api"   # 계속 잘리면 "client" 로 변경
TCGA_FILES_PER_CHUNK <- 2       # api 방식 청크당 파일 수 (작을수록 안전, 느려짐). 실패 시 1 로.

# ── 중간 산출물 / Intermediate outputs (파이프라인이 생성, 항상 저장소 내부) ──
PATH_LGG_COX_ALL <- file.path(DIR_DATA_PROC, "LGG_all_cox_zscore.csv")  # written by 02b
PATH_GBM_COX_ALL <- file.path(DIR_DATA_PROC, "GBM_all_cox_zscore.csv")  # written by 02b

# ── 다중 코호트 / Validation cohorts ─────────────────────────────────────────
COHORTS <- list(
  TCGA_LGG = list(expr = "LGG_TCGA_norm.csv",  clin = NULL, note = "written by 01a"),
  TCGA_GBM = list(expr = "GBM_TCGA_norm.csv",  clin = NULL, note = "written by 01a"),
  CGGA_693 = list(expr = "CGGA_693_norm.txt",
                  clin = "CGGA_693_clinical.txt"),
  CGGA_325 = list(expr = "CGGA_325_norm.txt",
                  clin = "CGGA_325_clinical.txt")
)
TRAIN_COHORT <- "TCGA_LGG"

# ── LASSO Cox / signature ────────────────────────────────────────────────────
LASSO_SEED   <- 42
LASSO_NFOLDS <- 10
LASSO_ALPHA  <- 1
LASSO_LAMBDA <- "lambda.min"
ROC_TIMES    <- c(12, 36, 60)

# ── 단일세포: donor 수준 분석 ────────────────────────────────────────────────
SC_DONOR_MIN_CELLS <- 20
AUCELL_SEED        <- 42
DEPTH_COVARIATE    <- "DEPTH" #"DEPTH" <- log1p("nCount_RNA")
CELLCHAT_SUBSAMPLE <- 20000
SC_BATCH_CANDIDATES <- c("dataset","author","study","batch","sample_source")

# ── 공간전사체 ───────────────────────────────────────────────────────────────
DIR_SPATIAL <- file.path(DIR_DATA_RAW, "spatial_data_visium")

# ── 생존분석 임계값 / Survival thresholds ────────────────────────────────────
# EPV(events-per-variable) >= 10 원칙에 따른 최소 이벤트 수.
COX_SCALE_EXPR    <- TRUE   # logCPM을 Z-score로 표준화 (HR 해석 일관성)
COX_MIN_EVENTS_LGG_OVERALL <- 50   # 5 covariates
COX_MIN_EVENTS_LGG_MUT     <- 40   # 4 covariates
COX_MIN_EVENTS_LGG_WT      <- 30   # 3 covariates
COX_MIN_EVENTS_GBM         <- 20   # 2~3 covariates
COX_SIG_P         <- 0.05   # forest plot 유의 필터

# ── GSEA 파라미터 / GSEA parameters ──────────────────────────────────────────
GSEA_MIN_SIZE <- 15
GSEA_MAX_SIZE <- 500
GSEA_PVAL     <- 0.05
GSEA_SEED     <- 42         # fgsea 재현성 (pseudobulk)

# ── scRNA / pseudobulk 파라미터 ──────────────────────────────────────────────
SCRNA_POS_THRESHOLD  <- 0   # SLC7A7 발현 > 0 이면 Positive
PB_MIN_CELLS         <- 10  # donor×status 조합 최소 세포 수
PB_FILTER_MIN_COUNT  <- 3   # edgeR::filterByExpr min.count
SUBSAMPLE_SEED       <- 42

# ── 그림 색상 팔레트 / Figure palette (논문 전반 통일) ───────────────────────
PALETTE_TWO   <- c("#D85A30", "#185FA5")               # High/Low, pos/neg 등 2군
PALETTE_GROUP <- c(Overall = "#D85A30", Mutant = "#BA7517", WT = "#185FA5")
PALETTE_STATUS <- c(Positive = "#E15759", Negative = "#4E79A7")
