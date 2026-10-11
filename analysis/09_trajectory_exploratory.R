# ============================================================================
# 09_trajectory_exploratory.R
# ----------------------------------------------------------------------------
# PURPOSE  Exploratory trajectory screen. Descriptive only. No statistical
#          test is applied and no claim in the paper rests on it: the ordering
#          is built across cells rather than within donors, it is not
#          independent of dataset of origin, and detection of SLC7A7 at the
#          cell level is depth-dependent.
#
# INPUT    BPCells matrix
# OUTPUT   Fig 6 panels
#
# Sources R/setup.R, which loads config/config.R. See docs/PIPELINE.md for
# the execution order and docs/OUTPUT_MAP.md for the mapping to additional
# files and figures. Comments below this block are the working notes kept
# from development and are in Korean.
# ============================================================================

# =============================================================================
# 09_trajectory_exploratory.R               [Fig 7]
# monocle3 pseudotime and slingshot cross-check. Descriptive only.
# ※ prolif 클러스터는 trajectory 에서 제외 (세포주기 신호가 활성화 연속체를 교란).
#    Methods 에 근거 명시. 비율/CCC 분석에는 포함.
# =============================================================================
cds <- readRDS(file.path(DIR_DATA_PROC, "cds_microglia_all.rds"))

source(here::here("R", "setup.R"))
suppressPackageStartupMessages({
  library(Seurat); library(monocle3); library(slingshot); library(SingleCellExperiment)

mg <- readRDS(file.path(DIR_DATA_PROC, "score_auc_microglia.rds"))

set.seed(SUBSAMPLE_SEED)
sm <- mg
# if ("annotation_level_4" %in% colnames(sm@meta.data))
#   sm <- subset(sm, subset = annotation_level_4 != "TAM-MG prolif")
cnt <- as(GetAssayData(sm, assay = "RNA", layer = "data"), "dgCMatrix")
# counts X, data O. normalized

cds <- new_cell_data_set(cnt, cell_metadata = sm@meta.data,
        gene_metadata = data.frame(gene_short_name = rownames(cnt), row.names = rownames(cnt)))
cds <- preprocess_cds(cds, num_dim = 30)
dn  <- table(colData(cds)$donor_id); cds <- cds[, colData(cds)$donor_id %in% names(dn[dn >= 20])]
colData(cds)$donor_id <- droplevels(as.factor(colData(cds)$donor_id))
cds <- align_cds(cds, alignment_group = "donor_id",
                 alignment_k = min(20, min(table(colData(cds)$donor_id)) - 1))
cds <- reduce_dimension(cds, preprocess_method = "Aligned")
cds <- cluster_cells(cds, resolution = 1e-4)
cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE)

hm <- intersect(c("P2RY12","TMEM119","CX3CR1","SALL1","MEF2C"), rownames(cds))
cds$homeo <- colMeans(as.matrix(exprs(cds)[hm, , drop = FALSE]))
root_cl <- names(sort(tapply(cds$homeo, monocle3::clusters(cds), mean), decreasing = TRUE))[1]
cds <- order_cells(cds, root_cells = colnames(cds)[monocle3::clusters(cds) == root_cl])


# donor batch confirm -----------------------------------------------------


plot_cells(
  cds,
  color_cells_by = "donor_id",
  show_trajectory_graph = TRUE
)


# scoring midules ---------------------------------------------------------


names <- c(
  "IFN_ALPHA",
  "IFN_GAMMA",
  "NFKB",
  "INFLAMMATORY",
  "TRANSPORT",
  "ARG_ENZYME",
  "HYPOXIA",
  "GLYCOLYSIS",
  "OXPHOS",
  "ROS",
  "IFN_REACTOME",
  "PPP_REACTOME",
  "POLYAMINE",
  "PAN_MYELOID",
  "HOMEOSTATIC",
  "IFN_ALPHA_disj",
  "IFN_GAMMA_disj",
  "NFKB_disj",
  "INFLAMMATORY_disj",
  "TRANSPORT_disj",
  "ARG_ENZYME_disj",
  "HYPOXIA_disj",
  "GLYCOLYSIS_disj",
  "OXPHOS_disj",
  "ROS_disj",
  "IFN_REACTOME_disj",
  "PPP_REACTOME_disj",
  "POLYAMINE_disj",
  "PAN_MYELOID_disj",
  "HOMEOSTATIC_disj"
)

for (name in names) {
  
  if (!name %in% colnames(mg@meta.data)) {
    warning(sprintf("'%s'가 mg@meta.data에 없습니다. 건너뜁니다.", name))
    next
  }
  
  colData(cds)[[name]] <- mg@meta.data[colnames(cds), name]
}


#annotation4에서 proliferate 제거 안 한 것

# colData(cds) <- cbind(
#   colData(cds),
#   U$mr3000[colnames(cds), ]
# )
saveRDS(cds, file.path(DIR_DATA_PROC, "cds_microglia_all.rds"))


for (name in names) {
  # plot
  p <- plot_cells(
    cds,
    color_cells_by = name,
    cell_size = 0.6,
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    label_roots = FALSE
  )
  
  # 저장
  ggsave(
    file.path(DIR_FIGURES, paste0("Fig7B_monocle3_", name, ".pdf")),
    plot = p,
    width = 6,
    height = 5
  )
}

for (name in c("TRANSPORT", "ARG_ENZYME", "IFN_ALPHA", "NFKB", "HYPOXIA")){
  p <- plot_cells(
    cds,
    color_cells_by = name,
    cell_size = 0.6,
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    label_roots = FALSE
  )
  p$layers[[2]]$aes_params$linewidth <- 2
  save_emf(p + theme_blank_frame(), paste0("Fig7B_monocle3_", name, ".emf"), legend_p = "right")
}


# ★ 연속값 투영 (이분 금지)
ggplot2::ggsave(file.path(DIR_FIGURES, "Fig7A_monocle3_gene_continuous.pdf"),
  plot_cells(cds, genes = GENE_OF_INTEREST, cell_size = .6, label_cell_groups = FALSE,
             label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE) +
    ggplot2::scale_color_viridis_c(), width = 6, height = 5)

save_emf(plot_cells(cds, genes = GENE_OF_INTEREST, cell_size = 1.2, label_cell_groups = FALSE,
                    label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE)+theme_blank_frame(), "Fig7A_monocle3_gene_continuous.emf", legend_p = "right")



ggplot2::ggsave(file.path(DIR_FIGURES, "Fig7B_monocle3_pseudotime.pdf"),
  plot_cells(cds, color_cells_by = "pseudotime", cell_size = .6, label_cell_groups = FALSE,
             label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE),
  width = 6, height = 5)
save_emf(plot_cells(cds, color_cells_by = "pseudotime", cell_size = 1.2, label_cell_groups = FALSE,
                    label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE)+theme_blank_frame(), "Fig7B_monocle3_pseudotime.emf", legend_p = "right")


for (name in c("SLC7A7", "IL1B", "CCL4", "CXCL10")) {
  # plot
  p <- plot_genes_in_pseudotime(
    cds[name, ],
    min_expr = 0.1,
    cell_size = 1.5
  )
  
  # 저장
  ggsave(
    file.path(DIR_FIGURES, paste0("Fig7D_monocle3_genebypseudotime_", name, ".pdf")),
    plot = p,
    width = 6,
    height = 5
  )
}

p <- plot_genes_in_pseudotime(
  cds["SLC7A7", ],
  min_expr = 0.1,
  cell_size = 3
)
p$layers[[2]]$aes_params$linewidth <- 2
p + theme_blank_frame()
save_emf(p + theme_blank_frame(), "Fig7D_genebypseudotime_SLC7A7.emf",legend_p = "right")
# # donor 수준 상관 (pseudoreplication 회피)
# pt <- data.frame(pt = colData(cds)$pseudotime,
#                  expr = as.numeric(exprs(cds)[GENE_OF_INTEREST, ]),
#                  donor = colData(cds)$donor_id) %>%
#   dplyr::filter(is.finite(pt)) %>% dplyr::group_by(donor) %>%
#   dplyr::summarise(pt = mean(pt), expr = mean(expr), .groups = "drop")
# print(stats::cor.test(pt$pt, pt$expr, method = "spearman"))
# write_result(pt, "Fig7_donor_pseudotime_vs_gene.csv")
# 
# ## Slingshot 교차검증
# sce <- as.SingleCellExperiment(NormalizeData(CreateSeuratObject(cnt, meta.data = sm@meta.data)))
# reducedDims(sce)$UMAP <- Embeddings(sm, "umap")[colnames(sce), ]
# sce$cluster <- sm$annotation_level_4[colnames(sce)]
# sce <- slingshot(sce, clusterLabels = "cluster", reducedDim = "UMAP", start.clus = "TAM-MG aging sig")
# saveRDS(sce, file.path(DIR_DATA_PROC, "slingshot_pseudotime.rds"))
#
