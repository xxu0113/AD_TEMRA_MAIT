# =====================================================================
# CSF scRNA-seq (GSE134577) — QC, clustering, MAIT module score
# Outputs (files): outputs/csf/
#   - csf_sct.rds
#   - csf_tsne_clusters.png
#   - csf_tsne_by_sample.png
#   - csf_violin_MAITscore.png
#   - csf_dot_MAITmarkers.png
#   - csf_cluster_mean_MAITscore.tsv
# Requirements in data/scRNA_CSF/<sample>/ :
#   matrix.mtx , barcodes.tsv , genes.tsv (or features.tsv)
# =====================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(data.table)
})

# --------------------------
# 0) Paths
# --------------------------
dir_csf <- "data/scRNA_CSF"
dir_out <- "outputs/csf"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

# --------------------------
# 1) Load all CSF samples
# --------------------------
all_dirs <- list.dirs(dir_csf, recursive = FALSE, full.names = TRUE)
has_10x <- function(d) {
  file.exists(file.path(d, "matrix.mtx")) &&
    file.exists(file.path(d, "barcodes.tsv")) &&
    (file.exists(file.path(d, "genes.tsv")) || file.exists(file.path(d, "features.tsv")))
}
sample_dirs <- Filter(has_10x, all_dirs)
stopifnot(length(sample_dirs) > 0)

objs <- lapply(sample_dirs, function(sd){
  m  <- Read10X(data.dir = sd)     # handles genes.tsv/features.tsv; gz or not
  id <- basename(sd)
  so <- CreateSeuratObject(m, project = id)
  so$sample_id <- id
  colnames(so) <- paste(id, colnames(so), sep = "_")  # make colnames unique
  so
})

csf <- Reduce(function(a,b) merge(a,b), objs)

# --------------------------
# 2) QC (paper settings for CSF)
# --------------------------
csf[["percent.mt"]] <- PercentageFeatureSet(csf, pattern = "^MT-")

csf <- subset(
  csf,
  subset = nFeature_RNA >= 200 & nFeature_RNA <= 1800 &
    nCount_RNA   <= 7500 & percent.mt   <= 10
)

# Keep genes expressed in ≥10 cells (Seurat v5 layer-safe)
rna <- csf[["RNA"]]
avail_layers <- Layers(rna)
layer_name <- if ("counts" %in% avail_layers) "counts" else
  if ("raw" %in% avail_layers)    "raw"    else avail_layers[1]
message("Using RNA layer: ", layer_name)

counts <- GetAssayData(rna, layer = layer_name)  # sparse dgCMatrix
keep_genes <- Matrix::rowSums(counts > 0) >= 10
csf <- csf[keep_genes, ]

message(sprintf("After QC: %d genes x %d cells", nrow(csf), ncol(csf)))

# --------------------------
# 3) Normalization (SCTransform) — minimal & robust
# --------------------------
DefaultAssay(csf) <- "RNA"

# Make sure regressors are numeric
csf$percent.mt  <- as.numeric(csf$percent.mt)
csf$nCount_RNA  <- as.numeric(csf$nCount_RNA)

# Use glmGamPoi if available
sct_method <- if (requireNamespace("glmGamPoi", quietly = TRUE)) "glmGamPoi" else "vst"
message("SCTransform method: ", sct_method, " | Regressors: percent.mt, nCount_RNA")

csf <- SCTransform(
  csf,
  vars.to.regress = c("percent.mt", "nCount_RNA"),  # ← no sample_id here
  method = sct_method,
  verbose = FALSE
)

# --------------------------
# 4) Dimensionality reduction + clustering
# --------------------------
csf <- RunPCA(csf, verbose = FALSE)
csf <- FindNeighbors(csf, dims = 1:4)      # 4 PCs
csf <- FindClusters(csf, resolution = 0.3) # res = 0.3
csf <- RunTSNE(csf, dims = 1:4)

# --------------------------
# 5) MAIT module score
# --------------------------
mait_genes <- list(c("KLRB1","SLC4A10","ZBTB16","TRAV1-2","IL18RAP","CCR6","CXCR6"))
csf <- AddModuleScore(csf, features = mait_genes, name = "MAIT_Score")

# --------------------------
# 6) Save object + plots (also print to Plots pane)
# --------------------------
saveRDS(csf, file.path(dir_out, "csf_sct.rds"))

p_tsne_clusters <- DimPlot(csf, reduction = "tsne", group.by = "seurat_clusters", label = TRUE)
p_tsne_sample   <- DimPlot(csf, reduction = "tsne", group.by = "sample_id")
p_violin        <- VlnPlot(csf, features = "MAIT_Score1", group.by = "seurat_clusters", pt.size = 0)
p_dot           <- DotPlot(csf, features = unlist(mait_genes)) + RotatedAxis()

# show interactively
print(p_tsne_clusters); print(p_tsne_sample); print(p_violin); print(p_dot)

# save to disk
ggsave(file.path(dir_out, "csf_tsne_clusters.png"), p_tsne_clusters, width = 6, height = 5, dpi = 300)
ggsave(file.path(dir_out, "csf_tsne_by_sample.png"),   p_tsne_sample,   width = 6, height = 5, dpi = 300)
ggsave(file.path(dir_out, "csf_violin_MAITscore.png"), p_violin,        width = 7, height = 4, dpi = 300)
ggsave(file.path(dir_out, "csf_dot_MAITmarkers.png"),  p_dot,           width = 7, height = 4, dpi = 300)

# --------------------------
# 7) Export per-cluster MAIT score summary
# --------------------------
Idents(csf) <- "seurat_clusters"
cluster_means <- data.frame(cluster = Idents(csf), score = csf$MAIT_Score1) |>
  dplyr::group_by(cluster) |>
  dplyr::summarise(mean_MAIT = mean(score, na.rm = TRUE), .groups = "drop")
cluster_means
fwrite(cluster_means, file.path(dir_out, "csf_cluster_mean_MAITscore.tsv"), sep = "\t")

message("Done. Outputs saved to: ", normalizePath(dir_out))

# I think cluster7 expresses some mait cell marker but not strong, but I still make analysis
# manually map GSM -> diagnosis (from GEO page)
lookup <- c(
  GSM3984199_CSF1 = "HC",
  GSM3984200_CSF2 = "AD",
  GSM3984201_CSF3 = "HC",
  GSM3984202_CSF4 = "HC",
  GSM3984203_CSF5 = "MCI",
  GSM3984204_CSF6 = "AD",
  GSM3984205_CSF8 = "HC",
  GSM3984206_CSF9 = "AD",
  GSM3984207_CSF11= "MCI",
  GSM3984208_CSF12= "HC",
  GSM3984209_CSF13= "HC",
  GSM3984210_CSF15= "MCI",
  GSM3984211_CSF16= "HC",
  GSM3984212_CSF17= "MCI",
  GSM3984213_CSF18= "AD",
  GSM3984214_CSF19= "MCI",
  GSM3984215_CSF20= "HC",
  GSM3984216_CSF23= "HC"
)
csf$Diagnosis <- unname(lookup[as.character(csf$sample_id)])
table(csf$Diagnosis)


Idents(csf) <- "seurat_clusters"

# raw counts
table(csf$Diagnosis[Idents(csf) == "7"])

# proportions per diagnosis group
prop.table(table(csf$Diagnosis, Idents(csf)), margin = 1)[,"7"]

mait_markers <- c("KLRB1","SLC4A10","ZBTB16","TRAV1-2","IL18RAP","CCR6","CXCR6")

FeaturePlot(csf, features = mait_markers, reduction = "tsne", pt.size = 0.3)
DotPlot(csf, features = mait_markers) + RotatedAxis()





