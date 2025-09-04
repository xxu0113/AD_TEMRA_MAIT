# =====================================================================
# BLOOD scRNA-seq (GSE134576) — QC, clustering, MAIT markers
# Outputs: plots show in RStudio Plots pane, tables print in console
# =====================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

dir_blood <- "data/scRNA_BLOOD"
dir_out   <- "outputs/blood"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

# ---------- 1) Load samples ----------
all_dirs <- list.dirs(dir_blood, recursive = FALSE, full.names = TRUE)
has_10x <- function(d) {
  file.exists(file.path(d,"matrix.mtx")) &&
    file.exists(file.path(d,"barcodes.tsv")) &&
    (file.exists(file.path(d,"genes.tsv")) || file.exists(file.path(d,"features.tsv")))
}
sample_dirs <- Filter(has_10x, all_dirs)
stopifnot(length(sample_dirs) > 0)

objs <- lapply(sample_dirs, function(sd){
  m  <- Read10X(sd)
  id <- basename(sd)
  so <- CreateSeuratObject(m, project = id)
  so$sample_id <- id
  colnames(so) <- paste(id, colnames(so), sep = "_")  # unique barcodes
  so
})
blood <- Reduce(function(a,b) merge(a,b), objs)
DefaultAssay(blood) <- "RNA"

# ---------- 2) QC ----------
blood[["percent.mt"]] <- PercentageFeatureSet(blood, pattern = "^MT-")
blood <- subset(
  blood,
  subset = nFeature_RNA >= 200 & nFeature_RNA <= 2500 &
    nCount_RNA   <= 10000 & percent.mt <= 10
)

# keep genes in >=10 cells
rna <- blood[["RNA"]]
avail <- Layers(rna)
layer <- if ("counts" %in% avail) "counts" else if ("raw" %in% avail) "raw" else avail[1]
counts <- GetAssayData(rna, layer = layer)
keep_genes <- Matrix::rowSums(counts > 0) >= 10
blood <- blood[keep_genes, ]
message(sprintf("After QC: %d genes x %d cells", nrow(blood), ncol(blood)))

# ---------- 3) SCTransform ----------
blood$percent.mt <- as.numeric(blood$percent.mt)
blood$nCount_RNA <- as.numeric(blood$nCount_RNA)
sct_method <- if (requireNamespace("glmGamPoi", quietly = TRUE)) "glmGamPoi" else "vst"
message("SCTransform method: ", sct_method)

blood <- SCTransform(
  blood,
  vars.to.regress = c("percent.mt","nCount_RNA"),
  method = sct_method,
  verbose = FALSE
)

# ---------- 4) PCA / neighbors / clustering / tSNE ----------
blood <- RunPCA(blood, verbose = FALSE)
blood <- FindNeighbors(blood, dims = 1:10)  # paper used 10 PCs for blood
blood <- FindClusters(blood, resolution = 0.4)
blood <- RunTSNE(blood, dims = 1:10)

# ---------- 5) Visualize MAIT markers ----------
mait_genes <- c("KLRB1","SLC4A10","ZBTB16","TRAV1-2","IL18RAP","CCR6","CXCR6")
present <- mait_genes[mait_genes %in% rownames(blood)]
message("Present markers: ", paste(present, collapse = ", "))

# FeaturePlots
for (g in present) {
  print(FeaturePlot(blood, features = g, reduction = "tsne", pt.size = 0.3))
  # ggsave(file.path(dir_out, paste0("blood_tsne_", g, ".png")), last_plot())
}

# DotPlot overview
Idents(blood) <- "seurat_clusters"
p_dot <- DotPlot(blood, features = present) + RotatedAxis()
print(p_dot)
ggsave(file.path(dir_out, "blood_dot_MAITmarkers.png"), p_dot)

# ---------- 6) Save object (commented) ----------
saveRDS(blood, file.path(dir_out, "blood_sct.rds"))

message("Done. All plots shown in Plots pane.")

# I think cluster5 is mait cells
#Make a lookup table inside R
lookup_blood <- c(
  GSM3956366_T1 = "AD",
  GSM3956367_T2 = "MCI",
  GSM3956368_T3 = "HC",
  GSM3956369_T4 = "HC",
  GSM3956370_T5 = "AD",
  GSM3956371_T6 = "AD",
  GSM3956372_T7 = "HC",
  GSM3956373_T8 = "HC",
  GSM3956374_T9 = "AD",
  GSM3956375_T10 = "HC",
  GSM3956376_T11 = "HC",
  GSM3956377_T12 = "HC",
  GSM3956378_T14 = "MCI"
)
#Add diagnosis to Seurat metadata
blood$Diagnosis <- unname(lookup_blood[as.character(blood$sample_id)])
table(blood$Diagnosis)   # sanity check

#count cluster5 by diganosis
Idents(blood) <- "seurat_clusters"

# Raw counts of cluster 5 cells by diagnosis
table(blood$Diagnosis[Idents(blood) == "5"])

# Or proportion of each diagnosis group in cluster 5
prop.table(table(blood$Diagnosis, Idents(blood)), margin = 1)[, "5"]