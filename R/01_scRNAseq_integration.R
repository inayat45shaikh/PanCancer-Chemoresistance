#!/usr/bin/env Rscript

# 01_scRNAseq_integration.R
# Pan-cancer scRNA-seq QC and Seurat integration
#
# Cleaned from the original analysis script for public release.
# Scientific parameters explicitly present in the original script are preserved.

suppressPackageStartupMessages({
  library(Seurat)
})

set.seed(1234)

# ----------------------------- CONFIGURATION -----------------------------
manifest_file <- "config/sample_manifest_template.csv"
output_dir <- "results/scRNAseq_integration"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

min_features <- 200
max_percent_mt <- 25
nfeatures_integration <- 3000
npcs <- 30
cluster_resolution <- 0.5
neighbor_k <- 10
# ------------------------------------------------------------------------

stopifnot(file.exists(manifest_file))
manifest <- read.csv(manifest_file, stringsAsFactors = FALSE)

required_cols <- c("sample_id", "path", "cancer_type", "response")
missing_cols <- setdiff(required_cols, colnames(manifest))
if (length(missing_cols) > 0) {
  stop("Manifest is missing required columns: ", paste(missing_cols, collapse = ", "))
}

missing_files <- manifest$path[!file.exists(manifest$path)]
if (length(missing_files) > 0) {
  stop(
    "The following input files were not found:\n",
    paste(missing_files, collapse = "\n"),
    "\nEdit ", manifest_file, " before running."
  )
}

message("Loading ", nrow(manifest), " Seurat objects...")
obj_list <- setNames(lapply(manifest$path, readRDS), manifest$sample_id)

prepare_object <- function(obj, sample_id, cancer_type, response) {
  if (!inherits(obj, "Seurat")) {
    stop(sample_id, " is not a Seurat object.")
  }
  if (!"RNA" %in% Assays(obj)) {
    stop(sample_id, " does not contain an RNA assay.")
  }

  DefaultAssay(obj) <- "RNA"
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  obj$sample_id <- sample_id
  obj$cancertype <- cancer_type
  obj$Response <- response

  n_before <- ncol(obj)
  obj <- subset(
    obj,
    subset = nFeature_RNA > min_features & percent.mt < max_percent_mt
  )
  message(sample_id, ": retained ", ncol(obj), "/", n_before, " cells after QC.")
  obj
}

obj_list <- Map(
  FUN = prepare_object,
  obj = obj_list,
  sample_id = manifest$sample_id,
  cancer_type = manifest$cancer_type,
  response = manifest$response
)

# SCTransform normalization and integration.
message("Running SCTransform...")
obj_list <- lapply(obj_list, function(x) {
  SCTransform(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(
  object.list = obj_list,
  nfeatures = nfeatures_integration
)

obj_list <- PrepSCTIntegration(
  object.list = obj_list,
  anchor.features = features,
  verbose = FALSE
)

anchors <- FindIntegrationAnchors(
  object.list = obj_list,
  normalization.method = "SCT",
  anchor.features = features,
  verbose = FALSE
)

integrated <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT",
  verbose = FALSE
)

DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = npcs, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = seq_len(npcs), verbose = FALSE)
integrated <- FindNeighbors(
  integrated,
  dims = seq_len(npcs),
  k.param = neighbor_k,
  verbose = FALSE
)
integrated <- FindClusters(
  integrated,
  resolution = cluster_resolution,
  verbose = FALSE
)

saveRDS(integrated, file.path(output_dir, "All_Resistant_Sensitive_integrated.rds"))

write.csv(
  as.data.frame(table(integrated$seurat_clusters, integrated$sample_id)),
  file.path(output_dir, "cluster_by_sample_counts.csv"),
  row.names = FALSE
)

pdf(file.path(output_dir, "UMAP_by_cancer_type.pdf"), width = 7, height = 6)
print(DimPlot(integrated, reduction = "umap", group.by = "cancertype", raster = FALSE))
dev.off()

pdf(file.path(output_dir, "UMAP_by_response.pdf"), width = 7, height = 6)
print(DimPlot(integrated, reduction = "umap", group.by = "Response", raster = FALSE))
dev.off()

message("Done. Integrated object written to: ", output_dir)
