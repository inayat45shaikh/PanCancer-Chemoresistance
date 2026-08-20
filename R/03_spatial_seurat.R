#!/usr/bin/env Rscript

# 03_spatial_seurat.R
# Seurat-based visualization and resistance-signature scoring for spatial data.

suppressPackageStartupMessages({
  library(Seurat)
})

set.seed(1234)

# ----------------------------- CONFIGURATION -----------------------------
manifest_file <- "config/spatial_manifest_template.csv"
signature_file <- "config/resistant_signature.txt"
output_dir <- "results/spatial_seurat"

celltype_column <- "CellType"
spatial_assay <- "Spatial"
scoring_assay <- "SCT"
# ------------------------------------------------------------------------

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(manifest_file)) stop("Missing spatial manifest: ", manifest_file)
if (!file.exists(signature_file)) {
  stop(
    "Missing resistance signature file: ", signature_file,
    "\nProvide one gene symbol per line."
  )
}

manifest <- read.csv(manifest_file, stringsAsFactors = FALSE)
required_cols <- c("sample_id", "path", "cancer_type")
missing_cols <- setdiff(required_cols, colnames(manifest))
if (length(missing_cols) > 0) {
  stop("Spatial manifest is missing: ", paste(missing_cols, collapse = ", "))
}

resistant_signature <- scan(signature_file, what = character(), quiet = TRUE)
resistant_signature <- unique(resistant_signature[nzchar(resistant_signature)])

for (i in seq_len(nrow(manifest))) {
  sample_id <- manifest$sample_id[i]
  input_file <- manifest$path[i]

  if (!file.exists(input_file)) {
    warning("Skipping ", sample_id, "; file not found: ", input_file)
    next
  }

  message("Processing spatial sample: ", sample_id)
  obj <- readRDS(input_file)

  if (!inherits(obj, "Seurat")) {
    warning("Skipping ", sample_id, "; object is not a Seurat object.")
    next
  }

  sample_out <- file.path(output_dir, sample_id)
  dir.create(sample_out, recursive = TRUE, showWarnings = FALSE)

  if (spatial_assay %in% Assays(obj)) {
    DefaultAssay(obj) <- spatial_assay
  }

  if (celltype_column %in% colnames(obj[[]]) && "umap" %in% Reductions(obj)) {
    pdf(file.path(sample_out, "UMAP_celltype.pdf"), width = 7, height = 6)
    print(DimPlot(obj, reduction = "umap", group.by = celltype_column, label = TRUE))
    dev.off()
  }

  if (celltype_column %in% colnames(obj[[]]) && length(Images(obj)) > 0) {
    pdf(file.path(sample_out, "Spatial_celltype.pdf"), width = 8, height = 7)
    print(SpatialDimPlot(obj, label = TRUE, group.by = celltype_column, label.size = 3))
    dev.off()
  }

  marker_genes <- intersect(c("EPCAM", "CD3D", "COL1A1"), rownames(obj))
  if (length(marker_genes) > 0 && length(Images(obj)) > 0) {
    pdf(file.path(sample_out, "Spatial_compartment_markers.pdf"), width = 12, height = 4)
    print(
      SpatialFeaturePlot(
        obj,
        features = marker_genes,
        ncol = length(marker_genes),
        image.alpha = 0.8,
        pt.size.factor = 2.5
      )
    )
    dev.off()
  }

  if (!scoring_assay %in% Assays(obj)) {
    warning(
      sample_id, ": assay '", scoring_assay,
      "' not found; resistance signature score was not calculated."
    )
  } else {
    obj <- AddModuleScore(
      obj,
      features = list(resistant_signature),
      assay = scoring_assay,
      name = "ResistantSignature",
      search = TRUE
    )

    if (length(Images(obj)) > 0) {
      pdf(file.path(sample_out, "Spatial_resistance_signature.pdf"), width = 7, height = 6)
      print(
        SpatialFeaturePlot(
          obj,
          features = "ResistantSignature1",
          image.alpha = 0.8,
          pt.size.factor = 1.5
        )
      )
      dev.off()
    }
  }

  for (gene in intersect(c("SRM", "MYC"), rownames(obj))) {
    if (length(Images(obj)) > 0) {
      pdf(file.path(sample_out, paste0("Spatial_", gene, ".pdf")), width = 7, height = 6)
      print(
        SpatialFeaturePlot(
          obj,
          features = gene,
          image.alpha = 0.8,
          pt.size.factor = 2.5
        )
      )
      dev.off()
    }
  }

  saveRDS(obj, file.path(sample_out, paste0(sample_id, "_scored.rds")))
}

message("Done. Spatial outputs written to: ", output_dir)
