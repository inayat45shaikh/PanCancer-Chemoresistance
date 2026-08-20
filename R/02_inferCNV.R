#!/usr/bin/env Rscript

# 02_inferCNV.R
# inferCNV analysis of the integrated pan-cancer scRNA-seq dataset.
#
# The original script used "Normal epithelial cells" as the reference group and
# inferCNV parameters cutoff = 1, denoise = TRUE, HMM = TRUE. Those settings
# are retained here. Verify that cutoff = 1 is appropriate for the input
# technology represented by your processed objects before reproducing results.

suppressPackageStartupMessages({
  library(Seurat)
  library(infercnv)
  library(ggplot2)
})

set.seed(1234)

# ----------------------------- CONFIGURATION -----------------------------
input_rds <- "results/scRNAseq_integration/All_Resistant_Sensitive_integrated.rds"
gene_order_file <- "config/gencode_infercnv_format.txt"
output_dir <- "results/infercnv"

celltype_column <- "CellTypeNew"
reference_group <- "Normal epithelial cells"

infercnv_cutoff <- 1
use_hmm <- TRUE
use_denoise <- TRUE
# ------------------------------------------------------------------------

if (!file.exists(input_rds)) stop("Missing input Seurat object: ", input_rds)
if (!file.exists(gene_order_file)) stop("Missing gene-order file: ", gene_order_file)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) stop("Input RDS is not a Seurat object.")
if (!"RNA" %in% Assays(obj)) stop("Input object has no RNA assay.")
if (!celltype_column %in% colnames(obj[[]])) {
  stop("Metadata column '", celltype_column, "' is missing.")
}

DefaultAssay(obj) <- "RNA"

cell_annotation <- obj[[celltype_column, drop = TRUE]]
if (!reference_group %in% unique(cell_annotation)) {
  stop(
    "Reference group '", reference_group, "' was not found in ",
    celltype_column, "."
  )
}

# inferCNV accepts an annotation data.frame/file. Here we build a temporary
# two-column annotation table so cell names remain explicit and reproducible.
annotation_df <- data.frame(
  cell = colnames(obj),
  group = cell_annotation,
  stringsAsFactors = FALSE
)
annotation_file <- file.path(output_dir, "cell_annotations.tsv")
write.table(
  annotation_df,
  annotation_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

counts <- GetAssayData(obj, assay = "RNA", slot = "counts")

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file = annotation_file,
  delim = "\t",
  gene_order_file = gene_order_file,
  ref_group_names = reference_group
)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = infercnv_cutoff,
  out_dir = output_dir,
  cluster_by_groups = TRUE,
  denoise = use_denoise,
  HMM = use_hmm
)

# Reproduce the original cell-level deviation score:
# fraction of genes with normalized inferCNV values outside [0.95, 1.05].
cnv_score <- apply(
  infercnv_obj@expr.data,
  2,
  function(x) mean(x < 0.95 | x > 1.05, na.rm = TRUE)
)

score_df <- data.frame(
  cell = names(cnv_score),
  inferCNV_score = as.numeric(cnv_score),
  stringsAsFactors = FALSE
)

write.csv(
  score_df,
  file.path(output_dir, "infercnv_cell_scores.csv"),
  row.names = FALSE
)

write.csv(
  infercnv_obj@expr.data,
  file.path(output_dir, "infercnv_expression_matrix.csv")
)

p <- ggplot(score_df, aes(x = inferCNV_score)) +
  geom_histogram(bins = 50) +
  labs(x = "inferCNV score", y = "Number of cells") +
  theme_classic()

ggsave(
  file.path(output_dir, "infercnv_score_histogram.pdf"),
  p,
  width = 6,
  height = 4
)

message("Done. inferCNV outputs written to: ", output_dir)
