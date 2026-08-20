#!/usr/bin/env Rscript

# 04_spacet_deconvolution.R
# SpaCET deconvolution using matched scRNA-seq resistance-state references.

suppressPackageStartupMessages({
  library(Seurat)
  library(SpaCET)
})

set.seed(1234)

# ----------------------------- CONFIGURATION -----------------------------
spacet_rds <- "data/spatial/example_SpaCET.rds"
reference_rds <- "results/reference/Resistant_only_malignant.rds"
output_dir <- "results/spacet"

resistance_state_column <- "ResistantScoreCategory"
core_no <- 56

# Optional external GMT file. Leave as NA_character_ if unused.
external_gmt <- NA_character_
# ------------------------------------------------------------------------

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(spacet_rds)) stop("Missing SpaCET object: ", spacet_rds)
if (!file.exists(reference_rds)) stop("Missing scRNA-seq reference: ", reference_rds)

spacet_obj <- readRDS(spacet_rds)
reference_obj <- readRDS(reference_rds)

if (!inherits(reference_obj, "Seurat")) {
  stop("reference_rds must contain a Seurat object.")
}
if (!resistance_state_column %in% colnames(reference_obj[[]])) {
  stop("Missing metadata column: ", resistance_state_column)
}

DefaultAssay(reference_obj) <- "RNA"
sc_counts <- as.matrix(GetAssayData(reference_obj, assay = "RNA", slot = "counts"))

sc_annotation <- reference_obj[[resistance_state_column, drop = TRUE]]
names(sc_annotation) <- colnames(reference_obj)

if (anyNA(sc_annotation)) stop("Resistance-state annotation contains NA values.")

# Flat lineage tree, matching the original analysis.
states <- unique(as.character(sc_annotation))
sc_lineage_tree <- as.list(states)
names(sc_lineage_tree) <- states

spacet_obj <- SpaCET.deconvolution.matched.scRNAseq(
  spacet_obj,
  sc_counts = sc_counts,
  sc_annotation = sc_annotation,
  sc_lineageTree = sc_lineage_tree,
  coreNo = core_no
)

pdf(file.path(output_dir, "Spatial_resistance_states.pdf"), width = 9, height = 7)
print(
  SpaCET.visualize.spatialFeature(
    spacet_obj,
    spatialType = "CellFraction",
    spatialFeatures = intersect(c("Low", "Medium", "High"), states),
    nrow = 2
  )
)
dev.off()

pdf(file.path(output_dir, "Spatial_SRM_MYC.pdf"), width = 9, height = 5)
print(
  SpaCET.visualize.spatialFeature(
    spacet_obj,
    spatialType = "GeneExpression",
    spatialFeatures = c("SRM", "MYC"),
    nrow = 2
  )
)
dev.off()

# Built-in SpaCET gene-set collections used in the original analysis.
for (gene_set_collection in c("Hallmark", "CancerCellState", "TLS")) {
  message("Calculating SpaCET gene-set scores: ", gene_set_collection)
  spacet_obj <- SpaCET.GeneSetScore(
    spacet_obj,
    GeneSets = gene_set_collection
  )
}

# Original custom High and Sensitive gene signatures.
custom_signatures <- list(
  High = c(
    "FABP5","STMN1","SCPEP1","HES6","AVP","RPS10","CPE","HMGN2","SEPT11",
    "PSIP1","TUBA1B","SEPT7","TUBA1A","TUBB4B","HIST1H4C","SOX2","NDUFA4L2",
    "CXCL17","PCNA","FAM96B","HMGB2","TYMS","GNB2L1","NEB","SNRPB","ASCL1",
    "TOP2A","MINOS1","COTL1","NUSAP1","PRDX2","PGF","SEC11C","CAMK2B","LMO2",
    "C8orf59","MRPS26","PTTG1","INSM1","C21orf59","RGS16","BRIX1","NKX2-1",
    "CDKN2C","TUBB","TMEM97","PLIN2","SIX1","KIF19","CA9","PAFAH1B3","EIF4A3",
    "AQP3","CENPV","TUBB2B","NOP56","AES","FAM111B","SOX4","BNIP3","SNRPB2",
    "TCF4","ILF2","CCT5","ENO1","USP1","C19orf70","UBE2T","NME1","ALDOA",
    "HMGB3","H2AFZ","NASP","GINS2","DNMT1","FABP7","SNRPG","NFIB","MYC","UBA2",
    "TUBA1C","HNRNPM","IRX2","CYB5A","SRM","COPS3","UQCC2","MAFB","ASF1B",
    "MTRNR2L2","MARCKSL1","GGCT","DDX39A","DYNLT1","SFPQ","SCG3","C6orf48",
    "CKB","GMNN","ACYP1"
  ),
  Sens = c(
    "NPC2","AHNAK2","CLDN1","ANXA4","RIOK3","MAL2","C15orf48","BLVRB","CDA",
    "ITGB6","PLAAT3","S100A6","TIMP2","SPTBN1","BIRC3","LY6D","S100A14",
    "ISG20","MET","IGFBP7","RPS17","SYT8","MYL12A","RHOC","MSLN","GPRC5A",
    "S100A16","NRP2","SELENOW","C4orf48","CD47","UBE2H","RNASEK","BHLHE40",
    "LMO7","GAS5","ITGB1","TFPI","AREG","CRIP1","F3","ERO1A","MT-ND4L",
    "PLAAT4","CD59","DHRS9","FXYD3","TSPO","RRAS","TXNIP","S100P","TMEM265",
    "CD9","APLP2","ITPRID2","LAMC2","RAB11FIP1","SDC4","CAV2","CLIC3","TIMP1",
    "NME2","EEF1G","MMP7","SEPTIN2","DSC2","CAST","POLD4","S100A11","TSPAN8",
    "ANXA1","CEACAM6","MUCL3","CAV1","KRT17","SEPTIN7","SNHG29","ATP6V0C",
    "FAM3C","IGFBP6","LAMA3","SLPI","AHNAK","KRT19","CIB1","EMP3","LCN2",
    "LGALS3","COL17A1","G0S2","ITGA6","DSG2","LINC01578","DCBLD2","C19orf33",
    "LAMB3","TM4SF1","SERPINA3","S100A10","S100A4"
  )
)

spacet_obj <- SpaCET.GeneSetScore(spacet_obj, GeneSets = custom_signatures)

pdf(file.path(output_dir, "Spatial_custom_resistance_signatures.pdf"), width = 9, height = 5)
print(
  SpaCET.visualize.spatialFeature(
    spacet_obj,
    spatialType = "GeneSetScore",
    spatialFeatures = c("High", "Sens")
  )
)
dev.off()

if (!is.na(external_gmt)) {
  if (!file.exists(external_gmt)) stop("External GMT file not found: ", external_gmt)
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("clusterProfiler is required to read an external GMT file.")
  }
  gmt_table <- clusterProfiler::read.gmt(external_gmt)
  gmt_list <- split(gmt_table$gene, gmt_table$term)
  spacet_obj <- SpaCET.GeneSetScore(spacet_obj, GeneSets = gmt_list)
}

saveRDS(spacet_obj, file.path(output_dir, "SpaCET_scored.rds"))
message("Done. SpaCET outputs written to: ", output_dir)
