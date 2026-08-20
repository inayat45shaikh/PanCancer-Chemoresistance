
install.packages("Seurat")

# install.packages("devtools")
devtools::install_github("data2intelligence/SpaCET")

library(Seurat)
library(SpaCET)

SpaCET_obj <- readRDS("HTA12_22_4_spatial_SpaCET.rds") # Pancreatic

SpaCET_obj <- readRDS("94D_SpaCET.rds") # TNBC

SpaCET_obj <- readRDS("GSM5708493_SpaCET.rds") # Ovarian

SpaCET_obj <- readRDS("GSM5708494_SpaCET.rds") # Ovarian

seurat_obj <- readRDS("cleaned_Resistant_only_malignant_2025_ResistantScore_Pseudotime_MSigDB_19052025.rds")
seurat_obj <- readRDS("cleaned_Resistant_only_malignant_2025_ResistantScore_pseudotime_16052025_.rds")

# Extract raw count matrix from the epithelial cell subset
sc_counts <- as.matrix(Seurat::GetAssayData(seurat_obj, layer = "counts"))


sc_annotation <- seurat_obj@meta.data$ResistantScoreCategory
names(sc_annotation) <- colnames(seurat_obj)



table(sc_annotation)

unique_celltypes <- unique(sc_annotation)

# Flat lineage tree (if unsure about hierarchy)
sc_lineageTree <- as.list(unique_celltypes)
names(sc_lineageTree) <- unique_celltypes


sc_annotation <- read.delim("annotation.txt", sep = "\t")



SpaCET_obj <- SpaCET.deconvolution.matched.scRNAseq(
  SpaCET_obj, 
  sc_counts=sc_counts, 
  sc_annotation=sc_annotation, 
  sc_lineageTree=sc_lineageTree, 
  coreNo=56
)

SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "CellFraction",
  spatialFeatures = c("Low","Medium","High"),
  nrow=2
)



# Markers for cancer clone A and B, acinar cell, and centroacinar like ductal cell
SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "GeneExpression",
  spatialFeatures = c("SRM","MYC"),
  nrow=2
)


######### Calculate hallmark score
# run gene set calculation
SpaCET_obj <- SpaCET.GeneSetScore(SpaCET_obj, GeneSets="Hallmark")

# show results
SpaCET_obj@results$GeneSetScore[1:6,1:6]

# show all gene sets
rownames(SpaCET_obj@results$GeneSetScore)

# visualize two gene sets
SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "GeneSetScore", 
  spatialFeatures = c("HALLMARK_WNT_BETA_CATENIN_SIGNALING","HALLMARK_TGF_BETA_SIGNALING")
)



######### Calculate cancer cell state score
# run gene set calculation
SpaCET_obj <- SpaCET.GeneSetScore(SpaCET_obj, GeneSets="CancerCellState")

# show all gene sets
rownames(SpaCET_obj@results$GeneSetScore)

# visualize two gene sets
SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "GeneSetScore", 
  spatialFeatures = c("CancerCellState_Cycle","CancerCellState_cEMT")
)


###### Calculate TLS score
# run gene set calculation
SpaCET_obj <- SpaCET.GeneSetScore(SpaCET_obj, GeneSets="TLS")

# visualize TLS
SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "GeneSetScore", 
  spatialFeatures = c("TLS")
)



######## Calculate other gene sets’ score
# 1)
gmt1 <- list(
  High = c("FABP5","STMN1","SCPEP1","HES6","AVP","RPS10","CPE","HMGN2","SEPT11","PSIP1","TUBA1B","SEPT7","TUBA1A","TUBB4B","HIST1H4C","SOX2","NDUFA4L2","CXCL17","PCNA","FAM96B","HMGB2","TYMS","GNB2L1","NEB","SNRPB","ASCL1","TOP2A","MINOS1","COTL1","NUSAP1","PRDX2","PGF","SEC11C","CAMK2B","LMO2","C8orf59","MRPS26","PTTG1","INSM1","C21orf59","RGS16","BRIX1","NKX2-1","CDKN2C","TUBB","TMEM97","PLIN2","SIX1","KIF19","CA9","PAFAH1B3","EIF4A3","AQP3","CENPV","TUBB2B","NOP56","AES","FAM111B","SOX4","BNIP3","SNRPB2","TCF4","ILF2","CCT5","ENO1","USP1","C19orf70","UBE2T","NME1","ALDOA","HMGB3","H2AFZ","NASP","GINS2","DNMT1","FABP7","SNRPG","NFIB","MYC","UBA2","TUBA1C","HNRNPM","IRX2","CYB5A","SRM","COPS3","UQCC2","MAFB","ASF1B","MTRNR2L2","MARCKSL1","GGCT","DDX39A","DYNLT1","SFPQ","SCG3","C6orf48","CKB","GMNN","ACYP1"), 
  Sens = c("NPC2","AHNAK2","CLDN1","ANXA4","RIOK3","MAL2","C15orf48","BLVRB","CDA","ITGB6","PLAAT3","S100A6","TIMP2","SPTBN1","BIRC3","LY6D","S100A14","ISG20","MET","IGFBP7","RPS17","SYT8","MYL12A","RHOC","MSLN","GPRC5A","S100A16","NRP2","SELENOW","C4orf48","CD47","UBE2H","RNASEK","BHLHE40","LMO7","GAS5","ITGB1","TFPI","AREG","CRIP1","F3","ERO1A","MT-ND4L","PLAAT4","CD59","DHRS9","FXYD3","TSPO","RRAS","TXNIP","S100P","TMEM265","CD9","APLP2","ITPRID2","LAMC2","RAB11FIP1","SDC4","CAV2","CLIC3","TIMP1","NME2","EEF1G","MMP7","SEPTIN2","DSC2","CAST","POLD4","S100A11","TSPAN8","ANXA1","CEACAM6","MUCL3","CAV1","KRT17","SEPTIN7","SNHG29","ATP6V0C","FAM3C","IGFBP6","LAMA3","SLPI","AHNAK","KRT19","CIB1","EMP3","LCN2","LGALS3","COL17A1","G0S2","ITGA6","DSG2","LINC01578","DCBLD2","C19orf33","LAMB3","TM4SF1","SERPINA3","S100A10","S100A4")
)
SpaCET_obj <- SpaCET.GeneSetScore(SpaCET_obj, GeneSets = gmt1)

# 2)
gmt2 <- read.gmt("Path_to_gmt_file")
SpaCET_obj <- SpaCET.GeneSetScore(SpaCET_obj, GeneSets = gmt2)


rownames(SpaCET_obj@results$GeneSetScore)


# visualize two gene sets
SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "GeneSetScore", 
  spatialFeatures = c("High","Sens")
)



