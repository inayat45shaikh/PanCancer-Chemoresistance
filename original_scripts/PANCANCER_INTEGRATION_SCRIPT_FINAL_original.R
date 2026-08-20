
install.packages("Seurat")

library(Seurat)
library(SingleR)
# Lung cancer Resistant
GSM4104144_SC16 <- readRDS("GSM4104144_SC16_clustered.rds")
GSM4104149_SC49 <- readRDS("GSM4104149_SC49_clustered.rds")
GSM4104155_SC55 <- readRDS("GSM4104155_SC55_clustered.rds")

GSM4104146_SC39 <- readRDS("GSM4104146_SC39_clustered.rds")
GSM4104154_SC4 <- readRDS("GSM4104154_SC4_clustered.rds")

SCLC1 <- readRDS("SCLC_Lung/SCLC1_Resistant_clustered.rds")
SCLC2 <- readRDS("SCLC_Lung/SCLC2_Resistant_clustered.rds")
SCLC5 <- readRDS("SCLC_Lung/SCLC5_Sensitive_clustered.rds")
SCLC7 <- readRDS("SCLC_Lung/SCLC7_Resistant_clustered.rds")

# TNBC cancer resistant and sensitive
HBCx_95_CAPAR_sensitive <- readRDS("HBCx_95_sensitive_clustered.rds")
HBCx_95_CAPAR_Resistant <- readRDS("HBCx_95_CAPAR_Resistant_clustered.rds")


# TNBC dataset kim et al, 2018

DimPlot(TNBC_Pre_Post_chemosensitive, reduction = "umap", group.by = "Drug")


# Pancreatic cancer
HTA12_1_1 <- readRDS("HTA12_1_1_clustered.rds")
HTA12_6 <- readRDS("HTA12_6_10_clustered.rds")

GSM6204114_P06 <- readRDS("GSE205013/GSM6204114_P06_clustered.rds")
GSM6204120_P12 <- readRDS("GSE205013/GSM6204120_P12_clustered.rds")

GSM6204111_P03 <- readRDS("GSE205013/GSM6204111_P03_clustered.rds")
GSM6204122_P14 <- readRDS("GSE205013/GSM6204122_P14_clustered.rds")

# Ovarian cancer Resistant
HGSOC_Ovarian_Resistant <- readRDS("HGSOC_Ovarian_Resistant_Clustered.rds")
HGSOC_Ovarian_Sensitive <- readRDS("HGSOC_Ovarian_Sensitive_Clustered.rds")





DefaultAssay(GSM4104144_SC16) <- "RNA"
DefaultAssay(GSM4104149_SC49) <- "RNA"
DefaultAssay(GSM4104155_SC55) <- "RNA"
DefaultAssay(GSM4104146_SC39) <- "RNA"
DefaultAssay(GSM4104154_SC4) <- "RNA"
DefaultAssay(SCLC1) <- "RNA"
DefaultAssay(SCLC2) <- "RNA"
DefaultAssay(SCLC5) <- "RNA"
DefaultAssay(SCLC7) <- "RNA"

DefaultAssay(HBCx_95_CAPAR_sensitive) <- "RNA"
DefaultAssay(HBCx_95_CAPAR_Resistant) <- "RNA"

DefaultAssay(TNBC_Pre_Post_chemoresistant) <- "RNA"
DefaultAssay(TNBC_Pre_Post_chemosensitive) <- "RNA"

DefaultAssay(BC159_T_3) <- "RNA"

DefaultAssay(HTA12_1_1) <- "RNA"
DefaultAssay(HTA12_6) <- "RNA"

DefaultAssay(GSM6204114_P06) <- "RNA"

DefaultAssay(GSM6204120_P12) <- "RNA"

DefaultAssay(GSM6204111_P03) <- "RNA"
DefaultAssay(GSM6204122_P14) <- "RNA"

DefaultAssay(HGSOC_Ovarian_Resistant) <- "RNA"
DefaultAssay(HGSOC_Ovarian_Sensitive) <- "RNA"




# merge all datasets, adding a cell ID to make sure cell names are unique



GSM4104144_SC16[["percent.mt"]]  <- PercentageFeatureSet(GSM4104144_SC16, pattern = "^MT-")
GSM4104149_SC49[["percent.mt"]]  <- PercentageFeatureSet(GSM4104149_SC49, pattern = "^MT-")

GSM4104155_SC55[["percent.mt"]]  <- PercentageFeatureSet(GSM4104155_SC55, pattern = "^MT-")

GSM4104146_SC39[["percent.mt"]]  <- PercentageFeatureSet(GSM4104146_SC39, pattern = "^MT-")

GSM4104154_SC4[["percent.mt"]]  <- PercentageFeatureSet(GSM4104154_SC4, pattern = "^MT-")

SCLC1[["percent.mt"]]  <- PercentageFeatureSet(SCLC1, pattern = "^MT-")
SCLC2[["percent.mt"]]  <- PercentageFeatureSet(SCLC2, pattern = "^MT-")
SCLC5[["percent.mt"]]  <- PercentageFeatureSet(SCLC5, pattern = "^MT-")
SCLC7[["percent.mt"]]  <- PercentageFeatureSet(SCLC7, pattern = "^MT-")

HBCx_95_CAPAR_sensitive[["percent.mt"]]  <- PercentageFeatureSet(HBCx_95_CAPAR_sensitive, pattern = "^MT-")
HBCx_95_CAPAR_Resistant[["percent.mt"]]  <- PercentageFeatureSet(HBCx_95_CAPAR_Resistant, pattern = "^MT-")

TNBC_Pre_Post_chemoresistant[["percent.mt"]]  <- PercentageFeatureSet(TNBC_Pre_Post_chemoresistant, pattern = "^MT-")
TNBC_Pre_Post_chemosensitive[["percent.mt"]]  <- PercentageFeatureSet(TNBC_Pre_Post_chemosensitive, pattern = "^MT-")

BC159_T_3[["percent.mt"]]  <- PercentageFeatureSet(BC159_T_3, pattern = "^MT-")

HTA12_1_1[["percent.mt"]]  <- PercentageFeatureSet(HTA12_1_1, pattern = "^MT-")
HTA12_6[["percent.mt"]]  <- PercentageFeatureSet(HTA12_6, pattern = "^MT-")

GSM6204114_P06[["percent.mt"]]  <- PercentageFeatureSet(GSM6204114_P06, pattern = "^MT-")
GSM6204120_P12[["percent.mt"]]  <- PercentageFeatureSet(GSM6204120_P12, pattern = "^MT-")

GSM6204111_P03[["percent.mt"]]  <- PercentageFeatureSet(GSM6204111_P03, pattern = "^MT-")
GSM6204122_P14[["percent.mt"]]  <- PercentageFeatureSet(GSM6204122_P14, pattern = "^MT-")


HGSOC_Ovarian_Resistant[["percent.mt"]]  <- PercentageFeatureSet(HGSOC_Ovarian_Resistant, pattern = "^MT-")
HGSOC_Ovarian_Sensitive[["percent.mt"]]  <- PercentageFeatureSet(HGSOC_Ovarian_Sensitive, pattern = "^MT-")






GSM4104144_SC16 <- subset(GSM4104144_SC16, subset = nFeature_RNA > 200 & percent.mt < 25)
GSM4104149_SC49 <- subset(GSM4104149_SC49, subset = nFeature_RNA > 200 & percent.mt < 25)

GSM4104155_SC55 <- subset(GSM4104155_SC55, subset = nFeature_RNA > 200 & percent.mt < 25)

GSM4104146_SC39 <- subset(GSM4104146_SC39, subset = nFeature_RNA > 200 & percent.mt < 25)

GSM4104154_SC4 <- subset(GSM4104154_SC4, subset = nFeature_RNA > 200 & percent.mt < 25)

SCLC1 <- subset(SCLC1, subset = nFeature_RNA > 200 & percent.mt < 25)
SCLC2 <- subset(SCLC2, subset = nFeature_RNA > 200 & percent.mt < 25)
SCLC5 <- subset(SCLC5, subset = nFeature_RNA > 200 & percent.mt < 25)
SCLC7 <- subset(SCLC7, subset = nFeature_RNA > 200 & percent.mt < 25)

HBCx_95_CAPAR_sensitive <- subset(HBCx_95_CAPAR_sensitive, subset = nFeature_RNA > 200 & percent.mt < 25)
HBCx_95_CAPAR_Resistant <- subset(HBCx_95_CAPAR_Resistant, subset = nFeature_RNA > 200 & percent.mt < 25)

TNBC_Pre_Post_chemoresistant <- subset(TNBC_Pre_Post_chemoresistant, subset = nFeature_RNA > 200 & percent.mt < 25)
TNBC_Pre_Post_chemosensitive <- subset(TNBC_Pre_Post_chemosensitive, subset = nFeature_RNA > 200 & percent.mt < 25)

BC159_T_3 <- subset(BC159_T_3, subset = nFeature_RNA > 200 & percent.mt < 25)

HTA12_1_1 <- subset(HTA12_1_1, subset = nFeature_RNA > 200 & percent.mt < 25)
HTA12_6 <- subset(HTA12_6, subset = nFeature_RNA > 200 & percent.mt < 25)

GSM6204114_P06 <- subset(GSM6204114_P06, subset = nFeature_RNA > 200 & percent.mt < 25)
GSM6204120_P12 <- subset(GSM6204120_P12, subset = nFeature_RNA > 200 & percent.mt < 25)

GSM6204111_P03 <- subset(GSM6204111_P03, subset = nFeature_RNA > 200 & percent.mt < 25)
GSM6204122_P14 <- subset(GSM6204122_P14, subset = nFeature_RNA > 200 & percent.mt < 25)

HGSOC_Ovarian_Resistant <- subset(HGSOC_Ovarian_Resistant, subset = nFeature_RNA > 200 & percent.mt < 25)
HGSOC_Ovarian_Sensitive <- subset(HGSOC_Ovarian_Sensitive, subset = nFeature_RNA > 200 & percent.mt < 25)





pbmc_list <- list()
pbmc_list[["GSM4104144_SC16"]] <- GSM4104144_SC16
pbmc_list[["GSM4104149_SC49"]] <- GSM4104149_SC49


pbmc_list[["GSM4104155_SC55"]] <- GSM4104155_SC55

pbmc_list[["GSM4104146_SC39"]] <- GSM4104146_SC39


pbmc_list[["GSM4104154_SC4"]] <- GSM4104154_SC4

pbmc_list[["SCLC1"]] <- SCLC1
pbmc_list[["SCLC2"]] <- SCLC2
pbmc_list[["SCLC5"]] <- SCLC5
pbmc_list[["SCLC7"]] <- SCLC7

pbmc_list[["HBCx_95_CAPAR_sensitive"]] <- HBCx_95_CAPAR_sensitive
pbmc_list[["HBCx_95_CAPAR_Resistant"]] <- HBCx_95_CAPAR_Resistant

pbmc_list[["TNBC_Pre_Post_chemoresistant"]] <- TNBC_Pre_Post_chemoresistant
pbmc_list[["TNBC_Pre_Post_chemosensitive"]] <- TNBC_Pre_Post_chemosensitive

pbmc_list[["BC159_T_3"]] <- BC159_T_3

pbmc_list[["HTA12_1_1"]] <- HTA12_1_1
pbmc_list[["HTA12_6"]] <- HTA12_6

pbmc_list[["GSM6204114_P06"]] <- GSM6204114_P06

pbmc_list[["GSM6204120_P12"]] <- GSM6204120_P12

pbmc_list[["GSM6204111_P03"]] <- GSM6204111_P03
pbmc_list[["GSM6204122_P14"]] <- GSM6204122_P14


pbmc_list[["HGSOC_Ovarian_Resistant"]] <- HGSOC_Ovarian_Resistant
pbmc_list[["HGSOC_Ovarian_Sensitive"]] <- HGSOC_Ovarian_Sensitive



for (i in 1:length(pbmc_list)) {
  pbmc_list[[i]] <- NormalizeData(pbmc_list[[i]], verbose = F)
  pbmc_list[[i]] <- FindVariableFeatures(pbmc_list[[i]], selection.method = "vst", nfeatures = 3000, verbose = F)
  pbmc_list[[i]] <- ScaleData(pbmc_list[[i]])
  pbmc_list[[i]] <- RunPCA(pbmc_list[[i]])
  pbmc_list[[i]] <- FindNeighbors(pbmc_list[[i]], dims = 1:30, reduction = "pca")
  pbmc_list[[i]] <- FindClusters(pbmc_list[[i]], resolution = 2, cluster.name = "unintegrated_clusters")
  
}




patient.list <- pbmc_list

#BiocManager::install('glmGamPoi')

for (i in 1:length(patient.list)) {
  patient.list[[i]] <- SCTransform(patient.list[[i]], verbose = FALSE)
  
}
patient.features <- SelectIntegrationFeatures(object.list = patient.list, nfeatures = 3000)
patient.list <- PrepSCTIntegration(object.list = patient.list, anchor.features = patient.features, 
                                   verbose = FALSE)
patient.anchors <- FindIntegrationAnchors(object.list = patient.list, normalization.method = "SCT", 
                                          anchor.features = patient.features, verbose = FALSE)

patient.integrated <- IntegrateData(anchorset = patient.anchors, normalization.method = "SCT", 
                                    verbose = FALSE)







features <- SelectIntegrationFeatures(object.list = pbmc_list, nfeatures=3000)

donor.list <- lapply(X = pbmc_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

patient.list <- PrepSCTIntegration(object.list = donor.list, anchor.features = patient.features, 
                                   verbose = FALSE)


anchors <- FindIntegrationAnchors(object.list = donor.list, #k.anchor=20,
                                  normalization.method = "SCT", # SCT
                                  #reference = c(26, 62, 66, 67, 68, 73), 
                                  reduction = "rpca", anchor.features = features)

pbmc_seurat  <- IntegrateData(anchorset = anchors, dims = 1:30)





DefaultAssay(pbmc_seurat) <- "integrated"
pbmc_seurat <- ScaleData(pbmc_seurat, verbose = F)
pbmc_seurat <- RunPCA(pbmc_seurat, npcs = 30, verbose = F)
pbmc_seurat <- RunUMAP(pbmc_seurat, reduction = "pca", dims = 1:30, verbose = F)

DimPlot(pbmc_seurat, reduction = "umap", split.by = "cancertype", raster=FALSE) #+ plot_annotation(title = "after integration (Seurat 3)")

pbmc_seurat <- FindNeighbors(pbmc_seurat, dims = 1:30, k.param = 10, verbose = F)
pbmc_seurat <- FindClusters(pbmc_seurat, verbose = F, resolution = 0.5)
DimPlot(pbmc_seurat,label = F, group.by = "cancertype", raster=FALSE) #+ NoLegend()

count_table <- table(pbmc_seurat@meta.data$seurat_clusters, pbmc_seurat@meta.data$orig.ident)
count_table
plot_integrated_clusters(pbmc_seurat) 



DimPlot(pbmc_seurat, reduction = "umap", group.by = "cancertype", 
        #split.by = "Response", 
        raster=FALSE)

saveRDS(pbmc_seurat, file = "All_Resist_Sensitive_integrated.rds")







chemo_resistance_genes <- read.delim("chemo_resistance_genes.txt")
CIN70.signature <- read.delim("CIN70.signature.txt")



DimPlot(seurat, reduction = "umap", group.by = "Cellstate", label = T, split.by = "Response")

DimPlot(seurat, reduction = "umap", group.by = "CellType", label = T)

DefaultAssay(seurat) <- "RNA"
#### Get extract mean expression of list of genes across cells
allo.smerged <- AddModuleScore(seurat, features = chemo_resistance_genes, 
                               assay = "RNA", name = "chemo_resistance_genes", search = TRUE)

#### Get extract mean expression of list of genes across cells
allo.smerged <- AddModuleScore(seurat, features = CIN70.signature, 
                               assay = "RNA", name = "CIN70.signature", search = TRUE)
							   
							   

library(ggplot2)
install.packages("ggpubr")
library(ggpubr)
library(tidyverse)

median.stat <- function(x){
  out <- quantile(x, probs = c(0.5))
  names(out) <- c("ymed")
  return(out) 
}


VlnPlot(object = allo.smerged, features = "chemo_resistance_genes1", group.by = "Response", 
        split.by = "Response", pt.size = 0)+ NoLegend() +
  geom_boxplot(width=0.2, fill="white") +
  stat_compare_means()

VlnPlot(object = allo.smerged, features = "CIN70.signature1", group.by = "Response", 
        split.by = "Response", pt.size = 0)+ NoLegend() +
  geom_boxplot(width=0.2, fill="white") +
  stat_compare_means()




VlnPlot(object = seurat, features = c("SRM","MYC"), 
        group.by = "ResistantScoreCategory", 
        #split.by = "Response", 
        pt.size = 0)
		