

# Install seurat 4 to visualise the gene expression
remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

library(Seurat)

seurat_obj <- readRDS("HTA12_22_4_spatial_annotated.rds")
seurat_obj <- readRDS("HTA12_24_6_spatial_annotated.rds")
seurat_obj <- readRDS("GSM5708494_spatial_annotated.rds")
seurat_obj <- readRDS("GSE210616_094D_spatial_annotated_final.rds")

seurat_obj <- readRDS("Lung_spatial_data_annotated.rds") # Lung

seurat_obj$CellType

DefaultAssay(seurat_obj) <- "Spatial"

# UMAP plot
DimPlot(seurat_obj, reduction = "umap", group.by = "CellType", label = TRUE)

# Cluster overlay on tissue
SpatialDimPlot(seurat_obj, label = TRUE, group.by = "CellType", label.size = 3)

# Feature expression on tissue
SpatialFeaturePlot(seurat_obj, features = c("EPCAM", "CD3D", "COL1A1"), ncol = 3, 
                   image.alpha=0.8, pt.size.factor=2.5)




DefaultAssay(seurat_obj) <- "Spatial"


#DefaultAssay(brain) <- "Spatial"
#### Get extract mean expression of list of genes across cells
allo.smerged <- AddModuleScore(seurat_obj, features = Resistant.gene.signature, 
                               assay = "SCT", name = "Resistant.gene.signature", search = TRUE)

SpatialFeaturePlot(allo.smerged, features = "Resistant.gene.signature1", #alpha = c(0.1, 1),
                   image.alpha=0.8, pt.size.factor=1.5)

SpatialFeaturePlot(allo.smerged, features = "SRM", #alpha = c(0.1, 1),
                   image.alpha=0.8, pt.size.factor=2.5)


SpatialFeaturePlot(allo.smerged, features = "MYC", #alpha = c(0.1, 1),
                   image.alpha=0.8, pt.size.factor=2.5)


















