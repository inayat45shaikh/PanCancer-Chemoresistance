
library(future)
parallelly::availableCores()
plan(multisession, workers = 30)
options(future.globals.maxSize = 350000 * 1024^2)

install.packages("stringr")                        # Install stringr package
library("stringr")   

#sudo apt-get install jags
install.packages("rjags")

BiocManager::install("SummarizedExperiment")

library("devtools")
devtools::install_github("broadinstitute/infercnv")

install.packages("Seurat")

library(Seurat)
library(infercnv)

output_dir = "output_PANCANCER"

options(Seurat.object.assay.version = "v3")


seuratobj<- readRDS("All_Resist_Sensitive_integrated.rds")




# convert a v5 assay to a v3 assay
seuratobj[["RNA"]] <- as(object = seuratobj[["RNA"]], Class = "Assay")

# convert a v3 assay to a v5 assay
pbmc3k[["RNA5"]] <- as(object = pbmc3k[["RNA3"]], Class = "Assay5")

count_mtx <- seuratobj@assays$RNA@counts


Idents(seuratobj) <- "CellTypeNew"

table(seuratobj@active.ident)

gene_order <- read.table("gencode_infercnv_formast.txt", row.names=1)

options(scipen = 100)

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=GetAssayData(seuratobj, "RNA"),
                                    annotations_file=as.matrix(seuratobj@active.ident),
                                    delim="\t",
                                    gene_order_file=gene_order,
                                    ref_group_names=c("Normal epithelial cells"))


# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="output_PANCANCER",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T
)


table(infercnv_obj@tumor_subclusters$subclusters)

scores=apply(infercnv_obj@expr.data,2,function(x){ sum(x < 0.95 | x > 1.05)/length(x) })

write.csv(scores, file="output/infercnv_obj_score.csv")

write.csv(infercnv_obj@expr.data, file="output_epith_vs_Normal/infercnv_obj_expression.csv")


hist(scores, xlab = "InferCNV score",
     col = "grey", border = "black")


abline(v=quantile(scores, .2), col='red', lwd=3, lty="dotted")

dat <- data.frame( x=scores, above=scores>.2 )
library(ggplot2)
qplot(x,data=dat,geom="histogram",fill=above)



ggplot(scores, aes(x = InferCNV)) +
  geom_histogram() +  facet_grid(vars(condition))


ggplot(scores, aes(x = InferCNV_Score)) +
  geom_histogram() +  facet_grid(vars(CellType))


ggplot(scores,
       aes(x = `InferCNV_Score`, fill = CellType))+
  geom_histogram(binwidth=1) +  facet_grid(vars(CellType))+
  geom_vline(data = file, aes(xintercept = .2),  color="red", linetype="dashed")+  theme_bw()
