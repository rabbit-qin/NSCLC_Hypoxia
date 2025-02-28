rm(list = ls())
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(readr)
library(sctransform)
library(ggsci)
library(paletteer)
library(glmGamPoi)

setwd("./NSCLC_SingleCell/")

experiment.data <- Read10X("./NSCLC_SingleCell/LUNG_jy8/")

experiment.aggregate <- CreateSeuratObject( experiment.data,
                                            project = "NSCLC", 
                                            min.cells = 3,
                                            min.features = 200)

dim(experiment.aggregate@meta.data)
View(experiment.aggregate@meta.data)

samInfo <- rep("Bone_Marrow_1",nrow(experiment.aggregate@meta.data))
Qexperiment.aggregate <- AddMetaData( object = experiment.aggregate,
                                      metadata = samInfo,
                                      col.name = "samInfo")
View(experiment.aggregate@meta.data)

experiment.aggregate[["percent.mt"]] <- PercentageFeatureSet( experiment.aggregate, 
                                                              pattern = "^MT-")

VlnPlot( experiment.aggregate, 
         features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
         ncol = 3)

gene.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$nFeature_RNA,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
rna.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$nCount_RNA,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
mt.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$percent.mt,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
freq.combine <- as.data.frame(cbind(gene.freq,rna.freq,mt.freq))
colnames(freq.combine) <- c("Count_Gene","Count_RNA","MT_percent")
View(freq.combine)
rm(gene.freq,rna.freq,mt.freq)

plot1 <- FeatureScatter(experiment.aggregate, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(experiment.aggregate, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2),legend = "none")
rm(plot1,plot2)

cat("Before filter :",nrow(experiment.aggregate@meta.data),"cells\n")
experiment.aggregate <- subset( experiment.aggregate, 
                                subset = nFeature_RNA > 101 & 
                                  nFeature_RNA < 6000 &
                                  nCount_RNA > 200 & 
                                  #nCount_RNA < 20000 &
                                  percent.mt < 10)
cat("After filter :",nrow(experiment.aggregate@meta.data),"cells\n")

experiment.aggregate <- SCTransform(experiment.aggregate, vars.to.regress = "percent.mt", verbose = FALSE)
experiment.aggregate <- RunPCA(experiment.aggregate, verbose = FALSE) 
DimPlot(experiment.aggregate,reduction = "pca")+NoLegend()
ElbowPlot(experiment.aggregate, ndims = 50, reduction = "pca")
experiment.aggregate <- RunTSNE(experiment.aggregate, dims = 1:20, verbose = FALSE)  

experiment.aggregate <- FindNeighbors(experiment.aggregate, dims = 1:20, verbose = FALSE) 
experiment.aggregate <- FindClusters(experiment.aggregate, verbose = FALSE) 

p1 <- DimPlot(experiment.aggregate, label = T)
p1 <- DotPlot(experiment.aggregate, features = c("CLDN18","FOLR1","AQP4","PEBP4","CLDN5","FLT1","CDH5","RAMP2","CAPS","TMEM190",
                                                 "PIFO","SNTN","COL1A1","DCN","COL1A2","C1R","CD79A","IGKC","IGLC3","IGHG3",
                                                 "CD3D","TRBC1","TRBC2","TRAC","LYZ","MARCO","CD68","FCGR3A","EPCAM"),cols = c("blue", "red"))+
  theme_classic()

marker <- c("COL1A2")
p0<-FeaturePlot(experiment.aggregate,cols = c("lightgrey", "red"),features = marker)+
  theme_bw() +
  theme(legend.position = 'none') +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
p0
ggsave("./NSCLC_SingleCell/Figure/BIOMARKERS/COL1A2.png",p0,width = 5,height = 5,device = "png",dpi = 1200,units = "in")

new.cluster.ids <- c("T_cell","T_cell","T_cell","B_cell","Myeloid","T_cell","T_cell","Cancer","Cancer","Myeloid","Myeloid","Myeloid",
                     "Cancer","T_cell","B_cell","Fibroblast","T_cell","Cancer","Endothelial","Alveolar","Unknown","Myeloid","T_cell",
                     "Cancer","T_cell","B_cell","Epithelial")
names(new.cluster.ids) <- levels(experiment.aggregate)
experiment.aggregate <- RenameIdents(experiment.aggregate, new.cluster.ids)
experiment.aggregate@meta.data$ident <- experiment.aggregate@active.ident
pal <- paletteer_d( "ggsci::nrc_npg")[c(1,3,4,7,5,2,6,8)]
p1 <- DimPlot(experiment.aggregate, label = F,pt.size = 0.6, repel = T) +
  theme_bw() +
  theme(legend.position = 'none') +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
p1
ggsave("./NSCLC_SingleCell/Figure/Clusters.png",p1,width = 5,height = 5,device = "png",dpi = 1200,units = "in")

p3 <- DimPlot(experiment.aggregate, group.by = "ident", label = F,pt.size = 0.6, repel = T) +
  scale_color_manual(values = c("Epithelial" = "blue")) +
  theme_bw() +
  theme(legend.position = 'none') +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
p3
ggsave("./NSCLC_SingleCell/Figure/subclusters/Epithelial.png",p3,width = 5,height = 5,device = "png",dpi = 1200,units = "in")

experiment.aggregate <- subset(experiment.aggregate, idents = "Unknown", invert = TRUE)
expression_matrix <- GetAssayData(experiment.aggregate, slot = "counts")
genes <- rownames(expression_matrix)
genes_df <- read.csv("./NSCLC_SingleCell/LUNG_jy8/genes.tsv", header = F, sep = "\t", stringsAsFactors = FALSE)
dropGenes <- genes[which(is.na(match(genes,genes_df[,2])))]
filtered_expression_matrix <- expression_matrix[!rownames(expression_matrix) %in% dropGenes, ]
genes <- rownames(filtered_expression_matrix)
which(is.na(match(genes,genes_df[,2])))
genesfile <- genes_df[(match(genes,genes_df[,2])),]
barcodes <- colnames(filtered_expression_matrix)
metadata <- experiment.aggregate@meta.data
write.table(metadata, file = "./NSCLC_SingleCell/filter_LUNG_jy8/metadata.txt", sep = "\t", row.names = TRUE)
saveRDS(metadata, file='./NSCLC_SingleCell/filter_LUNG_jy8/metadata.rds')
saveRDS(experiment.aggregate, file='./NSCLC_SingleCell/filter_LUNG_jy8/experiment_aggregate.rds')
write.table(barcodes, file = "./NSCLC_SingleCell/filter_LUNG_jy8/barcodes.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(genesfile, file = "./NSCLC_SingleCell/filter_LUNG_jy8/genes.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
writeMM(filtered_expression_matrix, file = "./NSCLC_SingleCell/filter_LUNG_jy8/matrix.mtx")

setwd("./NSCLC_SingleCell/NSCLCData/H5AD")
SaveH5Seurat(experiment.aggregate, filename = "scRNA_LUNG.h5seurat", overwrite = TRUE)
Convert("scRNA_LUNG.h5seurat", "scRNA_LUNG.h5ad", overwrite = TRUE)
