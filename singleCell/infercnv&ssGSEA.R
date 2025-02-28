rm(list=ls())

library(tidyr)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(GSVA)
library(dplyr)
library(readr)
library(limma)          
library(GSEABase)

setwd("./NSCLC_SingleCell/inferCNV/")
signatureInfercnv <- read.table("./cnv3Reasult/infercnv_meta.txt",sep = "\t",header = T,quote = "",check.names = F)
sc_cancer <- readRDS("./NSCLC_SingleCell/NSCLCData/rdsData/subtypesCelltypes.rds")
sc_cancer <- subset(sc_cancer, subset=ident=="Cancer")
sc_cancer@meta.data$ident <- sc_cancer@active.ident
table(sc_cancer@meta.data$ident)
gene_exp <- as.data.frame(sc_cancer@assays[["SCT"]]@counts)

emtGeneSet <- getGmt("./NSCLC_SingleCell/NSCLCData/MSigDBData/HALLMARK_HYPOXIA.v2022.1.Hs.gmt")
#emtGeneSet <- getGmt("./NSCLC_SingleCell/NSCLCData/MSigDBData/HALLMARK_GLYCOLYSIS.v2024.1.Hs.gmt")
#emtGeneSet <- getGmt("./NSCLC_SingleCell/NSCLCData/MSigDBData/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v2022.1.Hs.gmt")
expression <- apply(gene_exp, 2, as.numeric)
rownames(expression) <- rownames(gene_exp)
EsPre <- ssgseaParam(exprData=as.matrix(expression),geneSets=emtGeneSet)
Es <- gsva(EsPre)
Es <- data.frame(t(Es))
Es$type <- "NA"
Es$type <- signatureInfercnv[match(rownames(Es),signatureInfercnv$V1),3]
write.table(Es,file = paste0("./NSCLC_SingleCell/NSCLCData/csvData/","ssGSEA-HYPOXIA.csv"),quote = F,sep = ",")

mycon <- list(c('1','2'),c('2','3'),c('1','3'))
Es <- Es[order(Es$type, decreasing = F), ]
p <- ggboxplot(Es,x = 'type',y = 'HALLMARK_HYPOXIA',color = 'type',
               size = 2,add = 'nejm', palette=c('blue','red','gray')) +
     stat_compare_means(comparisons = mycon,method = "t.test",paired=F,bracket.size=0.6,step.increase=0.07) + 
     stat_compare_means(label.y = max(Es[,1])+0.3) +
     theme_bw() +
     theme(legend.position = 'none') +
     theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
     theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
ggsave("./NSCLC_SingleCell/Figure/infercnvSsgseaHYPOXIA.png",p,width = 5,height = 5,device = "png",dpi = 1200,units = "in")

metatype=signatureInfercnv[,c(1,3)]
colnames(metatype)=c("cell","infertype")
sc_cancer@meta.data$infer_type <- metatype[match(rownames(sc_cancer@meta.data),metatype$cell),"infertype"]
phe <- sc_cancer@meta.data
table(sc_cancer@meta.data$infer_type)
sc_cancer <- NormalizeData(sc_cancer) %>% FindVariableFeatures() %>% ScaleData(do.center = F)
sc_cancer <- RunPCA(sc_cancer)
set.seed(7)
sc_cancer.pca <- RunUMAP(sc_cancer, dims = 1:30) %>% FindNeighbors(dims = 1:30) %>% FindClusters()

sc_cancer.pca@meta.data$seurat_clusters <- sc_cancer@meta.data$infer_type
p_infecnv <- DimPlot(sc_cancer.pca,reduction = "umap",group.by = 'seurat_clusters',cols = c("blue", "red", "gray"),pt.size = 0.6,label = F,repel = T)+
             theme_bw() +
             theme(legend.position = 'none') +
             theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
             theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
ggsave("./NSCLC_SingleCell/Figure/infercnvUmap.png",p_infecnv,width = 5,height = 5,device = "png",dpi = 1200,units = "in")
write_rds(sc_cancer.pca,"./NSCLC_SingleCell/NSCLCData/rdsData/cancerInfercnv.rds")

hypoxia_cancer <- subset(sc_cancer.pca, subset=infer_type=="2")
normoxic_cancer <- subset(sc_cancer.pca, subset=infer_type=="1")
write_rds(hypoxia_cancer,"E:/RFiles/NSCLC_SingleCell/NSCLCData/rdsData/hypoxia_cancer.rds")
write_rds(normoxic_cancer,"E:/RFiles/NSCLC_SingleCell/NSCLCData/rdsData/normoxic_cancer.rds")
