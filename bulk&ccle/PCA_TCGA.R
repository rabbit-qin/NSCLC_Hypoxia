rm(list = ls())
library(dplyr)

expession_data <- read.csv(file = "./rawData/expression_TCGA.csv",sep = ",",row.names = 1,header = T)
hypoxia_geneset <- read.csv(file = "./NMFScore/NMF_geneset.csv",sep = ",",header = T)
cluster <- read.csv(file = "./NCScore/NC2NMF_cluster3TCGA.csv",sep = ",",header = T)

gene <- intersect(hypoxia_geneset$NMF.marker.genes,rownames(expession_data))
hypoxiaExp <- expession_data[gene,]
hypoxiaExp <- data.frame(t(hypoxiaExp))
hypoxiaExp$type <- cluster[match(rownames(hypoxiaExp),rownames(cluster)),6]
hypoxiaExp <- hypoxiaExp[-which(hypoxiaExp$type == "Mix"),]

express = hypoxiaExp[,1:(ncol(hypoxiaExp)-1)]
express =apply(express,2,as.numeric)
log2TPM = log2(express+1)
colnames(log2TPM) <- colnames(express)
express <- log2TPM
express_pca <-prcomp(express)
df_pcs  <-data.frame(express_pca$x,type = as.factor(hypoxiaExp[,ncol(hypoxiaExp)]))

library(ggplot2)
library(tidyverse)
library(scatterplot3d)
percentVar <- express_pca$sdev^2/sum(express_pca$sdev^2)
picture1 <- ggplot(df_pcs,aes(x=PC1,y=PC2,color=type))+
  geom_point()+
  scale_colour_manual(values=c("red","blue","green")) +
  stat_ellipse(level = 0.95, show.legend = T) + 
  theme_classic() +
  xlab(paste0("PC1:", round(percentVar[1]*100), "% variance")) +
  ylab(paste0("PC2:", round(percentVar[2]*100), "% variance")) +
  theme_bw() +
  theme(legend.position = 'none') +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
ggsave("./figure/clusterPcaTcga.png",picture1,width = 5,height = 5,device = "png",dpi = 600,units = "in")
