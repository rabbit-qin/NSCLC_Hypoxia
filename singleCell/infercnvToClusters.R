rm(list=ls())

options(stringsAsFactors = F)
library(phylogram)
library(gridExtra)
library(grid)
require(dendextend)
require(ggthemes)
library(tidyverse)
library(Seurat)
library(infercnv)
library(miscTools)

infercnv.dend <- read.dendrogram(file = "./infercnv/cnv3_png/infercnv.observations_dendrogram.txt")
infercnv.labels <- cutree(infercnv.dend,k = 3, order_clusters_as_data = FALSE)
table(infercnv.labels)

cb_palette <- c("#ed1299", "#09f9f5", "#246b93")
the_bars <- as.data.frame(cb_palette[1:4][infercnv.labels])
colnames(the_bars) <- "inferCNV_tree"
the_bars$inferCNV_tree <- as.character(the_bars$inferCNV_tree)

infercnv.labels <- as.data.frame(infercnv.labels)
setwd("./infercnv")
groupFiles <- './meta.csv'   
meta <- read.table(groupFiles,sep = ',')
meta$V1 <- rownames(meta)
infercnv.labels$V1 <- rownames(infercnv.labels)
meta1 <- merge(meta,infercnv.labels,by='V1')
write.table(meta1,"cnv3Reasult/infercnv_meta.txt",sep = "\t",row.names = F,col.names = T,quote = F)
