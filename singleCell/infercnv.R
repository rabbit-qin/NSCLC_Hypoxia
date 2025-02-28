rm(list=ls())

library(monocle)
library(Seurat)
library(ggplot2)
library(cowplot)
Sys.setenv(JAGS_HOME="./JAGS/JAGS-4.3.0")
library(rjags)
library(infercnv)
library(ComplexHeatmap)
library(ggpubr)
library(tidyverse)

CancerSCRNA <- readRDS("./subtypesCelltypes.rds")
CancerSCRNA@meta.data$ident <- CancerSCRNA@active.ident
dat=GetAssayData(CancerSCRNA, slot='counts',assay='SCT')
groupinfo=data.frame(v1=colnames(dat), v2=CancerSCRNA@active.ident )
groupinfo$v2=as.character(groupinfo$v2)

dat <- GetAssayData(CancerSCRNA,assay = "SCT",slot = "counts")
dat <- as.data.frame(dat)
expFile=system.file("extdata", "oligodendroglioma_expression_downsampled.counts.matrix.gz", package = "infercnv")
geneFile=system.file("extdata", "oligodendroglioma_annotations_downsampled.txt", package = "infercnv")
groupFiles=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv")

library(AnnoProbe)
geneInfor=annoGene(rownames(dat),"SYMBOL",'human')
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]      
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]

dat=dat[match(geneInfor[,1], rownames(dat)),] 
rownames(geneInfor) <- geneInfor$SYMBOL   
geneInfor <- geneInfor[,-1]
meta <- subset(CancerSCRNA@meta.data,select = c("ident"))
write.table(meta, file = "./meta.csv", sep = ",", col.names = T, row.names = T)

identical(colnames(dat),rownames(meta))  
identical(rownames(dat),rownames(geneInfor))
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dat,
                                    annotations_file=meta,
                                    delim="\t",
                                    gene_order_file=geneInfor,
                                    ref_group_names=c('Fibroblast','Endothelial','Epithelial')
                                    )
setwd("./infercnv/")
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,
                             out_dir="cnv3_png/", 
                             cluster_by_groups=F,
                             denoise=T, 
                             HMM=F,
                             output_format = "png",
                             num_threads=7)
