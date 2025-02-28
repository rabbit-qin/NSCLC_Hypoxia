rm(list = ls())
library(Seurat)
library(DropletUtils)

sample = "OTAR_LNGsp9476039"
rds210X_RDSFile = paste0("./NSCLC_STData/CNVOutput/HNRDSData/HNsTRNA(",sample,").rds")
Data = readRDS(rds210X_RDSFile)
rds210X_OutputFold = paste0("./NSCLC_STData/sTinteraction/preData/Interaction(",sample,")")
write10xCounts(rds210X_OutputFold,Data@assays$SCT@data,version="3")
Data=subset(Data,subset=cnvcluster!="NA")
write.table(tibble::rownames_to_column(as.data.frame(Data@meta.data),"barcode"),paste0(rds210X_OutputFold,"/metadata.txt"),
            quote = FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
