rm(list = ls())
library(Seurat)
library(ggplot2)
library(dplyr)

sample = "OTAR_LNGsp10391237"
Data = readRDS(paste0("./NSCLC_STData/CNVOutput/HNRDSData/HNsTRNA(",sample,").rds"))
Data=subset(Data,subset=cnvcluster!="NA")
EMData=subset(Data,subset=cnvcluster!="Other")

DEG <- FindMarkers(Data, ident.1 = "Normoxia",ident.2 = "Hypoxia", group.by="cnvcluster",test.use = "wilcox",slot="data",logfc.threshold=0)
DEG$cutoffvalue=mean(abs(DEG$avg_log2FC))+2*sd(abs(DEG$avg_log2FC))
set1=subset(Data,cnvcluster=="Normoxia")
set2=subset(Data,cnvcluster=="Hypoxia")
set1median=apply(set1@assays$SCT@data,1,median)
set2median=apply(set2@assays$SCT@data,1,median)
DEG$median1=as.double(set1median[rownames(DEG)])
DEG$median2=as.double(set2median[rownames(DEG)])

write.table(tibble::rownames_to_column(as.data.frame(DEG),"gene"),paste0("./NSCLC_STData/CNVOutput/DiffGenes/",sample,"DEGs.txt"),
            quote = FALSE,sep="\t",row.names=FALSE,col.names=TRUE)

SignificantData=subset(DEG,p_val_adj<0.05 & abs(DEG$avg_log2FC)>DEG$cutoff & abs(median2-median1)>0.05)
if(dim(SignificantData)[1]>0){
  for(i in 1:dim(SignificantData)[1]){
    gene=rownames(SignificantData)[i]
    p = SpatialFeaturePlot(EMData, features = gene,pt.size.factor=2) +scale_fill_gradient(high="red",low="white")
    ggsave(filename=paste0(sTDEgene_OutputFold,"/DEgeneplot/",gene,".png"),p,width = 8,height = 8,device = "png",dpi = 900,units = "in")
  }
}

rm(list = ls())
library(Seurat)
library(ggplot2)

c = c("OTAR_LNGsp10206157","OTAR_LNGsp10206158","OTAR_LNGsp10206159","OTAR_LNGsp10206160","OTAR_LNGsp10206165","OTAR_LNGsp10206166",
      "OTAR_LNGsp10206167","OTAR_LNGsp10206168","OTAR_LNGsp10391235","OTAR_LNGsp10391236","OTAR_LNGsp10391237","OTAR_LNGsp10391238")
sTcellinteraction = paste0("./NSCLC_STData/CNVOutput/DiffGenes/","OTAR_LNGsp10206157","DEGs.txt")
DiffGenes=read.table(sTcellinteraction,sep="\t",,row.names=NULL,header=TRUE,fill=FALSE,check.names=FALSE)
SignificantData=subset(DiffGenes,p_val_adj<0.05 & abs(median2-median1)>0.05)
SignData1=data.frame(cbind(SignificantData$gene,SignificantData$avg_log2FC,SignificantData$p_val_adj))
colnames(SignData1)=c("Genes","log2FC","Pjust") 
SignData1$sign=rep("OTAR_LNGsp10206157",nrow(SignData1))
for (sample in c[2:12]) {
  sTcellinteraction = paste0("./NSCLC_STData/CNVOutput/DiffGenes/",sample,"DEGs.txt")
  SignData2=read.table(sTcellinteraction,sep="\t",,row.names=NULL,header=TRUE,fill=FALSE,check.names=FALSE)
  SignData2=subset(SignData2,p_val_adj<0.05 & abs(median2-median1)>0.05)
  SignData3=data.frame(cbind(SignData2$gene,SignData2$avg_log2FC,SignData2$p_val_adj))
  colnames(SignData3)=c("Genes","log2FC","Pjust")
  SignData3$sign=rep(sample,nrow(SignData3))
  SignData1=rbind(SignData1,SignData3)
}
diffGenes=unique(SignData1$Genes)
write.table(SignData1,file = "./NSCLC_STData/statistic/diffGenes.txt",sep = "\t",col.names = T,row.names = F, quote = F)
write.table(diffGenes,file = "./NSCLC_STData/statistic/diffUniqueGenes.txt",sep = "\t",col.names = F,row.names = F, quote = F)
