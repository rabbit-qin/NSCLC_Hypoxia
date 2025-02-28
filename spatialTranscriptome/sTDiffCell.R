rm(list = ls())
library(Seurat)
library(ggplot2)
library(dplyr)

sample = "OTAR_LNGsp9476039"
Data = readRDS(paste0("./NSCLC_STData/CNVOutput/HNRDSData/HNsTRNA(",sample,").rds"))
Data=subset(Data,subset=cnvcluster!="NA")
EMData=subset(Data,subset=cnvcluster!="Other")
EMData = Data
NorCell=rownames(subset(EMData@meta.data,cnvcluster=="Normoxia"))
HypCell=rownames(subset(EMData@meta.data,cnvcluster=="Hypoxia"))
predictions=GetAssayData(EMData,assay='predictions') %>% as.matrix() %>% t()
intraspot=t(apply(predictions,2,function(x){Ttest=t.test(as.double(x[NorCell]),as.double(x[HypCell]),paired = FALSE);
est=as.double(Ttest$estimate);return(c(Ttest$p.value,est,est[2]-est[1]))}))
colnames(intraspot)=c("pval","mean1","mean2","Diff")
intraspot=as.data.frame(intraspot)
intraspot$p_val_adj=p.adjust(as.double(intraspot$pval),method = "BH")

sTDEcell_OutputFold = paste0("./NSCLC_STData/CNVOutput/DiffCell/",sample)
write.table(tibble::rownames_to_column(as.data.frame(intraspot),"cell"),paste0(sTDEcell_OutputFold,"/DEcells.txt"),
            quote = FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
SignificantData=subset(intraspot,p_val_adj<0.05 & abs(Diff)>0.05)
DefaultAssay(EMData)="predictions"
if(dim(SignificantData)[1]>0){
  for(i in 1:dim(SignificantData)[1]){
    cell=rownames(SignificantData)[i]
    p = SpatialFeaturePlot(EMData, features = cell,pt.size.factor=1.5)+
      scale_fill_gradient(high="red",low="white")
    ggsave(filename=paste0(sTDEcell_OutputFold,"/DEcellplot/",cell,".png"),p,width = 8,height = 8,device = "png",dpi = 900,units = "in")
  }
}


rm(list = ls())
library(Seurat)
library(ggplot2)

c = c("OTAR_LNGsp9476038","OTAR_LNGsp9476039","OTAR_LNGsp10206157","OTAR_LNGsp10206158","OTAR_LNGsp10206159","OTAR_LNGsp10206160","OTAR_LNGsp10206165",
      "OTAR_LNGsp10206166","OTAR_LNGsp10206167","OTAR_LNGsp10206168","OTAR_LNGsp10391235","OTAR_LNGsp10391236","OTAR_LNGsp10391237","OTAR_LNGsp10391238")
sTcellinteraction = paste0("./NSCLC_STData/CNVOutput/DiffCell/","OTAR_LNGsp9476038","/DEcells.txt")
LRscore=read.table(sTcellinteraction,sep="\t",,row.names=1,header=TRUE,fill=FALSE,check.names=FALSE)
LRscore=LRscore[,4:5]
colnames(LRscore)=c("OTAR_LNGsp9476038_Diff","OTAR_LNGsp9476038_Pjust") 
for (sample in c[2:14]) {
  sTcellinteraction = paste0("./NSCLC_STData/CNVOutput/DiffCell/",sample,"/DEcells.txt")
  LRscore2=read.table(sTcellinteraction,sep="\t",,row.names=1,header=TRUE,fill=FALSE,check.names=FALSE)
  LRscore=LRscore[intersect(rownames(LRscore),rownames(LRscore2)),]
  LRscore$sample=LRscore2[intersect(rownames(LRscore),rownames(LRscore2)),ncol(LRscore2)-1]
  colnames(LRscore)[ncol(LRscore)]=paste0(sample,"_Diff")
  LRscore$sample=LRscore2[intersect(rownames(LRscore),rownames(LRscore2)),ncol(LRscore2)]
  colnames(LRscore)[ncol(LRscore)]=paste0(sample,"_Pjust")
}
write.table(LRscore,file = "./NSCLC_STData/statistic/cellTypesToge.txt",sep = "\t",col.names = T,row.names = T, quote = F)


rm(list = ls())
library(pheatmap)
library(ggplot2)
library(ggpubr)

figure <- read.csv(file = "./NSCLC_STData/statistic/cellTypesToge.txt",sep="\t",header=T,row.names = 1)
figure <- figure[,seq(1,28,2)]
annotation_col <- data.frame(rownames(figure))
#numbers_matrix <- matrix("", nrow = nrow(figure), ncol = ncol(figure))
#numbers_matrix[!is.na(figure)] <- "*"

p <- pheatmap(figure,
              color = colorRampPalette(c("skyblue", "white", "firebrick3"))(100),
              cluster_row = FALSE,cluster_col = FALSE,
              na_col = "grey",
              #display_numbers = numbers_matrix,
              #border_color = "firebrick3"
              cellwidth = 15, cellheight = 12)
p
ggsave("./NSCLC_STData/statistic/figure/cellTypesTogeFigureNEW2.png",p,width=4.2,height=3.5,dpi=900)
