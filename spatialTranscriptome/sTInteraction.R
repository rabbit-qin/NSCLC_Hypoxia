rm(list = ls())
library(Seurat)
library(ggplot2)

sample = "OTAR_LNGsp9476039"
sTcellinteraction_RDSFile = paste0("./NSCLC_STData/CNVOutput/HNRDSData/HNsTRNA(",sample,").rds")
Data = readRDS(sTcellinteraction_RDSFile)
Data=subset(Data,subset=cnvcluster!="NA")
EMData=subset(Data,subset=cnvcluster!="Other")

sTcellinteraction_LRscore = paste0("./NSCLC_STData/sTinteraction/stlearnOutputCell/",sample,"/lr_interaction_score.txt")
LRscore=read.table(sTcellinteraction_LRscore,sep=",",,row.names=1,header=TRUE,fill=FALSE,check.names=FALSE)
LRscore=LRscore[colnames(EMData),]

cnvcluster=EMData@meta.data$cnvcluster
DEGindex1=which(cnvcluster=="Normoxia")
DEGindex2=which(cnvcluster=="Hypoxia")
DEGresult=t(apply(LRscore,2,function(x){Ttest=t.test(as.double(x[DEGindex1]),as.double(x[DEGindex2]),paired = FALSE);
est=as.double(Ttest$estimate);
return(c(Ttest$p.value,est,est[2]-est[1]))}))
colnames(DEGresult)=c("p","mean1","mean2","dif")
DEGresult=cbind(data.frame("LR"=rownames(DEGresult)),as.data.frame(DEGresult))
DEGresult$p_val_adj=p.adjust(as.double(DEGresult$p),method = "BH")
sTcellinteraction_OutputFold = paste0("./NSCLC_STData/sTinteraction/sTInteractionOutput/",sample)
write.table(DEGresult,paste0(sTcellinteraction_OutputFold,"/DEcellinteraction.txt"),quote = FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
SignificantData=subset(DEGresult,p_val_adj<0.05 & abs(dif)>0.05)
write.table(SignificantData,paste0(sTcellinteraction_OutputFold,"/SignificantInteraction.txt"),quote = FALSE,sep="\t",row.names=FALSE,col.names=TRUE)

EMData@assays$SCT@data=t(LRscore)
for(i in 1:dim(SignificantData)[1]){
  LR=rownames(SignificantData)[i]
  p = SpatialFeaturePlot(EMData, features = LR,pt.size.factor=1.5) +scale_fill_gradient(high="red",low="white")
  ggsave(filename=paste0(sTcellinteraction_OutputFold,"/DEcellinteractionplot/",LR,".png"),p,width = 8,height = 8,device = "png",dpi = 900,units = "in")
}


rm(list = ls())
library(Seurat)
library(ggplot2)

c = c("OTAR_LNGsp9476038","OTAR_LNGsp9476039","OTAR_LNGsp10206157","OTAR_LNGsp10206158","OTAR_LNGsp10206159","OTAR_LNGsp10206160","OTAR_LNGsp10206165",
      "OTAR_LNGsp10206166","OTAR_LNGsp10206167","OTAR_LNGsp10206168","OTAR_LNGsp10391237","OTAR_LNGsp10391238")
sTcellinteraction = paste0("./NSCLC_STData/sTinteraction/sTInteractionOutputCellChat/","OTAR_LNGsp9476038","/SignificantInteraction.txt")
LRscore=read.table(sTcellinteraction,sep="\t",,row.names=1,header=TRUE,fill=FALSE,check.names=FALSE)
colnames(LRscore)[ncol(LRscore)-1]="OTAR_LNGsp9476038"  #ncol(LRscore)-1
for (sample in c[2:12]) {
  sTcellinteraction = paste0("./NSCLC_STData/sTinteraction/sTInteractionOutputCellChat/",sample,"/SignificantInteraction.txt")
  LRscore2=read.table(sTcellinteraction,sep="\t",,row.names=1,header=TRUE,fill=FALSE,check.names=FALSE)
  LRscore=LRscore[intersect(rownames(LRscore),rownames(LRscore2)),]
  LRscore$sample=LRscore2[intersect(rownames(LRscore),rownames(LRscore2)),ncol(LRscore2)-1]  #ncol(LRscore2)-1
  colnames(LRscore)[ncol(LRscore)]=sample
}
write.table(LRscore,file = "./NSCLC_STData/statistic/stlearnTogeDiffCellChat.txt",sep = "\t",col.names = T,row.names = T, quote = F)


rm(list = ls())
library(pheatmap)
library(ggplot2)
library(ggpubr)

figure <- read.csv(file = "./NSCLC_STData/statistic/stlearnTogeDiff.txt",sep="\t",header=T,row.names = 1)
figure <- figure[,-c(1:3)]
figure <- figure[,-2]

figure <- figure[order(rownames(figure),decreasing = F),]
annotation_col <- data.frame(rownames(figure))

colors <- colorRampPalette(c("skyblue", "white", "firebrick3"))(200)
min_val <- -max(abs(na.omit(as.matrix(figure))))
max_val <- max(abs(na.omit(as.matrix(figure))))
breaks <- c(seq(min_val, -0.5, length.out = 50),
            seq(-0.5, 0.5, length.out = 101)[-c(1, 101)],
            seq(0.5, max_val, length.out = 50))
numbers_matrix <- matrix("", nrow = nrow(figure), ncol = ncol(figure))
numbers_matrix[!is.na(figure)] <- "*"

p <- pheatmap(figure,
              color = colors,
              breaks = breaks,
              cluster_row = FALSE,cluster_col = FALSE,
              display_numbers = numbers_matrix,
              #border_color = "firebrick3"
              cellwidth = 12, cellheight = 15)
p
ggsave("./NSCLC_STData/statistic/figure/stlearnTogeDiffNEW2.png",p,width=4.2,height=4.8,dpi=900)
