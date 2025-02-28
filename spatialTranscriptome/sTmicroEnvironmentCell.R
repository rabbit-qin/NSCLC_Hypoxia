rm(list = ls())
library(Seurat)
library(dplyr)
library(distances)
library(tibble)
library(purrr)
library(ggplot2)
library(cowplot)

sample="OTAR_LNGsp9476038"
sTmicroenvironmentcell_RDSFile=paste0("./NSCLC_STData/CNVOutput/HNRDSData/HNsTRNA(",sample,").rds")
spData=readRDS(sTmicroenvironmentcell_RDSFile)
predictions=GetAssayData(spData,assay='predictions') %>% as.matrix() %>% t()
spatial_coord=data.frame(spData@images[[names(spData@images)]]@coordinates)
spatial_coord$imagerow_scaled=spatial_coord$imagerow*spData@images[[names(spData@images)]]@scale.factors$hires
spatial_coord$imagecol_scaled=spatial_coord$imagecol*spData@images[[names(spData@images)]]@scale.factors$hires
index=intersect(rownames(predictions),rownames(spatial_coord))
spatial_coord=spatial_coord[index,]
predictions=predictions[index,]
TumorData=subset(spData,subset=cnvcluster!="NA")
TumorData=subset(TumorData,cells=index)
TumorData=subset(TumorData,subset=cnvcluster!="Other")
EpiCell=rownames(subset(TumorData@meta.data,cnvcluster=="Normoxia"))
MesCell=rownames(subset(TumorData@meta.data,cnvcluster=="Hypoxia"))
knn=7
location=spatial_coord[,c('imagerow','imagecol')]
nnmatrix=RANN::nn2(location,k=knn)$nn.idx

micropredictions=apply(nnmatrix,1,function(x){apply(predictions[x[2:knn],],2,sum)}) %>% t()
rownames(micropredictions)=rownames(predictions)
interspot=t(apply(micropredictions,2,function(x){Ttest=t.test(as.double(x[EpiCell]),as.double(x[MesCell]),paired = FALSE);
est=as.double(Ttest$estimate);
return(c(Ttest$p.value,est,est[2]-est[1]))}))
colnames(interspot)=c("pval","mean1","mean2","Diff")
interspot=as.data.frame(interspot)
interspot$p_val_adj=p.adjust(as.double(interspot$pval),method = "BH")
sTmicroenvironmentcell_OutputFold=paste0("./NSCLC_STData/sTmicroEnvironmentCell/",sample)
write.table(tibble::rownames_to_column(as.data.frame(interspot),"cell"),paste0(sTmicroenvironmentcell_OutputFold,"/DEcells.txt"),
            quote = FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
SignificantData=subset(interspot,p_val_adj<0.05 & abs(Diff)>0.05)
DefaultAssay(TumorData)="predictions"
micropredictions=micropredictions[colnames(TumorData),]
TumorData[["predictions"]] = CreateAssayObject(t(micropredictions))
if(dim(SignificantData)[1]>0){
  for(i in 1:dim(SignificantData)[1]){
    cell=rownames(SignificantData)[i]
    p=SpatialFeaturePlot(TumorData, features = cell,pt.size.factor=1.5)+
      scale_fill_gradient(high="red",low="white")
    ggsave(filename=paste0(sTmicroenvironmentcell_OutputFold,"/DEcellplot/",cell,".png"),p,width = 8,height = 8,device = "png",dpi = 900,units = "in")
  }
}


rm(list = ls())
library(Seurat)
library(ggplot2)

c = c("OTAR_LNGsp9476038","OTAR_LNGsp9476039","OTAR_LNGsp10206157","OTAR_LNGsp10206158","OTAR_LNGsp10206159","OTAR_LNGsp10206160","OTAR_LNGsp10206165",
      "OTAR_LNGsp10206166","OTAR_LNGsp10206167","OTAR_LNGsp10206168","OTAR_LNGsp10391237","OTAR_LNGsp10391238")
sTcellinteraction = paste0("./NSCLC_STData/sTmicroEnvironmentCell/","OTAR_LNGsp9476038","/DEcells.txt")
LRscore=read.table(sTcellinteraction,sep="\t",,row.names=1,header=TRUE,fill=FALSE,check.names=FALSE)
colnames(LRscore)[ncol(LRscore)]="OTAR_LNGsp9476038"
for (sample in c[2:12]) {
  sTcellinteraction = paste0("./NSCLC_STData/sTmicroEnvironmentCell/",sample,"/DEcells.txt")
  LRscore2=read.table(sTcellinteraction,sep="\t",,row.names=1,header=TRUE,fill=FALSE,check.names=FALSE)
  LRscore=LRscore[intersect(rownames(LRscore),rownames(LRscore2)),]
  LRscore$sample=LRscore2[intersect(rownames(LRscore),rownames(LRscore2)),ncol(LRscore2)]
  colnames(LRscore)[ncol(LRscore)]=sample
}
write.table(LRscore,file = "./NSCLC_STData/statistic/sTmicroEnvironmentCellPvalue.csv",sep = ",",col.names = T,row.names = T, quote = F)


rm(list = ls())
library(pheatmap)
library(ggplot2)
library(ggpubr)

figure <- read.csv(file = "./NSCLC_STData/statistic/sTmicroEnvironmentCellDiff.csv",sep=",",header=T,row.names = 1)
figure <- figure[,-c(1:3)]
figure <- figure[,-2]
figure <- figure[order(rownames(figure),decreasing = F),]
annotation_col <- data.frame(rownames(figure))

colors <- colorRampPalette(c("skyblue", "white", "firebrick3"))(100)
min_val <- -max(abs(na.omit(as.matrix(figure))))
max_val <- max(abs(na.omit(as.matrix(figure))))
numbers_matrix <- matrix("", nrow = nrow(figure), ncol = ncol(figure))
numbers_matrix[!is.na(figure)] <- "*"

p <- pheatmap(figure,
              color = colors,
              cluster_row = FALSE,cluster_col = FALSE,
              display_numbers = numbers_matrix,
              cellwidth = 15, cellheight = 12)
p
ggsave("./NSCLC_STData/statistic/figure/sTmicroEnvironmentCellDiffNEW.png",p,width=3.8,height=3.2,dpi=900)
