rm(list = ls())
library(Seurat)
library(ggplot2)
library(dplyr)

sample = "OTAR_LNGsp9476039"
Data = readRDS(paste0("./NSCLC_STData/CNVOutput/HNRDSData/HNsTRNA(",sample,").rds"))
Data=subset(Data,subset=cnvcluster!="NA")
EMData=subset(Data,subset=cnvcluster!="Other")


SET = "EMT"
sTDEgeneset_GeneSetFile = "./NSCLC_SingleCell/NSCLCData/MSigDBData/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v2022.1.Hs.gmt"
GeneSets=read.table(sTDEgeneset_GeneSetFile,sep="\t",header=FALSE,fill=TRUE,check.names=FALSE)

for(i in 1:dim(GeneSets)[1]){
  eachgene=as.character(GeneSets[i,3:dim(GeneSets)[2]])
  eachgene=eachgene[which(eachgene!="")]
  genes=list("genes"=eachgene[eachgene %in% rownames(EMData)])
  ctrl=100
  ErrorInfo=NA
  nbins=24
  while(class(ErrorInfo)=="try-error" || is.na(ErrorInfo)){
    ErrorInfo=try(AddModuleScore(object=EMData,features=genes,name="geneset",assay = "SCT",slot = "data",ctrl=ctrl,nbins=nbins))
    if(class(ErrorInfo)=="try-error" && ctrl>0 && 
       attr(ErrorInfo,'condition')$message=="cannot take a sample larger than the population when 'replace = FALSE'"){
      ctrl=ctrl-1
    }else if(class(ErrorInfo)=="try-error"&& nbins>0 && attr(ErrorInfo,'condition')$message=="Insufficient data values to produce 24 bins."){
      nbins=nbins-1
    }else{
      break
    }
  }
  if(class(ErrorInfo)!="try-error"){
    eachgenesetscore=AddModuleScore(object=EMData,features=genes,name="geneset",assay = "SCT",slot = "data",ctrl=ctrl,nbins=nbins)
    if(i==1){
      GeneSetscore=eachgenesetscore@meta.data$geneset1
    }else{
      GeneSetscore=rbind(GeneSetscore,eachgenesetscore@meta.data$geneset1)
    }
  }
}
sTDEgeneset_OutputFold = paste0("./NSCLC_STData/CNVOutput/GeneSetScores/",sample)
if(exists('GeneSetscore')){
  GeneSetscore = data.frame(t(GeneSetscore))
  colnames(GeneSetscore)=rownames(EMData@meta.data)
  rownames(GeneSetscore)=as.character(GeneSets[,1])
  write.table(tibble::rownames_to_column(as.data.frame(GeneSetscore),"geneset"),paste0(sTDEgeneset_OutputFold,"/genesetscore",SET,".txt"),
              quote = FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
  
  cnvcluster=EMData@meta.data$cnvcluster
  DEGindex1=which(cnvcluster=="Normoxia")
  DEGindex2=which(cnvcluster=="Hypoxia")
  DEGresult=t(apply(GeneSetscore,1,function(x){Ttest=t.test(as.double(x[DEGindex1]),as.double(x[DEGindex2]),paired = FALSE);
  est=as.double(Ttest$estimate);return(c(Ttest$p.value,est,est[2]-est[1]))}))
  colnames(DEGresult)=c("p","mean1","mean2","dif")
  DEGresult=cbind(data.frame("geneset"=rownames(DEGresult)),as.data.frame(DEGresult))
  DEGresult$p_val_adj=p.adjust(as.double(DEGresult$p),method = "BH")
  SignificantData=subset(DEGresult,p_val_adj<0.05)
  EMData@assays$SCT@data=as.matrix(GeneSetscore)
  if(dim(SignificantData)[1]>0){
    for(i in 1:dim(SignificantData)[1]){
      geneset=rownames(SignificantData)[i]
      p = SpatialFeaturePlot(EMData, features = geneset,pt.size.factor=1.5) +scale_fill_gradient(high="red",low="white")
      ggsave(filename=paste0(sTDEgeneset_OutputFold,"/DEpathwayplot/",geneset,".png"),p,width = 8,height = 8,device = "png",dpi = 900,units = "in")
    }
  }
}else{
  DEGresult=data.frame(geneset="",p=0,mean1=0,mean2=0,dif=0,p_val_adj=0)
  DEGresult=DEGresult[-1,]
}
write.table(DEGresult,paste0(sTDEgeneset_OutputFold,"/DEpathway",SET,".txt"),quote = FALSE,sep="\t",row.names=FALSE,col.names=TRUE)


rm(list = ls())
library(Seurat)
library(ggplot2)

c = c("OTAR_LNGsp9476038","OTAR_LNGsp9476039","OTAR_LNGsp10206157","OTAR_LNGsp10206158","OTAR_LNGsp10206159","OTAR_LNGsp10206160","OTAR_LNGsp10206165",
      "OTAR_LNGsp10206166","OTAR_LNGsp10206167","OTAR_LNGsp10206168","OTAR_LNGsp10391235","OTAR_LNGsp10391236","OTAR_LNGsp10391237","OTAR_LNGsp10391238")
sTcellinteraction = paste0("./NSCLC_STData/CNVOutput/GeneSetScores/","OTAR_LNGsp9476038","/DEpathway.txt")
LRscore=read.table(sTcellinteraction,sep="\t",,row.names=1,header=TRUE,fill=FALSE,check.names=FALSE)
LRscore=LRscore[,4:5]
colnames(LRscore)=c("OTAR_LNGsp9476038_Diff","OTAR_LNGsp9476038_Pjust") 
for (sample in c[2:14]) {
  sTcellinteraction = paste0("./NSCLC_STData/CNVOutput/GeneSetScores/",sample,"/DEpathway.txt")
  LRscore2=read.table(sTcellinteraction,sep="\t",,row.names=1,header=TRUE,fill=FALSE,check.names=FALSE)
  LRscore=LRscore[intersect(rownames(LRscore),rownames(LRscore2)),]
  LRscore$sample=LRscore2[intersect(rownames(LRscore),rownames(LRscore2)),ncol(LRscore2)-1]
  colnames(LRscore)[ncol(LRscore)]=paste0(sample,"_Diff")
  LRscore$sample=LRscore2[intersect(rownames(LRscore),rownames(LRscore2)),ncol(LRscore2)]
  colnames(LRscore)[ncol(LRscore)]=paste0(sample,"_Pjust")
}
write.table(LRscore,file = "./NSCLC_STData/statistic/HALLMAKERToge.txt",sep = "\t",col.names = T,row.names = T, quote = F)

rm(list = ls())
library(pheatmap)
library(ggplot2)
library(ggpubr)

figure <- read.table(file = "./NSCLC_STData/statistic/HALLMAKERToge.txt",sep="\t")
figure <- figure[,seq(1,28,2)]
rownames(figure) <- c("2-EMT","3-GLY","1-HYP")
figure[2,4:6] <- NA
figure$OTAR_LNGsp10391237_Diff[2] <- NA
figure <- figure[,-(11:12)]
figure <- figure[order(rownames(figure),decreasing = F),]
annotation_col <- data.frame(rownames(figure))

min_val <- -min((na.omit(as.matrix(figure))))
max_val <- max(abs(na.omit(as.matrix(figure))))
numbers_matrix <- matrix("", nrow = nrow(figure), ncol = ncol(figure))
numbers_matrix[!is.na(figure)] <- "*"

p <- pheatmap(figure,
              color = colors,
              cluster_row = FALSE,
              cluster_col = FALSE,
              na_col = "grey",
              display_numbers = numbers_matrix,
              cellwidth = 15, cellheight = 12)
p
ggsave("./NSCLC_STData/statistic/figure/HALLMAKERTogeNEW2.png",p,width=4,height=3,dpi=900)
