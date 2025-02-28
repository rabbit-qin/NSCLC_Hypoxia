rm(list = ls())
library(GSVA)
library(GSEABase)

expession_data <- read.csv(file = "./rawData/expression_TCGA.csv", sep = ",", row.names = NULL, header = T)
GenSet=getGmt("./NCScore/Hypoxia_scores.gmt")

index <- order(rowMeans(expession_data[,-1]),decreasing = T)
expr_ordered <- expession_data[index,]
keep <- !duplicated(expr_ordered[,1])
expession_data <- expr_ordered[keep,]
rownames(expession_data) <- expession_data[,1]
expession_data <- expession_data[,-1]

setnames=names(GenSet)
OutputData=t(apply(expession_data,1,function(x){y=x;index=which(x>median(x));y[index]=1;y[-index]=-1;return(y)}))
Score=OutputData[1:length(setnames),]
for(i in 1:length(setnames)){
  EachSet=setnames[i]
  genes=geneIds(GenSet[EachSet])[[EachSet]]
  if(length(genes)==1 && genes==""){
    Score[i,]="NA"
  }else{
    AllGenes=rownames(OutputData)
    index=match(genes,AllGenes)
    index=index[which(!is.na(index))]
    Score[i,]=as.vector(colSums(OutputData[index,]))
  }
}
rownames(Score)=setnames
scores <- t(data.frame("GeneSet"=rownames(Score),Score))
scores <- data.frame(scores[-1,])
write.table(scores,"./NCScore/Hypoxia_scoresTCGA.csv",quote = FALSE,sep=",",row.names=T,col.names=TRUE)

library(ggplot2)
library(ggpubr)
scores <- data.frame(apply(scores, 2, as.numeric))
b <- ggplot(scores, aes(x = scores$Winter, y = scores$Ragnum)) + 
     geom_point(size=1.5,shape=16) + 
     geom_smooth(method = "lm",formula=y~x,size = 2) + theme_bw() + 
     theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
     theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) + 
     stat_cor(method = "pearson")
target <- "Winter&Ragnum"
ggsave(paste0("./NCScore/figure/",target,"_correlationTCGA.png"),b,width = 4,height = 4,device = "png",dpi = 900,units = "in")
