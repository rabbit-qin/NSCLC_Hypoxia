rm(list = ls())
library(tidyr)
library(survival)
library(survminer)
library(GSVA)
library(GSEABase)

bulk_pdui <- read.csv(file = "./rawData/PDUIWithSignatureTCGA-20240724.csv",sep = ",",header = T,row.names = NULL)
bulk_pdui <- bulk_pdui[,-which(bulk_pdui[nrow(bulk_pdui),] == "Mix")]
bulk_pdui <- separate(bulk_pdui,row.names,into= c("Num","APA_event","Position","Location"),sep= "\\|",remove = T)
genename <- c("NM_002081","NM_005245","NM_017955","NM_006904","NM_002748","NM_003472","NM_004237","NM_002268","NM_018098","NM_021972",
              "NM_130834","NM_199141")
bulk_select <- bulk_pdui[match(genename,bulk_pdui[,1]),]

for (i in 1:(nrow(bulk_select))) {
  ids_i <- na.omit(which(is.na(bulk_select[i,])))
  bulk_select[i,ids_i] <- median(as.numeric(bulk_select[i,7:ncol(bulk_select)]),na.rm = T)
}
bulk_select[13,1:ncol(bulk_select)] <- c("Num","Gene","Position","Location",rep("NA",1015))
for (i in 7:(ncol(bulk_select))) {
  value <- 
    (4.040)*as.numeric(bulk_select[1,i])-4.343*as.numeric(bulk_select[2,i])+3.359*as.numeric(bulk_select[3,i])-5.223*as.numeric(bulk_select[4,i])+2.868*as.numeric(bulk_select[5,i])-
    4.987*as.numeric(bulk_select[6,i])+4.738*as.numeric(bulk_select[7,i])-33.227*as.numeric(bulk_select[8,i])+5.845*as.numeric(bulk_select[9,i])-5.367*as.numeric(bulk_select[10,i])-
    5.581*as.numeric(bulk_select[11,i])-3.653*as.numeric(bulk_select[12,i])
  bulk_select[13,i] <- value
}

U = mean(as.numeric(bulk_select[13,7:ncol(bulk_select)]))
sd = sd(as.numeric(bulk_select[13,7:ncol(bulk_select)]))
bulk_select[13,7:ncol(bulk_select)] = (as.numeric(bulk_select[13,7:ncol(bulk_select)])-U)/sd

library(corrplot)
correlation <- as.data.frame(t(bulk_select[,7:ncol(bulk_select)]))
colnames(correlation) <- bulk_select[,2]
cor(apply(correlation, 2, as.numeric))%>% corrplot.mixed()

bulk_select <- bulk_select[ncol(correlation),]

if(F){
  ids <- na.omit(match(colnames(bulk_select),colnames(bulk_pdui)))
  bulk_select[2,] <- bulk_pdui[nrow(bulk_pdui),ids]
  bulk_select <- data.frame(t(bulk_select[,-c(1:6)]))
  colnames(bulk_select) <- c("score","Type")
  
  bulk_select$Type <- factor(bulk_select$Type,levels = c("Hypoxia","Normoxic"))
  mycon <- list(c("Hypoxia","Normoxic"))
  bulk_select[,1] <- as.numeric(bulk_select[,1])
  bulk_select <- bulk_select[order(bulk_select$Type),]
  p=ggboxplot(bulk_select,x=colnames(bulk_select)[2],y=colnames(bulk_select)[1],color=colnames(bulk_select)[2],
              size=1.6,add='nejm',palette=c('red','blue')) + 
    stat_compare_means(comparisons=mycon,method="t.test") +
    ylim(c(-3,3)) +
    theme_bw() +
    theme(legend.position = 'none') +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
  p
  ggsave(filename="./NSCLC_Result02/Figure/bulk_boxplot.png",p,width = 4.5,height = 4.5,device = "png",dpi = 900,units = "in")
}

pdui_data <- read.csv(file = "./rawData/PDUIWithSignatureCCLE-20240724.csv",sep = ",",header = T,row.names = NULL)
pdui_data <- pdui_data[,-which(pdui_data[nrow(pdui_data),] == "Mix")]
pdui_data <- separate(pdui_data,row.names,into= c("Num","APA_event","Position","Location"),sep= "\\|",remove = T)
genename <- c("NM_002081","NM_005245","NM_017955","NM_006904","NM_002748","NM_003472","NM_004237","NM_002268","NM_018098","NM_021972",
              "NM_130834","NM_199141")
ccle_select <- pdui_data[match(genename,pdui_data[,1]),]
for (i in 1:(nrow(ccle_select))) {
  ids_i <- na.omit(which(is.na(ccle_select[i,])))
  ccle_select[i,ids_i] <- median(as.numeric(ccle_select[i,5:ncol(ccle_select)]),na.rm = T)
}

ccle_select[13,1:ncol(ccle_select)] <- c("Num","Gene","Position","Location",rep("NA",49))
for (i in 5:(ncol(ccle_select))) {
  value <- (4.040)*as.numeric(ccle_select[1,i])-4.343*as.numeric(ccle_select[2,i])+3.359*as.numeric(ccle_select[3,i])-5.223*as.numeric(ccle_select[4,i])+2.868*as.numeric(ccle_select[5,i])-
    4.987*as.numeric(ccle_select[6,i])+4.738*as.numeric(ccle_select[7,i])-33.227*as.numeric(ccle_select[8,i])+5.845*as.numeric(ccle_select[9,i])-5.367*as.numeric(ccle_select[10,i])-
    5.581*as.numeric(ccle_select[11,i])-3.653*as.numeric(ccle_select[12,i])
  ccle_select[13,i] <- value
}
rownames(ccle_select) <- ccle_select[,1]
ccle_select <- ccle_select[,-c(1,2,3,4)]
U = mean(as.numeric(ccle_select[13,]))
sd = sd(as.numeric(ccle_select[13,]))
ccle_select[13,] = (as.numeric(ccle_select[13,])-U)/sd

ccle_select <- ccle_select[ncol(correlation),]
ids <- na.omit(match(colnames(ccle_select),colnames(pdui_data)))
ccle_select[2,] <- pdui_data[nrow(pdui_data),ids]
ccle_select <- data.frame(t(ccle_select))
colnames(ccle_select) <- c("score","Type")

ccle_select$Type <- factor(ccle_select$Type,levels = c("Hypoxia","Normoxic"))
mycon <- list(c("Hypoxia","Normoxic"))
ccle_select[,1] <- as.numeric(ccle_select[,1])
ccle_select <- ccle_select[order(ccle_select$Type),]

p=ggboxplot(ccle_select,x=colnames(ccle_select)[2],y=colnames(ccle_select)[1],color=colnames(ccle_select)[2],
            size=1.8,add='nejm',palette=c('red','blue')) +
  stat_compare_means(comparisons=mycon,method="t.test",label.y=max(as.numeric(ccle_select[,1]))+0.01) +
  ylim(c(-2.5,2.5)) +
  theme_bw() +
  theme(legend.position = 'none') +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
p
ggsave(filename="./NSCLC_Result02/Figure/ccle_boxplot.png",p,width = 4.5,height = 4.5,device = "png",dpi = 900,units = "in")
