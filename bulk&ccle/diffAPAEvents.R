rm(list = ls())
library(tidyr)
library(ggplot2)

Hypoxia_scoes_cluster <- read.csv(file = "./NCScore/NC2NMF_cluster3TCGA.csv",sep = ",",header = T)
PDUI <- read.csv(file = "./rawData/TCGA_PDUI.csv",sep = ",",header = T)
PDUI[nrow(PDUI)+1,] <- c("Type","Chr","Location",rep("NA",1013))

PDUI <- data.frame(t(PDUI))
Hypoxia_scoes_cluster$names <- rownames(Hypoxia_scoes_cluster)
PDUI$X8973[na.omit(match(Hypoxia_scoes_cluster$names[which(Hypoxia_scoes_cluster$type == "Hypoxia")],rownames(PDUI)))] <- "Hypoxia"
PDUI$X8973[na.omit(match(Hypoxia_scoes_cluster$names[which(Hypoxia_scoes_cluster$type == "Normoxic")],rownames(PDUI)))] <- "Normoxic"
PDUI$X8973[na.omit(match(Hypoxia_scoes_cluster$names[which(Hypoxia_scoes_cluster$type == "Mix")],rownames(PDUI)))] <- "Mix"

colnames(PDUI) <- PDUI[1,]
PDUI <- PDUI[-1,]
PDUI <- data.frame(t(PDUI))
type <- PDUI[nrow(PDUI),1:ncol(PDUI)]
ids_Normoxic <- which(type == "Normoxic")
ids_Hypoxia <- which(type == "Hypoxia")

ng <- rownames(PDUI)
ng <-as.data.frame(ng)
colnames(ng) <- "name"
genename <- separate(ng,name,into = c('name',"gene","chromosome","signal"),sep = "\\|")
PDUI_Normoxic <- PDUI[1:(nrow(PDUI)-1),ids_Normoxic]
PDUI_Hypoxia <- PDUI[1:(nrow(PDUI)-1),ids_Hypoxia]
exp = apply(PDUI[1:(nrow(PDUI)-1),3:ncol(PDUI)],2,as.numeric)

k <- 1
ids_1 <- matrix(0,nrow = nrow(PDUI),ncol = 1)
for (i in 1:nrow(exp)) {
  if (sum(exp[i,] == 1,na.rm = T)>400){
    ids_1[k,1] <- i
    k <- k+1
  }
}
if(sum(ids_1) == 0){
  PDUI1 <- PDUI
}else {
  ids_1 <- ids_1[which(ids_1 != 0),]
  PDUI_Normoxic <- PDUI_Normoxic[-ids_1,]
  PDUI_Hypoxia <- PDUI_Hypoxia[-ids_1,]
  PDUI1 <- PDUI[-ids_1,]
  genename <- genename[-ids_1,]
}

j <- 1
n <- 1
ids <- matrix(0,nrow = nrow(PDUI_Normoxic),ncol = 1)
p <- matrix(0,nrow = nrow(PDUI_Normoxic),ncol = 1)
for (i in 1:nrow(PDUI_Normoxic)){
  if((i != 0)){
    p_value <- t.test(as.numeric(PDUI_Normoxic[i,]),as.numeric(PDUI_Hypoxia[i,]))
    ids[j,1] <- i
    j <- j+1
    p[n,1] <- p_value[["p.value"]]
    n <- n+1                  
  }
}

p <- data.frame(p)
ids_1 <- data.frame(ids[which(ids != 0),])
p_value <- matrix(1,nrow = nrow(PDUI_Normoxic),ncol = 1)
for (i in 1:length(ids_1[,1])) {
  p_value[ids_1[i,1],1] <- p[i,1]
}

PDUI_Normoxic <- apply(PDUI_Normoxic,2,as.numeric)
PDUI_Hypoxia <- apply(PDUI_Hypoxia,2,as.numeric)
fold_Normoxic <- data.frame(apply(PDUI_Normoxic, 1, function(x){mean(x,na.rm=T)}))
fold_Hypoxia <- data.frame(apply(PDUI_Hypoxia, 1, function(x){mean(x,na.rm=T)}))
proximal_fold_change <- matrix(0,nrow = nrow(PDUI_Normoxic),ncol = 1)
distal_fold_change <- matrix(0,nrow = nrow(PDUI_Hypoxia),ncol = 1)
for (i in 1:length(fold_Normoxic[,1])) {
  distal_fold_change[i,1] <- log(fold_Hypoxia[i,1]/fold_Normoxic[i,1],2)
  proximal_fold_change[i,1] <- log(fold_Normoxic[i,1]/fold_Hypoxia[i,1],2)
}
log10FDR <- -log10(p_value)
figure <- data.frame(p_value,log10FDR,distal_fold_change,proximal_fold_change)
figure$gene <- genename[1:(nrow(genename)-1),2]

figure$Change="No Change"
figure$Change[figure$p_value<=0.05&figure$distal_fold_change>0.05]="Upregulated in Hypoxia"
figure$Change[figure$p_value<=0.05&figure$proximal_fold_change>0.05]="Upregulated in Normoxic"
table(figure$Change)

picture1 <- ggplot(figure, aes(y=log10FDR, x=distal_fold_change, color=Change))+
  geom_point()+
  scale_colour_manual(values=c("black","red","blue")) +
  geom_hline(yintercept = -log10(0.05), color='black', lty=4, lwd=1)+
  geom_vline(xintercept = 0.05,color="black", lty=4, lwd=1) +
  geom_vline(xintercept = -0.05,color="black", lty=4, lwd=1) +
  theme_classic() +
  theme(legend.position = 'none') +
  xlab("PDUI Log2 Fold-Change form Hypoxia to Normoxic") + 
  ylab("Log10P_value")+
  xlim(-3.6,3.6)+
  theme_bw() +
  theme(legend.position = 'none') +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
ggsave("./figure/differencialApaEvents.png",picture1,width = 5,height = 5,device = "png",dpi = 600,units = "in")
rownames(figure) <- genename[1:(nrow(genename)-1),1]
write.table(figure,file = "./APAEventsDiff/differencialApaEventsTCGA.csv",sep = ",",row.names = T,col.names = T)
