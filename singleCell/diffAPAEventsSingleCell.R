rm(list = ls())

library(tidyr)
library(dplyr)

HypoxiaPDUI <- read.csv(file = "./rawData/Hypoxia_PDUI_singleCell.csv",sep = ",",header = T,row.names = 1)
NormoxicPDUI <- read.csv(file = "./rawData/Normoxic_PDUI_singleCell.csv",sep = ",",header = T,row.names = 1)

apaEvents <- data.frame(intersect(rownames(HypoxiaPDUI),rownames(NormoxicPDUI)))
HypoxiaPDUIToge <- HypoxiaPDUI[match(apaEvents[,1],rownames(HypoxiaPDUI)),]
NormoxicPDUIToge <- NormoxicPDUI[match(apaEvents[,1],rownames(NormoxicPDUI)),]

j <- 1
n <- 1
ids <- matrix(0,nrow = nrow(HypoxiaPDUIToge),ncol = 1)
p <- matrix(0,nrow = nrow(HypoxiaPDUIToge),ncol = 1)
for (i in 1:nrow(HypoxiaPDUIToge)){
  if((i != 0)){
    p_value <- t.test(as.numeric(HypoxiaPDUIToge[i,]),as.numeric(NormoxicPDUIToge[i,]))
    ids[j,1] <- i
    j <- j+1
    p[n,1] <- p_value[["p.value"]]
    n <- n+1                  
  }
}

p <- as.data.frame(p)
ids_1 <- data.frame(ids[which(ids != 0),])
p_value <- matrix(1,nrow = nrow(HypoxiaPDUIToge),ncol = 1)
for (i in 1:length(ids_1[,1])) {
  p_value[ids_1[i,1],1] <- p[i,1]
}

fold_Hypoxia <- data.frame(apply(HypoxiaPDUIToge, 1, function(x){mean(x,na.rm=T)}))
fold_Normoxic <- data.frame(apply(NormoxicPDUIToge, 1, function(x){mean(x,na.rm=T)}))
proximal_fold_change <- matrix(0,nrow = nrow(NormoxicPDUIToge),ncol = 1)
distal_fold_change <- matrix(0,nrow = nrow(HypoxiaPDUIToge),ncol = 1)
for (i in 1:length(fold_Normoxic[,1])) {
  distal_fold_change[i,1] <- log(fold_Normoxic[i,1]/fold_Hypoxia[i,1],2)
  proximal_fold_change[i,1] <- log(fold_Hypoxia[i,1]/fold_Normoxic[i,1],2)
}
log10FDR <- -log10(p_value)
figure <- data.frame(p_value,log10FDR,fold_Normoxic,fold_Hypoxia,distal_fold_change,proximal_fold_change)
figure$gene <- rownames(figure)

figure$Change="No Change"
figure$Change[figure$p_value<=0.05&figure$distal_fold_change > 0.05]="Upregulated in Normoxic"
figure$Change[figure$p_value<=0.05&figure$proximal_fold_change > 0.05]="Upregulated in Hypoxia"

figure <- separate(figure,gene, into= c("Num","Gene","Chr","Stand"),sep= "\\|")
colnames(figure) <- c("p_value","log10p_value","Nor","Hyp","distal_fold_change","proximal_fold_change","Num","Gene","Chr","Stand","Change")
write.table(figure,file = "./differencialApaEvents.csv",sep = ",",row.names = T,col.names = T)
