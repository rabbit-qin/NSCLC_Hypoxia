rm(list=ls())

library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(tidyr)

HypoxiaPDUI <- read.csv(file = "./rawData/Hypoxia_PDUI_singleCell.csv",sep = ",",header = T,row.names = 1)
NormoxicPDUI <- read.csv(file = "./rawData/Normoxic_PDUI_singleCell.csv",sep = ",",header = T,row.names = 1)
apaEvents <- data.frame(intersect(rownames(HypoxiaPDUI),rownames(NormoxicPDUI)))
HypoxiaPDUIToge <- HypoxiaPDUI[match(apaEvents[,1],rownames(HypoxiaPDUI)),]
NormoxicPDUIToge <- NormoxicPDUI[match(apaEvents[,1],rownames(NormoxicPDUI)),]
differencialApa = read.csv("./differencialApaEvents.csv",sep = ",",row.names = 1,header = T)
differencialApa = differencialApa[which(differencialApa$p_value <= 0.05),]
Hypoxia_scPDUI <- HypoxiaPDUIToge[match(rownames(differencialApa),rownames(HypoxiaPDUIToge)),]
Normoxic_scPDUI <- NormoxicPDUIToge[match(rownames(differencialApa),rownames(NormoxicPDUIToge)),]

pi <- apply(Normoxic_scPDUI,2,as.numeric)
expr_meanNormoxic <- matrix(0,nrow = nrow(pi),ncol = 1)
for (i in 1:nrow(pi)) {
  expr_meanNormoxic[i,1] <- mean(pi[i,],na.rm = T)
}
es <- apply(Hypoxia_scPDUI,2,as.numeric)
expr_meanHypoxia <- matrix(0,nrow = nrow(es),ncol = 1)
for (i in 1:nrow(pi)) {
  expr_meanHypoxia[i,1] <- mean(es[i,],na.rm = T)
}

ids_1 <- which(expr_meanNormoxic>=0.5)
ids_2 <- which(expr_meanHypoxia>=0.5)
expr_meanNormoxic_1 <- expr_meanNormoxic[ids_1,]
expr_meanHypoxia_1 <- expr_meanHypoxia[ids_2,]
  
library(sm)
png(filename='./figure/density.png',units="in",width=5, height=5,res=1200)
par(lwd = 3)
a1 <- density(expr_meanNormoxic_1)
a2 <- density(expr_meanHypoxia_1)
picture2 <- hist(expr_meanHypoxia_1,freq = F,breaks = 60,col = rgb(1,0,0,0.5),main = "Hyp-Nor differential apa distribution",
                 sub = "Two-sample Kolmogorov-Smirnov test,p-value = 0.05596",
                 xlim = c(0.5,1),xlab = "mean PDUI value",ylim = c(0,11))
hist(expr_meanNormoxic_1,freq = F,breaks = 60,col = rgb(0,0,1,0.5),add = T)
abline(v=mean(expr_meanNormoxic_1),col = "blue",lty = 2)
abline(v=mean(expr_meanHypoxia_1),col = "red",lty = 2)
lines(a1,col = "blue")
lines(a2,col = "red")
legend(0.8,10,title = "Type",c("Normoxic","Hypoxia"),lty = c(1,1),col = c("blue","red"))
dev.off()
