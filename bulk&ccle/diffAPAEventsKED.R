rm(list = ls())
library(tidyr)
library(ggplot2)

Hypoxia_scoes_cluster <- read.csv(file = "./NCScore/NC2NMF_cluster3TCGA.csv",sep = ",",header = T)
PDUI <- read.csv(file = "./rawData/TCGA_PDUI.csv",sep = ",",header = T)
diffAPA <- read.csv(file = "./APAEventsDiff/differencialApaEventsTCGA.csv",sep = ",",header = T)
diffAPA <- diffAPA[-which(diffAPA$Change == "No Change"),]
PDUI <- separate(PDUI,event_id,into = c('name',"gene","chromosome","signal"),sep = "\\|")
PDUI <- PDUI[na.omit(match(rownames(diffAPA),PDUI[,1])),]
PDUI[nrow(PDUI)+1,] <- c("Type","gene","chromosome","signal","Chr","Location",rep("NA",1013))

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
expr_Normoxic <- PDUI[1:(nrow(PDUI)-1),ids_Normoxic]
expr_Hypoxia <- PDUI[1:(nrow(PDUI)-1),ids_Hypoxia]

pi <- apply(expr_Normoxic,2,as.numeric)
expr_meanNormoxic <- matrix(0,nrow = nrow(pi),ncol = 1)
for (i in 1:nrow(pi)) {
  expr_meanNormoxic[i,1] <- mean(pi[i,],na.rm = T)
}
es <- apply(expr_Hypoxia,2,as.numeric)
expr_meanHypoxia <- matrix(0,nrow = nrow(es),ncol = 1)
for (i in 1:nrow(pi)) {
  expr_meanHypoxia[i,1] <- mean(es[i,],na.rm = T)
}
fig_number <- data.frame(expr_meanNormoxic,expr_meanHypoxia)
fig_number <- fig_number[-nrow(fig_number),]

library(sm)
png(filename='./figure/differencialApaDensityTCGA.png',units="in",width=5, height=5,res=1200)
par(lwd = 3)
a1 <- density(fig_number$expr_meanNormoxic)
a2 <- density(fig_number$expr_meanHypoxia)
hist(fig_number$expr_meanNormoxic,freq = F,breaks = 60,col = rgb(0,0,0.6,0.5),main = "Hypoxia-Normoxic differential apa distribution",
     sub = "Two-sample Kolmogorov-Smirnov test,p-value < 2.2e-16",
     xlim = c(0,1),ylim = c(0,2),xlab = "mean PDUI value")
hist(fig_number$expr_meanHypoxia,freq = F,breaks = 60,col = rgb(1,0,0,0.5),add = T)
abline(v=mean(fig_number$expr_meanHypoxia),col = "red",lty = 1)
abline(v=mean(fig_number$expr_meanNormoxic),col = "blue",lty = 2)
lines(a1,col = "blue")
lines(a2,col = "red")
legend(0,2.0,title = "Type",c("Hypoxia","Normoxic"),lty = c(1,2),col = c("red","blue"))
dev.off()
