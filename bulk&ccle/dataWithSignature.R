rm(list = ls())
library(tidyr)

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
write.table(PDUI,file = "./rawData/PDUIWithSignatureTCGA-20240724.csv",sep = ",",row.names = T,col.names = T)

express <- read.csv(file = "./rawData/expression_TCGA.csv",sep = ",",header = T)
express <- data.frame(t(express))
express$Type <- "NA"
express$Type[na.omit(match(Hypoxia_scoes_cluster$names[which(Hypoxia_scoes_cluster$type == "Hypoxia")],rownames(express)))] <- "Hypoxia"
express$Type[na.omit(match(Hypoxia_scoes_cluster$names[which(Hypoxia_scoes_cluster$type == "Normoxic")],rownames(express)))] <- "Normoxic"
express$Type[na.omit(match(Hypoxia_scoes_cluster$names[which(Hypoxia_scoes_cluster$type == "Mix")],rownames(express)))] <- "Mix"
express <- data.frame(t(express))
express[nrow(express),1] <- 'Type'
write.table(express,file = "./rawData/ExpressWithSignatureTCGA-20240724.csv",sep = ",",row.names = F,col.names = T)


rm(list = ls())
library(tidyr)

Hypoxia_scoes_cluster <- read.csv(file = "./NCScore/NC2NMF_cluster3CCLE.csv",sep = ",",header = T)
PDUI <- read.csv(file = "./rawData/CCLE_PDUI.csv",sep = ",",header = T)
PDUI[nrow(PDUI)+1,] <- c("Type","Chr","Location",rep("NA",1013))
PDUI <- data.frame(t(PDUI))
Hypoxia_scoes_cluster$names <- rownames(Hypoxia_scoes_cluster)
Hypoxia_scoes_cluster$names <- gsub("genes.results_count","wigPDUI",Hypoxia_scoes_cluster$names)
PDUI$X22065[na.omit(match(Hypoxia_scoes_cluster$names[which(Hypoxia_scoes_cluster$type == "Hypoxia")],rownames(PDUI)))] <- "Hypoxia"
PDUI$X22065[na.omit(match(Hypoxia_scoes_cluster$names[which(Hypoxia_scoes_cluster$type == "Normoxic")],rownames(PDUI)))] <- "Normoxic"
PDUI$X22065[na.omit(match(Hypoxia_scoes_cluster$names[which(Hypoxia_scoes_cluster$type == "Mix")],rownames(PDUI)))] <- "Mix"
colnames(PDUI) <- PDUI[1,]
PDUI <- PDUI[-1,]
PDUI <- data.frame(t(PDUI))
write.table(PDUI,file = "./rawData/PDUIWithSignatureCCLE-20240724.csv",sep = ",",row.names = T,col.names = T)

express <- read.csv(file = "./rawData/expression_CCLE.csv",sep = ",",header = T)
express <- data.frame(t(express))
express$Type <- "NA"
Hypoxia_scoes_cluster$names1 <- rownames(Hypoxia_scoes_cluster)
express$Type[na.omit(match(Hypoxia_scoes_cluster$names1[which(Hypoxia_scoes_cluster$type == "Hypoxia")],rownames(express)))] <- "Hypoxia"
express$Type[na.omit(match(Hypoxia_scoes_cluster$names1[which(Hypoxia_scoes_cluster$type == "Normoxic")],rownames(express)))] <- "Normoxic"
express$Type[na.omit(match(Hypoxia_scoes_cluster$names1[which(Hypoxia_scoes_cluster$type == "Mix")],rownames(express)))] <- "Mix"
express <- data.frame(t(express))
express[nrow(express),1] <- 'Type'
write.table(express,file = "./rawData/ExpressWithSignatureCCLE-20240724.csv",sep = ",",row.names = F,col.names = T)
