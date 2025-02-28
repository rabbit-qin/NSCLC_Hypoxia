rm(list = ls())

cluster <- read.csv(file = "./NMFScore/cluster2Tcga.csv",sep = ",",header = T)
Hypoxia_scores <- read.csv(file = "./NCScore/Hypoxia_scoresTCGA.csv",sep = ",",header = T)

Hypoxia_scores$cluster2 <- cluster[match(rownames(Hypoxia_scores),rownames(cluster)),2]
table(Hypoxia_scores$cluster2)
Hypoxia_scores$signature <- "Mix"
Hypoxia_scores$signature[which(Hypoxia_scores$Buffa > 0)] <- "Hypoxia"
Hypoxia_scores$signature[which(Hypoxia_scores$Buffa < 0)] <- "Normoxic"
table(Hypoxia_scores$signature)

Hypoxia_scores$type <- "Mix"
Hypoxia_scores$type[intersect(which(Hypoxia_scores$cluster2 == "hypoxia"),which(Hypoxia_scores$signature == "Hypoxia"))] <- "Hypoxia"
Hypoxia_scores$type[intersect(which(Hypoxia_scores$cluster2 == "normoxic"),which(Hypoxia_scores$signature == "Normoxic"))] <- "Normoxic"
table(Hypoxia_scores$type)

write.table(Hypoxia_scores,"./NCScore/NC2NMF_cluster3TCGA.csv",quote = T,sep=",",row.names=T,col.names=TRUE)
