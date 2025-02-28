rm(list = ls())
library(tidyr)

PDUI <- read.csv(file = "./rawData/PDUIWithSignatureTCGA-20240724.csv",sep = ",",header = T)
PDUI_Hypoxia <- PDUI[,which(PDUI[nrow(PDUI),] == "Hypoxia")]
PDUI_Normoxic <- PDUI[,which(PDUI[nrow(PDUI),] == "Normoxic")]

APA_Events <- ""
hypoxia_Select <- data.frame(t(PDUI_Hypoxia[which(rownames(PDUI_Hypoxia) == APA_Events),]))
normoxic_Select <- data.frame(t(PDUI_Normoxic[which(rownames(PDUI_Normoxic) == APA_Events),]))
hypoxia_Select$type <- "hypoxia"
normoxic_Select$type <- "normoxic"
figure <- rbind(hypoxia_Select,normoxic_Select)
figure[,1] <- as.numeric(figure[,1])

mycon <- list(unique(figure[,2]))
p=ggboxplot(figure,x=colnames(figure)[2],y=colnames(figure)[1],color=colnames(figure)[2],
            size=2,add='nejm',palette=c('red','blue')) +
  stat_compare_means(comparisons=mycon,method="t.test",label.y=max(as.numeric(figure[,1]))+0.05) +
  theme_bw() +
  theme(legend.position = 'none') +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("./APAEventsFigure/NM_199141-CARM1.png",p,width = 3,height = 3.5,device = "png",dpi = 900,units = "in")


#CCLE
rm(list = ls())
library(tidyr)

PDUI <- read.csv(file = "./rawData/PDUIWithSignatureCCLE-20240724.csv",sep = ",",header = T)
PDUI_Hypoxia <- PDUI[,which(PDUI[nrow(PDUI),] == "Hypoxia")]
PDUI_Normoxic <- PDUI[,which(PDUI[nrow(PDUI),] == "Normoxic")]

APA_Events <- ""
hypoxia_Select <- data.frame(t(PDUI_Hypoxia[which(rownames(PDUI_Hypoxia) == APA_Events),]))
normoxic_Select <- data.frame(t(PDUI_Normoxic[which(rownames(PDUI_Normoxic) == APA_Events),]))
hypoxia_Select$type <- "hypoxia"
normoxic_Select$type <- "normoxic"
figure <- rbind(hypoxia_Select,normoxic_Select)
figure[,1] <- as.numeric(figure[,1])

mycon <- list(unique(figure[,2]))
p=ggboxplot(figure,x=colnames(figure)[2],y=colnames(figure)[1],color=colnames(figure)[2],
            size=2,add='nejm',palette=c('red','blue')) +
  stat_compare_means(comparisons=mycon,method="t.test",label.y=max(as.numeric(figure[,1]))+0.05) +
  theme_bw() +
  theme(legend.position = 'none') +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")

ggsave("./APAEventsFigure/NM_199141-CARM1.png",p,width = 3,height = 3.5,device = "png",dpi = 900,units = "in")
