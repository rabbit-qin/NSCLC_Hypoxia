rm(list = ls())
library(tidyr)

express_file <- read.csv(file = "./rawData/ExpressWithSignatureTCGA-20240724.csv",sep = ",",header = T,row.names = 1)
express_file <- data.frame(t(express_file))
express_normoxic <- express_file[which(express_file$Type == "Normoxic"),]
express_hypoxia <- express_file[which(express_file$Type == "Hypoxia"),]

genename <- "SELENBP1"
ids <- which(colnames(express_normoxic) == genename)
figure <- data.frame(c(express_normoxic[,ids],express_hypoxia[,ids]))
figure$type <- c(rep("Normoxic",times = nrow(express_normoxic)),rep("Hypoxia",times = nrow(express_hypoxia)))
figure[,1] <- as.numeric(figure[,1])
mycon=list(paste0(unique(figure[,2])))
colnames(figure)[1] <- "Express"

figure$type <- factor(figure$type,levels = c("Hypoxia","Normoxic"))
p=ggboxplot(figure,x=colnames(figure)[2],y=colnames(figure)[1],color=colnames(figure)[2],
            size=2,add='nejm',palette=c('red','blue')) +
  stat_compare_means(comparisons=mycon,method="t.test",label.y=max(as.numeric(figure[,1]))+0.5) +
  theme_bw() +
  theme(legend.position = 'none') +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
ggsave(paste0("./figure/ExpressDiff/",genename,"_TCGA.png"),p,width = 3,height = 3.5,device = "png",dpi = 900,units = "in")


#CCLE
rm(list = ls())
library(tidyr)

express_file <- read.csv(file = "./rawData/ExpressWithSignatureCCLE-20240724.csv",sep = ",",header = T,row.names = 1)
express_file <- data.frame(t(express_file))
express_normoxic <- express_file[which(express_file$Type == "Normoxic"),]
express_hypoxia <- express_file[which(express_file$Type == "Hypoxia"),]

genename <- "SELENBP1"
ids <- which(colnames(express_normoxic) == genename)
figure <- data.frame(c(express_normoxic[,ids],express_hypoxia[,ids]))
figure$type <- c(rep("Normoxic",times = nrow(express_normoxic)),rep("Hypoxia",times = nrow(express_hypoxia)))
figure[,1] <- as.numeric(figure[,1])
mycon=list(paste0(unique(figure[,2])))
colnames(figure)[1] <- "Express"

figure$type <- factor(figure$type,levels = c("Hypoxia","Normoxic"))
p=ggboxplot(figure,x=colnames(figure)[2],y=colnames(figure)[1],color=colnames(figure)[2],
            size=2,add='nejm',palette=c('red','blue')) +
  stat_compare_means(comparisons=mycon,method="t.test",label.y=max(as.numeric(figure[,1]))+0.5) +
  theme_bw() +
  theme(legend.position = 'none') +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
p
ggsave(paste0("./figure/ExpressDiff/",genename,"_CCLE.png"),p,width = 3,height = 3.5,device = "png",dpi = 900,units = "in")