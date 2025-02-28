rm(list = ls())

express_file <- read.csv(file = "./rawData/ExpressWithSignatureTCGA-20240724.csv",sep = ",",header = T,row.names = NULL)
express_file <- express_file[,-which(express_file[nrow(express_file),] == "Mix")]

genename <- "CARM1"
express_select <- express_file[which(express_file[,1] == genename),]
targetgene <- "SELENBP1"
express_select_target <- express_file[which(express_file[,1] == targetgene),]
figure <- data.frame(t(express_select[1,2:ncol(express_select)]),t(express_select_target[1,2:ncol(express_select_target)]))
figure <- data.frame(apply(figure, 2, as.numeric))

max(figure[,1]);min(figure[,1])
max(figure[,2]);min(figure[,2])
b <- ggplot(figure, aes(x = figure[,1], y = figure[,2])) + 
  geom_point(size=1.5,shape=16) + 
  geom_smooth(method = "lm",formula=y~x,size = 2) + theme_bw() + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))  + stat_cor(method = "pearson")
ggsave(paste0("./figure/ExpressCorrelation/",targetgene,"_bulkexpression.png"),b,width = 3,height = 3.5,device = "png",dpi = 900,units = "in")


######CCLE
rm(list = ls())

express_ccle <- read.csv(file = "./rawData/ExpressWithSignatureCCLE-20240724.csv",sep = ",",header = T) 
express_ccle <- express_ccle[,-which(express_ccle[nrow(express_ccle),] == "Mix")]

genename <- "CARM1"
express_select <- express_ccle[which(express_ccle[,1] == genename),]
targetgene <- "SELENBP1"
express_select_target <- express_ccle[which(express_ccle[,1] == targetgene),]
figure <- data.frame(t(express_select[1,2:ncol(express_select)]),t(express_select_target[1,2:ncol(express_select_target)]))
figure <- data.frame(apply(figure, 2, as.numeric))

max(figure[,1]);min(figure[,1])
max(figure[,2]);min(figure[,2])
b <- ggplot(figure, aes(x = figure[,1], y = figure[,2])) + 
  geom_point(size=1.5,shape=16) +
  geom_smooth(method = "lm",formula=y~x,size = 2) + theme_bw() + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))  + stat_cor(method = "pearson")
ggsave(paste0("./figure/ExpressCorrelation/",targetgene,"_ccleExpression.png"),b,width = 3,height = 3.5,device = "png",dpi = 900,units = "in")
