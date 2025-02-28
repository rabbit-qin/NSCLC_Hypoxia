rm(list = ls())

express_file <- read.csv(file = "./rawData/ExpressWithSignatureTCGA-20240724.csv",sep = ",",header = T,row.names = NULL)
CARM1_mirna_targetgene <- "CARM1"
genename <- "NM_199141|CARM1|chr19|+"
pdui_file <- read.csv(file = "./rawData/PDUIWithSignatureTCGA-20240724.csv",sep = ",",header = T,row.names=NULL)
pdui_file <- pdui_file[,-which(pdui_file[nrow(pdui_file),] == "Mix")]
express_file <- express_file[,c(1,na.omit(match(colnames(pdui_file),colnames(express_file))))]

pdui_select <- pdui_file[which(pdui_file[,1] == genename),]
express_targetgene <- express_file[na.omit(match(CARM1_mirna_targetgene,express_file[,1])),]
ids_pdui <- match(colnames(express_targetgene),colnames(pdui_select))
pdui_select <- pdui_select[,c(1,na.omit(ids_pdui))]
ids_express <- match(colnames(pdui_select),colnames(express_targetgene))
express_targetgene <- express_targetgene[,c(1,na.omit(ids_express))]
colnames(pdui_select)[1] <- "X"
fig <- data.frame(t(rbind(pdui_select,express_targetgene)))
colnames(fig) <- fig[1,]
fig <- fig[-1,]
fig <- fig[-which(is.na(fig[,1])),]
figure <- data.frame(apply(fig, 2, as.numeric))

b <- ggplot(figure, aes(x = figure[,1], y = figure[,2])) + 
  geom_point(size=1.5,shape=16) + 
  geom_smooth(method = "lm",formula=y~x,size = 2) + theme_bw() + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))  + stat_cor(method = "pearson")
ggsave(paste0("./figure/PDUIExpressCorrelation/",targetgene,"_bulkPDUI.png"),b,width = 3,height = 3.5,device = "png",dpi = 900,units = "in")


######CCLE
rm(list = ls())

express_ccle <- read.csv(file = "./rawData/ExpressWithSignatureCCLE-20240724.csv",sep = ",",header = T)
pdui_ccle <- read.csv(file = "./rawData/PDUIWithSignatureCCLE-20240724.csv",sep = ",",header = T,row.names = NULL)
CARM1_mirna_targetgene <- "CARM1"
genename <- "NM_199141|CARM1|chr19|+"
colnames(express_ccle) <- gsub(".genes.results_count","",colnames(express_ccle))
colnames(pdui_ccle) <- gsub(".wigPDUI","",colnames(pdui_ccle))
pdui_ccle <- pdui_ccle[,-which(pdui_ccle[nrow(pdui_ccle),] == "Mix")]
express_ccle <- data.frame(express_ccle[,c(na.omit(match(colnames(pdui_ccle),colnames(express_ccle))))])

pdui_select_ccle <- pdui_ccle[which(pdui_ccle[,1] == genename),]
express_targetgene_ccle <- express_ccle[na.omit(match(CARM1_mirna_targetgene,express_ccle[,1])),]
ids_pdui_ccle <- match(colnames(express_targetgene_ccle),colnames(pdui_select_ccle))
pdui_select_ccle <- pdui_select_ccle[,c(na.omit(ids_pdui_ccle))]
fig <- data.frame(t(rbind(pdui_select_ccle,express_targetgene_ccle)))
colnames(fig) <- fig[1,]
fig <- fig[-1,]
figure <- data.frame(apply(fig, 2, as.numeric))

b <- ggplot(figure, aes(x = figure[,1], y = figure[,2])) + 
  geom_point(size=1.5,shape=16) + 
  geom_smooth(method = "lm",formula=y~x,size = 2) + theme_bw() + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))  + stat_cor(method = "pearson")
ggsave(paste0("./figure/PDUIExpressCorrelation/",targetgene,"_cclePDUI.png"),b,width = 3,height = 3.5,device = "png",dpi = 900,units = "in")
