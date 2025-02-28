rm(list = ls())
library(tidyr)
library(survival)
library(survminer)
library(GSVA)
library(GSEABase)

bulk_pdui <- read.csv(file = "./rawData/PDUIWithSignatureTCGA-20240724.csv",sep = ",",header = T,row.names = NULL)
bulk_pdui <- bulk_pdui[,-which(bulk_pdui[nrow(bulk_pdui),] == "Mix")]
bulk_pdui <- separate(bulk_pdui,row.names,into= c("Num","APA_event","Position","Location"),sep= "\\|",remove = T)
genename <- c("NM_002081","NM_005245","NM_017955","NM_006904","NM_002748","NM_003472","NM_004237","NM_002268","NM_018098","NM_021972",
              "NM_130834","NM_199141")
bulk_select <- bulk_pdui[match(genename,bulk_pdui[,1]),]

for (i in 1:(nrow(bulk_select))) {
  ids_i <- na.omit(which(is.na(bulk_select[i,])))
  bulk_select[i,ids_i] <- median(as.numeric(bulk_select[i,7:ncol(bulk_select)]),na.rm = T)
}
bulk_select[13,1:ncol(bulk_select)] <- c("Num","Gene","Position","Location",rep("NA",1015))
for (i in 7:(ncol(bulk_select))) {
  value <- 
    (4.040)*as.numeric(bulk_select[1,i])-4.343*as.numeric(bulk_select[2,i])+3.359*as.numeric(bulk_select[3,i])-5.223*as.numeric(bulk_select[4,i])+2.868*as.numeric(bulk_select[5,i])-
    4.987*as.numeric(bulk_select[6,i])+4.738*as.numeric(bulk_select[7,i])-33.227*as.numeric(bulk_select[8,i])+5.845*as.numeric(bulk_select[9,i])-5.367*as.numeric(bulk_select[10,i])-
    5.581*as.numeric(bulk_select[11,i])-3.653*as.numeric(bulk_select[12,i])
  bulk_select[13,i] <- value
}

U = mean(as.numeric(bulk_select[13,7:ncol(bulk_select)]))
sd = sd(as.numeric(bulk_select[13,7:ncol(bulk_select)]))
bulk_select[13,7:ncol(bulk_select)] = (as.numeric(bulk_select[13,7:ncol(bulk_select)])-U)/sd
bulk_select <- as.data.frame(t(bulk_select[nrow(bulk_select),7:ncol(bulk_select)]))

expession_data <- read.csv(file = "./rawData/ExpressWithSignatureTCGA-20240724.csv",sep = ",",row.names = 1,header = T)
expession_data <- expession_data[-nrow(expession_data),]
gene_exp <- apply(expession_data, 2, as.numeric)
rownames(gene_exp) <- rownames(expession_data)
hypoxiaGeneSet <- getGmt("./MSigDBData/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v2022.1.Hs.gmt")
EsPre <- ssgseaParam(exprData = as.matrix(gene_exp), geneSets = hypoxiaGeneSet)
Es <- gsva(EsPre)
Es <- data.frame(t(Es))

ids <- na.omit(match(rownames(bulk_select),rownames(Es)))
Es <- data.frame(Es[ids,])

figure <- cbind(bulk_select,Es)
figure <- data.frame(apply(figure, 2, as.numeric))
figure[,2] <- (figure[,2]-mean(figure[,2]))/sd(figure[,2])
b <- ggplot(figure, aes(x = figure[,1], y = figure[,2])) + 
  geom_point(size=1.5,shape=16) + 
  #xlim(3,18.5) + 
  ylim(-4,3.1) +
  geom_smooth(method = "lm",formula=y~x,size = 2) + theme_bw() + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) + stat_cor(method = "pearson")
b
ggsave(filename="./NSCLC_Result02/Figure/TCGA_GSVA&ES_EMT.png",b,width = 4.5,height = 4.5,device = "png",dpi = 1200,units = "in")

######CCLE
pdui_data <- read.csv(file = "./rawData/PDUIWithSignatureCCLE-20240724.csv",sep = ",",header = T,row.names = NULL)
pdui_data <- pdui_data[,-which(pdui_data[nrow(pdui_data),] == "Mix")]
pdui_data <- separate(pdui_data,row.names,into= c("Num","APA_event","Position","Location"),sep= "\\|",remove = T)
genename <- c("NM_002081","NM_005245","NM_017955","NM_006904","NM_002748","NM_003472","NM_004237","NM_002268","NM_018098","NM_021972",
              "NM_130834","NM_199141")
ccle_select <- pdui_data[match(genename,pdui_data[,1]),]
for (i in 1:(nrow(ccle_select))) {
  ids_i <- na.omit(which(is.na(ccle_select[i,])))
  ccle_select[i,ids_i] <- median(as.numeric(ccle_select[i,5:ncol(ccle_select)]),na.rm = T)
}

ccle_select[13,1:ncol(ccle_select)] <- c("Num","Gene","Position","Location",rep("NA",49))
for (i in 5:(ncol(ccle_select))) {
  value <- (4.040)*as.numeric(ccle_select[1,i])-4.343*as.numeric(ccle_select[2,i])+3.359*as.numeric(ccle_select[3,i])-5.223*as.numeric(ccle_select[4,i])+2.868*as.numeric(ccle_select[5,i])-
    4.987*as.numeric(ccle_select[6,i])+4.738*as.numeric(ccle_select[7,i])-33.227*as.numeric(ccle_select[8,i])+5.845*as.numeric(ccle_select[9,i])-5.367*as.numeric(ccle_select[10,i])-
    5.581*as.numeric(ccle_select[11,i])-3.653*as.numeric(ccle_select[12,i])
  ccle_select[13,i] <- value
}
rownames(ccle_select) <- ccle_select[,1]
ccle_select <- ccle_select[,-c(1,2,3,4)]
U = mean(as.numeric(ccle_select[13,]))
sd = sd(as.numeric(ccle_select[13,]))
ccle_select[13,] = (as.numeric(ccle_select[13,])-U)/sd
colnames(ccle_select) <- gsub("\\.wigPDUI","",colnames(ccle_select))
ccle_select <- as.data.frame(t(ccle_select[nrow(ccle_select),]))

express_file <- read.csv(file = "./rawData/ExpressWithSignatureCCLE-20240724.csv",sep = ",",header = T,row.names = 1)
express_file <- express_file[-nrow(express_file),]
gene_exp <- apply(express_file, 2, as.numeric)
rownames(gene_exp) <- rownames(express_file)
hypoxiaGeneSet <- getGmt("./MSigDBData/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v2022.1.Hs.gmt")
EsPre <- ssgseaParam(exprData = as.matrix(gene_exp), geneSets = hypoxiaGeneSet)
Es <- gsva(EsPre)
Es <- data.frame(t(Es))
rownames(Es) <- gsub("\\.genes.results_count","",rownames(Es))

ids <- na.omit(match(rownames(ccle_select),rownames(Es)))
Es <- data.frame(Es[ids,])

figure <- cbind(ccle_select,Es)
figure <- data.frame(apply(figure, 2, as.numeric))
figure[,2] <- (figure[,2]-mean(figure[,2]))/sd(figure[,2])
b <- ggplot(figure, aes(x = figure[,1], y = figure[,2])) + 
  geom_point(size=1.5,shape=16) + 
  #xlim(-3,3) + 
  ylim(-2.3,2.1) +
  geom_smooth(method = "lm",formula=y~x,size = 2) + theme_bw() + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))  + stat_cor(method = "pearson")
b
ggsave(filename="./NSCLC_Result02/Figure/CCLE_GSVA&ES_EMT.png",b,width = 4.5,height = 4.5,device = "png",dpi = 1200,units = "in")
