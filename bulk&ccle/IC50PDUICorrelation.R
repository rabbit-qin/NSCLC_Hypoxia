rm(list = ls())
library(pRRophetic)
library(ggplot2)
library(ggpubr)
library(tidyr)

expression_data <- read.csv(file = "./rawData/ExpressWithSignatureTCGA-20240724.csv",sep = ",",header = T)
rownames(expression_data) <- expression_data[,1]
expression_data <- expression_data[,-1]
pdui_data <- read.csv(file = "./rawData/PDUIWithSignatureTCGA-20240724.csv",sep = ",",header = T)
pdui_data <- pdui_data[,-which(pdui_data[nrow(pdui_data),] == "Mix")]
expression_data <- expression_data[,na.omit(match(colnames(pdui_data),colnames(expression_data)))]
expression <- apply(expression_data, 2, as.numeric)
rownames(expression) <- rownames(expression_data)

predictedPtype <- pRRopheticPredict(testMatrix=expression, 
                                    drug="Cisplatin",
                                    tissueType = "all", 
                                    batchCorrect = "eb",
                                    selection=1,
                                    dataset = "cgp2014")
boxplot(predictedPtype)
drug_score <- data.frame(predictedPtype)
pdui_data[,1] <- rownames(pdui_data)
pdui <- separate(pdui_data,loci.x,into= c("Num","APA_event","Position","Location"),sep= "\\|",remove = T)
gene_apa <- "NM_199141"
pdui_select <- pdui[which(pdui[,1] == gene_apa),]

drug_score$pdui <- as.numeric(pdui_select[,match(rownames(drug_score),colnames(pdui_select))])
drug_score <- drug_score[-which(is.nan(drug_score[,1])),]
max(drug_score[,1]);min(drug_score[,1])
max(drug_score[,2]);min(drug_score[,2])

b <- ggplot(drug_score, aes(x = drug_score[,2], y = drug_score[,1])) + 
     geom_point(size=1.5,shape=16) + 
     geom_smooth(method = "lm",formula=y~x,size = 2) + theme_bw() + 
     theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
     theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))  + stat_cor(method = "pearson")
ggsave(paste0("./",gene_apa,"_drug_Cisplatin_correlation.png"),b,width = 3,height = 3.5,device = "png",dpi = 900,units = "in")
