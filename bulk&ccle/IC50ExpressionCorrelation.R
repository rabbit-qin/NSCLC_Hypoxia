rm(list = ls())
library(pRRophetic)
library(ggplot2)
library(ggpubr)

expression_data <- read.csv(file = "./rawData/ExpressWithSignatureTCGA-20240724.csv",sep = ",",header = T)
expression_data <- expression_data[,-which(expression_data[nrow(expression_data),] == "Mix")]
rownames(expression_data) <- expression_data[,1]
expression_data <- expression_data[,-1]
expression <- apply(expression_data, 2, as.numeric)
rownames(expression) <- rownames(expression_data)

predictedPtype <- pRRopheticPredict(testMatrix=expression, 
                                    drug="Cisplatin",
                                    tissueType = "all", 
                                    batchCorrect = "eb",
                                    selection=1,
                                    dataset = "cgp2014")
genename <- "CARM1"
boxplot(predictedPtype)
drug_score <- data.frame(predictedPtype)
drug_score$express <- as.numeric(expression_data[which(rownames(expression_data) == genename),match(rownames(drug_score),colnames(expression_data))])
drug_score <- drug_score[-which(is.nan(drug_score[,1])),]

max(drug_score[,1]);min(drug_score[,1])
max(drug_score[,2]);min(drug_score[,2])
b <- ggplot(drug_score, aes(x = drug_score[,2], y = drug_score[,1])) + 
  geom_point(size=1.5,shape=16) + 
  geom_smooth(method = "lm",formula=y~x,size = 2) + theme_bw() + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))  + stat_cor(method = "pearson")
ggsave(paste0("./",genename,"_drug_Cisplatin_correlation.png"),b,width = 3,height = 3.5,device = "png",dpi = 900,units = "in")
