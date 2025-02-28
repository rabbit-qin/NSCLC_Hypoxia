rm(list = ls())
library(survival)
library(survminer)

Hypoxia_scoes <- read.csv(file = "./rawData/expression_TCGA.csv",sep = ",",header = T,row.names = 1)
clinical_data <- read.csv(file = "./rawData/NSCLC_clinical.CSV",sep = ",",header = T,row.names = 1)
clinical_data <- clinical_data[which(as.character(clinical_data$treatments_pharmaceutical_treatment_or_therapy) == "no"),]
clinical_data <- clinical_data[which(as.character(clinical_data$treatments_radiation_treatment_or_therapy) == "no"),]

rownames(Hypoxia_scoes) <- gsub("\\.","\\-",rownames(Hypoxia_scoes))
Hypoxia_scoes$label <- substr(rownames(Hypoxia_scoes),1,12)
ids_clin <- na.omit(match(Hypoxia_scoes$label,rownames(clinical_data)))
clinical_data <- clinical_data[ids_clin,]
ids_type <- match(rownames(clinical_data),Hypoxia_scoes$label)
clinical_data <- clinical_data[-which(is.na(ids_type)),]
clinical_data$type <- Hypoxia_scoes$type[na.omit(ids_type)]
clinical_data <- clinical_data[-which(clinical_data$type == "Mix"),]
clinical_data <- clinical_data[-which(is.na(clinical_data[,3])),]

for (i in 1:nrow(clinical_data)) {
  if(clinical_data[i,10] == "Dead"){
    clinical_data[i,3] <- clinical_data[i,12]
  }
}
if(sum(which(is.na(clinical_data[,10]))) != 0){
  clinical_data <- clinical_data[-which(is.na(clinical_data[,10])),]
}

clinical_data$vital_status[clinical_data$vital_status == "Alive"] <- 0
clinical_data$vital_status[clinical_data$vital_status == "Dead"] <- 1
clinical_data$vital_status <- as.numeric(clinical_data$vital_status)

surv_model <- surv_fit(Surv(days_to_last_follow_up,vital_status) ~ type, data = clinical_data)
png(paste0("./figure/clinical.png"),width=5,height=5,unit="in",res=900)
theme_temp <- theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill="white"),plot.title = element_text(hjust = 0.5))+
  theme(legend.position = 'none') +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
ggsurvplot(surv_model, data = clinical_data, 
           conf.int = TRUE, 
           pval = TRUE,  
           fun = "pct", 
           size = 1.5, 
           linetype = "strata", 
           palette = c("red","blue"),
           legend = c(0.19, 0.2),
           surv.median.line = "hv",
           ggtheme = theme_temp)   
dev.off()
