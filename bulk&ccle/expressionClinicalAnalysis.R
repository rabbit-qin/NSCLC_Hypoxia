rm(list = ls())
library(survival)
library(survminer)

clinical <- read.csv(file = "./rawData/NSCLC_clinical.CSV",sep = ",",header = T,row.names = 1)
clinical <- clinical[which(as.character(clinical$treatments_pharmaceutical_treatment_or_therapy) == "no"),]
clinical <- clinical[which(as.character(clinical$treatments_radiation_treatment_or_therapy) == "no"),]
express_data <- read.csv(file = "./rawData/ExpressWithSignatureTCGA-20240724.csv",sep = ",",header = T,row.names = NULL)
express_data <- express_data[,-which(express_data[nrow(express_data),] == "Mix")]

genename <- "SELENBP1"
express_select <- express_data[which(express_data[,1] == genename),]

colnames(express_select) <- gsub("\\.","\\-",colnames(express_select))
colnames(express_select) <- sapply(colnames(express_select),function(x){substr(x,1,12)})
rownames(clinical) <- sapply(rownames(clinical),function(x){substr(x,1,12)})
ids_clin <- na.omit(match(unique(colnames(express_select)),rownames(clinical)))
clinical_data <- clinical[ids_clin,]
ids_expr <- na.omit(match(rownames(clinical_data),colnames(express_select)))
clinical_data$expression <- as.numeric(express_select[1,ids_expr])
clinical_data <- clinical_data[-which(is.na(clinical_data[,3])),]

for (i in 1:nrow(clinical_data)) {
  if(clinical_data[i,10] == "Dead"){
    clinical_data[i,3] <- clinical_data[i,12]
  }
}

clinical_data$signature <- "NA"
clinical_data$signature[clinical_data$expression >= median(as.numeric(clinical_data$expression))] <- "High express"
clinical_data$signature[clinical_data$expression < median(as.numeric(clinical_data$expression))] <- "Low express"
clinical_data$vital_status[clinical_data$vital_status == "Alive"] <- 0
clinical_data$vital_status[clinical_data$vital_status == "Dead"] <- 1
clinical_data$vital_status <- as.numeric(clinical_data$vital_status)
best_threshold_surv <- surv_cutpoint(clinical_data,
                                     time = "days_to_last_follow_up", 
                                     event = "vital_status", 
                                     variables = "expression", 
                                     minprop = 0.4,  
                                     progressbar = TRUE) 

summary(best_threshold_surv)
best_threshold_data <- surv_categorize(best_threshold_surv)
surv_obj <- Surv(time = best_threshold_data$days_to_last_follow_up, event = best_threshold_data$vital_status)
surv_fit <- survfit(surv_obj ~ best_threshold_data$expression)

png(paste0("./Figure/ExpressClinical/",genename,"_ExprSurv.png"),width=5,height=5,unit="in",res=1200)
theme_temp <- theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill="white"),plot.title = element_text(hjust = 0.5))+
  theme(legend.position = 'none') +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
ggsurvplot(surv_fit, data = best_threshold_data, 
           conf.int = TRUE, 
           pval = TRUE, 
           fun = "pct", 
           size = 1.5, 
           linetype = c("solid","dotted"), 
           palette = c("mediumblue","red"),  
           legend = c(0.25, 0.2),  
           legend.title = genename, 
           surv.median.line = "hv",
           ggtheme = theme_temp)   
dev.off()
