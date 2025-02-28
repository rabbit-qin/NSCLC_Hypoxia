rm(list = ls())
library(tidyr)
library(survival)
library(survminer)

clinical <- read.csv(file = "./rawData/NSCLC_clinical.CSV",sep = ",",header = T,row.names = 1)
clinical <- clinical[which(as.character(clinical$treatments_pharmaceutical_treatment_or_therapy) == "no"),]
clinical <- clinical[which(as.character(clinical$treatments_radiation_treatment_or_therapy) == "no"),]
bulk_pdui <- read.csv(file = "./rawData/PDUIWithSignatureTCGA-20240724.csv",sep = ",",header = T,row.names = NULL)
bulk_pdui <- bulk_pdui[,-which(bulk_pdui[nrow(bulk_pdui),] == "Mix")]
bulk_pdui <- separate(bulk_pdui,row.names,into= c("Num","APA_event","Position","Location"),sep= "\\|",remove = T)
colnames(bulk_pdui) <- substr(colnames(bulk_pdui),start = 1,stop = 12)

genename <- "CARM1"
bulk_select <- bulk_pdui[which(bulk_pdui[,2] == genename),]
colnames(bulk_select) <- gsub("\\.","\\-",colnames(bulk_select))
if (sum(is.na(bulk_select[1,])) != 0) {
  bulk_select <- bulk_select[,-which(is.na(bulk_select[1,]))]
}
rownames(clinical) <- gsub("\\.","\\-",rownames(clinical))
ids_clin <- unique(na.omit(match(colnames(bulk_select),rownames(clinical))))
clinical_data <- clinical[ids_clin,]
ids_expr <- na.omit(match(rownames(clinical_data),colnames(bulk_select)))
clinical_data$pdui <- as.numeric(bulk_select[1,ids_expr])

clinical_data <- clinical_data[-which(is.na(clinical_data[,3])),]
for (i in 1:nrow(clinical_data)) {
  if(clinical_data[i,10] == "Dead"){
    clinical_data[i,3] <- clinical_data[i,12]
  }
}

clinical_data$signature <- "NA"
clinical_data$signature[clinical_data$pdui >= median(as.numeric(clinical_data$pdui))] <- "High pdui"
clinical_data$signature[clinical_data$pdui < median(as.numeric(clinical_data$pdui))] <- "Low pdui"
clinical_data$vital_status[clinical_data$vital_status == "Alive"] <- 0
clinical_data$vital_status[clinical_data$vital_status == "Dead"] <- 1

clinical_data$vital_status <- as.numeric(clinical_data$vital_status)
best_threshold_surv <- surv_cutpoint(clinical_data,
                                     time = "days_to_last_follow_up",
                                     event = "vital_status",
                                     variables = "pdui",
                                     minprop = 0.4,
                                     progressbar = TRUE)

summary(best_threshold_surv)
best_threshold_data <- surv_categorize(best_threshold_surv)
surv_obj <- Surv(time = best_threshold_data$days_to_last_follow_up, event = best_threshold_data$vital_status)
surv_fit <- survfit(surv_obj ~ best_threshold_data$pdui)

png(paste0("./",genename,"_PDUI.png"),width=5,height=5,unit="in",res=1200)
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
           palette = c("blue", "red"),
           legend = c(0.15, 0.2),
           legend.title = genename,
           legend.labs = c("High PDUI", "Low PDUI"), 
           ggtheme = theme_temp)
dev.off()
