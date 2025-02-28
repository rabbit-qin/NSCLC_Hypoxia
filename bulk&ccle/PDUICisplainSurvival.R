rm(list = ls())
library(survival)
library(survminer)
library(tidyr)

NSCLC.drug <- read.csv(file = "./rawData/NSCLC_clinical_drug.csv",header = T,sep = ",")
NSCLC.drug$drug_name <- tolower(NSCLC.drug$drug_name)
drug.freq <- data.frame(table(NSCLC.drug$drug_name))
drug.freq <- drug.freq[-1,]
k <- 1
ids <- matrix(0,1,nrow = nrow(drug.freq))
for (i in 1:nrow(drug.freq)) {
  NSCLC.drug.select <- NSCLC.drug[which(NSCLC.drug$drug_name == drug.freq[i,1]),]
  if (length(unique(NSCLC.drug.select$bcr_patient_barcode)) >= 10) {
    ids[k,1] <- i
    k <- k+1
  }
}
drug.freq <- drug.freq[ids[which(ids[,1] != 0),1],]

clinical_data <- read.csv(file = "./rawData/NSCLC_clinical.CSV",sep = ",",header = T,row.names = 1)
clinical_data$drug.note <- 0
clinical_data$drug.note[match(unique(NSCLC.drug$bcr_patient_barcode),rownames(clinical_data))] <- 1
for (i in 1:nrow(drug.freq)) {
  NSCLC.drug.select <- NSCLC.drug[which(NSCLC.drug$drug_name == drug.freq[i,1]),]
  clinical_data[,ncol(clinical_data)+1] <- matrix(1,1,nrow = nrow(clinical_data))
  clinical_data[na.omit(match(unique(NSCLC.drug.select$bcr_patient_barcode),rownames(clinical_data))),ncol(clinical_data)] <- 0
}
colnames(clinical_data)[(ncol(clinical_data)-nrow(drug.freq)+1):ncol(clinical_data)] <- as.character(drug.freq$Var1)
clinical_data <- clinical_data[which(clinical_data$drug.note == 1),]

bulk_pdui <- read.csv(file = "./rawData/PDUIWithSignatureTCGA-20240724.csv",sep = ",",header = T,row.names = NULL)
bulk_pdui <- separate(bulk_pdui,row.names,into= c("Num","APA_event","Position","Location"),sep= "\\|",remove = T)
colnames(bulk_pdui) <- substr(colnames(bulk_pdui),start = 1,stop = 12)
genename <- "CARM1"
bulk_select <- bulk_pdui[which(bulk_pdui[,2] == genename),]

colnames(bulk_select) <- gsub("\\.","\\-",colnames(bulk_select))
if (sum(is.na(bulk_select[1,])) != 0) {
  bulk_select <- bulk_select[,-which(is.na(bulk_select[1,]))]
}
rownames(clinical_data) <- gsub("\\.","\\-",rownames(clinical_data))
ids_clin <- unique(na.omit(match(colnames(bulk_select),rownames(clinical_data))))
clinical_data <- clinical_data[ids_clin,]
ids_expr <- na.omit(match(rownames(clinical_data),colnames(bulk_select)))
clinical_data$pdui <- as.numeric(bulk_select[1,ids_expr])
clinical_data <- clinical_data[-which(is.na(clinical_data[,3])),]

for (i in 1:nrow(clinical_data)) {
  if(clinical_data[i,10] == "Dead"){
    clinical_data[i,3] <- clinical_data[i,12]
  }
}

clinical_data <- clinical_data[which(clinical_data$cisplatin == 0),]
clinical_data$signature <- "NA"
clinical_data$signature[clinical_data$pdui >= mean(as.numeric(clinical_data$pdui))] <- "High PDUI"
clinical_data$signature[clinical_data$pdui < mean(as.numeric(clinical_data$pdui))] <- "Low PDUI"
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
best_threshold_data$pdui <- factor(best_threshold_data$pdui, levels = c("low", "high"))
surv_obj <- Surv(time = best_threshold_data$days_to_last_follow_up, event = best_threshold_data$vital_status)
surv_fit <- survfit(surv_obj ~ best_threshold_data$pdui)

png(paste0("./clinicalCisplatinPDUI/PDUI_cisplatin_",genename,".png"),width=3.5,height=3.7,unit="in",res=900)
ggsurvplot(surv_fit, data = best_threshold_data, 
           pval = TRUE,
           fun = "pct",
           size = 1.5,
           linetype = c("dotdash", "solid"),
           palette = c("red","blue"),
           legend = c(0.25, 0.3),
           legend.title = genename,
           surv.median.line = "hv")
dev.off()
