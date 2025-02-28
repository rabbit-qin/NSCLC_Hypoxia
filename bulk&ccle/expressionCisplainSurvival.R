rm(list = ls())
library(survival)
library(survminer)

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

express_data <- read.csv(file = "./rawData/ExpressWithSignatureTCGA-20240724.csv",sep = ",",header = T,row.names = NULL)
genename <- "CARM1"
express_select <- express_data[which(express_data[,1] == genename),]

colnames(express_select) <- gsub("\\.","\\-",colnames(express_select))
colnames(express_select) <- sapply(colnames(express_select),function(x){substr(x,1,12)})
ids_clin <- na.omit(match(colnames(express_select),rownames(clinical_data)))
clinical <- clinical_data[ids_clin,]
ids_expr <- na.omit(match(rownames(clinical),colnames(express_select)))
clinical$expression <- as.numeric(express_select[1,ids_expr])
clinical <- clinical[-which(is.na(clinical[,3])),]

for (i in 1:nrow(clinical)) {
  if(clinical[i,10] == "Dead"){
    clinical[i,3] <- clinical[i,12]
  }
}

clinical <- clinical[which(clinical$cisplatin == 0),]
clinical$signature <- "NA"
clinical$signature[clinical$expression >= median(as.numeric(clinical$expression))] <- "High expression"
clinical$signature[clinical$expression < median(as.numeric(clinical$expression))] <- "Low expression"
clinical$vital_status[clinical$vital_status == "Alive"] <- 0
clinical$vital_status[clinical$vital_status == "Dead"] <- 1
clinical$vital_status <- as.numeric(clinical$vital_status)

surv_model <- surv_fit(Surv(days_to_last_follow_up,vital_status) ~ signature, data = clinical)

best_threshold_surv <- surv_cutpoint(clinical,
                                     time = "days_to_last_follow_up",
                                     event = "vital_status",
                                     variables = "expression",
                                     minprop = 0.4,
                                     progressbar = TRUE)
summary(best_threshold_surv)
best_threshold_data <- surv_categorize(best_threshold_surv)
surv_obj <- Surv(time = best_threshold_data$days_to_last_follow_up, event = best_threshold_data$vital_status)
surv_fit <- survfit(surv_obj ~ best_threshold_data$expression)

png(paste0("./clinicalCisplatinExpress/Express_cisplatin_",genename,".png"),width=3.5,height=3.7,unit="in",res=900)
ggsurvplot(surv_fit, data = best_threshold_data, 
           pval = TRUE,
           fun = "pct",
           size = 1.5,
           linetype = "strata",
           palette = c("blue","red"),
           legend = c(0.3, 0.3),
           legend.title = genename,
           surv.median.line = "hv")
dev.off()
