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

clinical_data <- clinical_data[order(clinical_data$pdui),]
clinical_data$signature <- "NA"
clinical_data$signature[clinical_data$pdui >= 0.84] <- "High PDUI"
clinical_data$signature[clinical_data$pdui < 0.84] <- "Low PDUI"
clinical_data$vital_status[clinical_data$vital_status == "Alive"] <- 0
clinical_data$vital_status[clinical_data$vital_status == "Dead"] <- 1
clinical_data$vital_status <- as.numeric(clinical_data$vital_status)
clinical_data_1 <- clinical_data[which(clinical_data$signature == "Low PDUI"),]

covariates <- c("alimta", "carboplatin", "cisplatin", "docetaxel", "etoposide", "gemcitabine", "navelbine", "paclitaxel", "taxol", "taxotere", "vinorelbine")                              
univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(days_to_last_follow_up,vital_status) ~', x)))
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = clinical_data_1)})

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         HR <-signif(x$coef[2], digits=2);
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res <- as.data.frame(res)
res$not.teatment <- apply(clinical_data_1[,17:27],2,function(x){sum(x)})
res$treatment <- nrow(clinical_data_1)-res$not.teatment
write.table(res,file = "./cox/Low_PDUI.csv",sep = ",",row.names = T,col.names = T,quote = F)

surv_model <- surv_fit(Surv(days_to_last_follow_up,vital_status) ~ cisplatin, data = clinical_data_1)
png("./clinicalHightLowPDUI/CARM1_Low_PDUI.png",width=3.5,height=3.7,unit="in",res=900)
ggsurvplot(surv_model, data = clinical_data_1, 
           pval = TRUE,
           fun = "pct",
           size = 1.5,
           linetype = "strata",
           palette = c("blue","red"),
           legend = c(0.3, 0.3),
           legend.title = genename,
           surv.median.line = "hv")   
dev.off()
