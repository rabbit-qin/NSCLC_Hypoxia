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
bulk_select <- as.data.frame(bulk_select[nrow(bulk_select),7:ncol(bulk_select)])

clinical <- read.csv(file = "./rawData/NSCLC_clinical.CSV",sep = ",",header = T,row.names = 1)
clinical <- clinical[which(as.character(clinical$treatments_pharmaceutical_treatment_or_therapy) == "no"),]
clinical <- clinical[which(as.character(clinical$treatments_radiation_treatment_or_therapy) == "no"),]

colnames(bulk_select) <- gsub("\\.","\\-",colnames(bulk_select))
colnames(bulk_select) <- sapply(colnames(bulk_select),function(x){substr(x,1,12)})
rownames(clinical) <- sapply(rownames(clinical),function(x){substr(x,1,12)})
ids_clin <- na.omit(match(unique(colnames(bulk_select)),rownames(clinical)))
clinical_data <- clinical[ids_clin,]
ids_expr <- na.omit(match(rownames(clinical_data),colnames(bulk_select)))
clinical_data$Scroes <- as.numeric(bulk_select[1,ids_expr])

clinical_data <- clinical_data[-which(is.na(clinical_data[,3])),]

for (i in 1:nrow(clinical_data)) {
  if(clinical_data[i,10] == "Dead"){
    clinical_data[i,3] <- clinical_data[i,12]
  }
}

clinical_data$signature <- "NA"
clinical_data$signature[clinical_data$Scroes >= median(as.numeric(clinical_data$Scroes))] <- "High scroes"
clinical_data$signature[clinical_data$Scroes < median(as.numeric(clinical_data$Scroes))] <- "Low scroes"

clinical_data$vital_status[clinical_data$vital_status == "Alive"] <- 0
clinical_data$vital_status[clinical_data$vital_status == "Dead"] <- 1

clinical_data$vital_status <- as.numeric(clinical_data$vital_status)

best_threshold_surv <- surv_cutpoint(clinical_data,
                                     time = "days_to_last_follow_up",  
                                     event = "vital_status",           
                                     variables = "Scroes",             
                                     minprop = 0.4,                    
                                     progressbar = TRUE)              
summary(best_threshold_surv)

best_threshold_data <- surv_categorize(best_threshold_surv)

surv_obj <- Surv(time = best_threshold_data$days_to_last_follow_up, event = best_threshold_data$vital_status)

surv_fit <- survfit(surv_obj ~ best_threshold_data$Scroes)

png(paste0("./NSCLC_Result02/Figure/Scores_Survival.png"),width=5,height=5,unit="in",res=1200)
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
           linetype = c("dotted","solid"), #"strata", 
           palette = c("red","mediumblue"), 
           legend = c(0.25, 0.2), 
           #legend.title = genename,  
           #risk.table = TRUE,
           #legend.labs = c("High PDUI", "Low PDUI"), 
           surv.median.line = "hv",
           #ylim=c(15,100),xlim=c(0,8000),
           ggtheme = theme_temp)   
dev.off()
