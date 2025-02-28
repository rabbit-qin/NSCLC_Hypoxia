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
express_data <- read.csv(file = "./rawData/ExpressWithSignatureTCGA-20240724.csv",sep = ",",header = T,row.names = NULL)
express_data <- express_data[,-which(express_data[nrow(express_data),] == "Mix")]

clinical$new_stage = ifelse(clinical$ajcc_pathologic_stage %in% c( 'Stage I','Stage IA','Stage IB'),'s1',
                            ifelse(clinical$ajcc_pathologic_stage %in% c('Stage II' ,'Stage IIA','Stage IIB','Stage IIC'),'s2',
                                   ifelse(clinical$ajcc_pathologic_stage %in% c('Stage III','Stage IIIA','Stage IIIB','Stage IIIC'),'s2',
                                          ifelse(clinical$ajcc_pathologic_stage %in% c( 'Stage IV','Stage IVA' ,'Stage IVB','Stage IVC'),'s2','other'
                                          ) ) ) )
table(clinical$new_stage)

colnames(bulk_select) <- gsub("\\.","\\-",colnames(bulk_select))
colnames(bulk_select) <- sapply(colnames(bulk_select),function(x){substr(x,1,12)})
rownames(clinical) <- sapply(rownames(clinical),function(x){substr(x,1,12)})
ids_clin <- na.omit(match(unique(colnames(bulk_select)),rownames(clinical)))
clinical_data <- clinical[ids_clin,]
ids_expr <- na.omit(match(rownames(clinical_data),colnames(bulk_select)))
clinical_data$expression <- as.numeric(bulk_select[1,ids_expr])
clinical_data=clinical_data[clinical_data$new_stage != 'other',]

figure <- clinical_data[,c(ncol(clinical_data)-1,ncol(clinical_data))]

figure$new_stage <- factor(figure$new_stage,levels = c("s2","s1"))
figure <- figure[order(figure$new_stage,decreasing = TRUE),]
mycon <- list(c("s1","s2"))
p=ggboxplot(figure,x=colnames(figure)[1],y=colnames(figure)[2],color=colnames(figure)[1],fill=colnames(figure)[2],
            size=1.5,add='nejm',palette=c('red','blue')) + #jitter
  stat_compare_means(comparisons=mycon,method="t.test") +
  ylim(c(-3,2.8)) + #,label.y=max(as.numeric(figure[,1]))+0.01
  theme_bw() +
  theme(legend.position = 'none') +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
p

ggsave(paste0("./NSCLC_Result02/Figure/Scores_Stage.png"),p,width = 4.5,height = 4.5,device = "png",dpi = 1200,units = "in")
