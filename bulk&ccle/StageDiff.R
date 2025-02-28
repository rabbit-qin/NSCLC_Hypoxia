rm(list = ls())
library(survival)
library(survminer)

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

genename <- "SELENBP1"
express_select <- express_data[which(express_data[,1] == genename),]
colnames(express_select) <- gsub("\\.","\\-",colnames(express_select))
colnames(express_select) <- sapply(colnames(express_select),function(x){substr(x,1,12)})
rownames(clinical) <- sapply(rownames(clinical),function(x){substr(x,1,12)})
ids_clin <- na.omit(match(unique(colnames(express_select)),rownames(clinical)))
clinical_data <- clinical[ids_clin,]
ids_expr <- na.omit(match(rownames(clinical_data),colnames(express_select)))
clinical_data$expression <- as.numeric(express_select[1,ids_expr])
clinical_data=clinical_data[clinical_data$new_stage != 'other',]

figure <- clinical_data[,c(ncol(clinical_data)-1,ncol(clinical_data))]
figure$new_stage <- factor(figure$new_stage,levels = c("s2","s1"))
mycon <- list(c("s1","s2"))
p=ggboxplot(figure,x=colnames(figure)[1],y=colnames(figure)[2],color=colnames(figure)[1],fill=colnames(figure)[2],
            size=1.5,add='nejm',palette=c('red','blue')) + 
  stat_compare_means(comparisons=mycon,method="t.test") +
  theme_bw() +
  theme(legend.position = 'none') +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
ggsave(paste0("./figure/ExpressStage/",genename,"_Express.png"),p,width = 3,height = 3.5,device = "png",dpi = 900,units = "in")
