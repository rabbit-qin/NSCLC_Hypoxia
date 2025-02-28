rm(list=ls())
library(dplyr)
library(tidyr)
library(pROC)

Genes <- read.csv(file = "./rawData/TCGA_APAEventsSelection.csv",sep = ",",header = F,row.names = NULL)
epi_up <- Genes[3:nrow(Genes),1]
bulk_pdui <- read.csv(file = "./rawData/PDUIWithSignatureTCGA-20240724.csv",sep = ",",header = T,row.names = NULL)

bulk_pdui <- separate(bulk_pdui,row.names,into=c("num","gene","chr","trand"),sep="\\|")
rownames(bulk_pdui) <- bulk_pdui[,1]
bulk_pdui <- bulk_pdui[,-c(1,4,5)]

ids <- na.omit(match(epi_up[,1],rownames(bulk_pdui)))
epi_up_pdui <- bulk_pdui[c(ids,nrow(bulk_pdui)),]
ids_mes <- na.omit(which(bulk_pdui[nrow(bulk_pdui),] == "Hypoxia"))
ids_epi <- na.omit(which(bulk_pdui[nrow(bulk_pdui),] == "Normoxic"))
logs_pdui <- epi_up_pdui[,c(ids_epi,ids_mes)]

logs_pdui <- t(logs_pdui)
for (i in 1:(ncol(logs_pdui)-1)) {
  ids_i <- na.omit(which(is.na(logs_pdui[,i])))
  logs_pdui[ids_i,i] <- median(as.numeric(logs_pdui[,i]),na.rm = T)
}
logs_pdui <- as.data.frame(logs_pdui)
logs_pdui$Type[which(logs_pdui$Type == "Normoxic")] <- 0
logs_pdui$Type[which(logs_pdui$Type == "Hypoxia")] <- 1
colnames(logs_pdui)[ncol(logs_pdui)] <- "Type"
logs_pdui$Type <- factor(logs_pdui$Type)

library(corrplot)
cor(apply(logs_pdui[,1:(ncol(logs_pdui)-1)], 2, as.numeric))%>% corrplot(method="color")
logs_pdui[,1:(ncol(logs_pdui)-1)] <- apply(logs_pdui[,1:(ncol(logs_pdui)-1)], 2, as.numeric)
kappa(cor(apply(logs_pdui[,1:(ncol(logs_pdui)-1)], 2, as.numeric)), exact=T)

set.seed(106)    
require(caret)
folds <- createFolds(y=logs_pdui$Type,k=10) 

max=0
num=0
k=1
auc_test=c(rep(0,10))
auc_train=c(rep(0,10))
train_accuracy=c(rep(0,10))
test_accuracy=c(rep(0,10))
element_remarkable=data.frame(matrix(0,300,5))
for(i in 1:10){
  fold_test <- logs_pdui[folds[[i]],] 
  fold_train <- logs_pdui[-folds[[i]],] 
  
  print("***group***") 
  fold_pre <- glm(Type ~ NM_002081+NM_199141+NM_001040876+NM_002268+NM_004336+NM_031966+NM_001039490+NM_001242+NM_017955+NM_001202523+
                    NM_001130823+NM_022073+NM_198244+NM_002354+NM_004443+NM_005245+NM_003088+NM_001130442+NM_005348+NM_012289+
                    NM_002748+NM_130834+NM_006904+NM_005982+NM_021972+NM_001040060+NM_000358+NM_001135699+NM_006561+NM_001029989+
                    NM_003472+NM_012112+NM_004237+NM_001243142+NM_018098+NM_007027+NM_000853+NM_005562,
                  family=binomial(link='logit'),data=fold_train)
  element <- data.frame(summary(fold_pre)[["coefficients"]])
  element <- element[which(element[,4]<=0.05),]
  element$apa_event <- rownames(element)
  element_remarkable[k:(k+nrow(element)-2),] <- element[-1,]
  k=k+nrow(element)-1
  
  print(i);print("***train dataset***")
  fold_predict2 <- predict(fold_pre,type='response',newdata=fold_train)
  
  modelroc_train <- roc(as.numeric(fold_train$Type),fold_predict2)
  auc_train[i] <- auc(modelroc_train)
  plot(modelroc_train, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),print.thres=TRUE)
  
  F2<-2*(modelroc_train$sensitivities*modelroc_train$specificities)/(modelroc_train$sensitivities+modelroc_train$specificities)
  e2<-cbind(modelroc_train$thresholds,F2)
  best2<-subset(e2,e2[,2]==max(e2[,2]))[,1]
  best2<-as.numeric(best2)
  
  fold_predict2 =ifelse(fold_predict2>best2,1,0)
  fold_train$predict = fold_predict2
  fold_error2 = (as.numeric(fold_train[,(ncol(fold_train)-1)])-1)-fold_train[,ncol(fold_train)]
  fold_accuracy2 = (nrow(fold_train)-sum(abs(fold_error2)))/nrow(fold_train) 
  print(fold_accuracy2)
  
  print("***test dataset***")
  fold_predict <- predict(fold_pre,type='response',newdata=fold_test)
  
  modelroc_test <- roc(as.numeric(fold_test$Type),fold_predict)
  auc_test[i] <- auc(modelroc_test)
  plot(modelroc_test, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),print.thres=TRUE)
  
  fold_predict = ifelse(fold_predict>best2,1,0)
  fold_test$predict = fold_predict
  fold_error = (as.numeric(fold_test[,(ncol(fold_test)-1)])-1)-fold_test[,ncol(fold_test)]
  fold_accuracy = (nrow(fold_test)-sum(abs(fold_error)))/nrow(fold_test) 
  print(fold_accuracy)
  
  train_accuracy[i]=fold_accuracy2
  test_accuracy[i]=fold_accuracy
}

if(F){
  table(element_remarkable[,5])
  element_remarkable <- element_remarkable[-which(element_remarkable[,5] == 0),]
  element_select <- data.frame(table(element_remarkable[,5]))
  element_select$parameter <- rep("NA",nrow(element_select))
  for (i in 1:length(table(element_remarkable[,5]))) {
    select <- element_remarkable[which(element_remarkable[,5] == as.character(element_select[i,1])),]
    element_select$parameter[i] <- mean(select[,1])
  }
  element_select$Type <- "NA"
  element_select$Type[which(element_select$parameter > 0)] <- "Positive"
  element_select$Type[which(element_select$parameter < 0)] <- "Negative"
  colnames(element_select) <- c("APA_event","Freq","Parameter","Type")
  element_select$Parameter <- as.numeric(element_select$Parameter)
  element_select <- element_select[order(element_select$Freq),]
  element_select$Freq[which(element_select$Parameter > 0)] <- (-1)*(element_select$Freq[which(element_select$Parameter > 0)])
  parameter <- ggplot(element_select, aes(x=reorder(APA_event,Freq), y=Freq, label=Freq))+
    geom_bar(stat='identity', aes(fill=Type), width=0.7)+
    scale_fill_manual(name="Type",
                      labels = c("Positive", "Negative"),
                      values = c("Positive"="#3484BE", "Negative"="#f8766d"))+
    labs(subtitle="", title= "")+
    coord_flip()+theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
  parameter
  ggsave(filename="./NSCLC_Result02/Figure/bulk_parameter1.png",parameter,width = 5.5,height = 5.5,device = "png",dpi = 600,units = "in")
  element_select$genes <- bulk_pdui$gene[na.omit(match(element_select$APA_event,rownames(bulk_pdui)))]
}

paste0(mean(train_accuracy),"±",sd(train_accuracy))
paste0(mean(test_accuracy),"±",sd(test_accuracy))
paste0(mean(auc_train),"±",sd(auc_train))
paste0(mean(auc_test),"±",sd(auc_test))
train_accuracy;test_accuracy;auc_train;auc_test
num=10

aucValues <- c(rep(0,10))
for(k in 1:10){
  num = k
  testi <- logs_pdui[folds[[num]],]
  traini <- logs_pdui[-folds[[num]],]
  prei <- glm(Type ~ NM_002081+NM_199141+NM_001040876+NM_002268+NM_004336+NM_031966+NM_001039490+NM_001242+NM_017955+NM_001202523+
                NM_001130823+NM_022073+NM_198244+NM_002354+NM_004443+NM_005245+NM_003088+NM_001130442+NM_005348+NM_012289+
                NM_002748+NM_130834+NM_006904+NM_005982+NM_021972+NM_001040060+NM_000358+NM_001135699+NM_006561+NM_001029989+
                NM_003472+NM_012112+NM_004237+NM_001243142+NM_018098+NM_007027+NM_000853+NM_005562,
              family=binomial(link='logit'),data=traini)
  predicti <- predict.glm(prei,type='response',newdata=traini)
  
  
  train.probs <- predict(prei,type="response")
  modelroc <- roc(as.numeric(traini$Type),train.probs)
  F3<-2*(modelroc$sensitivities*modelroc$specificities)/(modelroc$sensitivities+modelroc$specificities)
  e3<-cbind(modelroc$thresholds,F3)
  best3<-subset(e3,e3[,2]==max(e3[,2]))[,1]
  best3<-as.numeric(best3)

  plot(modelroc,                    
       print.auc=TRUE,                   
       print.auc.x=0.5,print.auc.y=0.5,  
       auc.polygon=TRUE,                 
       auc.polygon.col="skyblue",        
       grid= FALSE,                      
       legacy.axes=TRUE)                 
  
  predicti = ifelse(predicti>best3,1,0)
  traini$predict = predicti
  errori = (as.numeric(traini[,(ncol(traini)-1)])-1)-traini[,ncol(traini)]
  accuracyi = (nrow(traini)-sum(abs(errori)))/nrow(traini)
  
  testi_predict2 <- predict(prei,type='response',newdata=testi)
  modelroc_test <- roc(as.numeric(testi$Type),testi_predict2)
  plot(modelroc_test,                    
       print.auc=TRUE,                   
       print.auc.x=0.5,print.auc.y=0.5,  
       auc.polygon=TRUE,                 
       auc.polygon.col="skyblue",        
       grid= FALSE,                      
       legacy.axes=TRUE)                 
  if(F){
    plot.roc(modelroc_test,
             lty = 2,
             lwd = 3,
             add=TRUE,                       
             col = "red",                     
             #auc.polygon.col="green",        
             #print.auc=TRUE,                 
             print.auc.col = "red",           
             print.auc.x=0.5,print.auc.y=0.4) 
    legend(0.35,0.30,                         
           bty = "n",                        
           title="", 
           legend=c("",""), 
           col=c("#252A34","#FF2E63"), 
           lwd=2)  
  }
  
  predicti2 <- predict.glm(prei,type='response',newdata=testi)
  predicti2 =ifelse(predicti2>best3,1,0)
  testi$predict = predicti2
  errori2 = (as.numeric(testi[,(ncol(testi)-1)])-1)-testi[,ncol(testi)]
  accuracyi2 = (nrow(testi)-sum(abs(errori2)))/nrow(testi) 
  
  accuracyi;num;accuracyi2
  logit.step <- step(prei, direction = "both")
  summary(logit.step)
  summary(prei)
  element <- data.frame(coefficients(prei))
  
  library(car)
  vif(prei) 
  
  pdui_ccle <- read.csv(file = "./rawData/PDUIWithSignatureCCLE-20240724.csv",sep = ",",header = T,row.names = NULL)
  pdui_ccle <- separate(pdui_ccle,row.names,into=c("num","gene","chr","trand"),sep="\\|")
  rownames(pdui_ccle) <- pdui_ccle[,1]
  pdui_ccle <- pdui_ccle[,-c(1,3,4,5)]
  ids_mes <- na.omit(which(bulk_pdui[nrow(bulk_pdui),] == "Hypoxia"))
  ids_epi <- na.omit(which(bulk_pdui[nrow(bulk_pdui),] == "Normoxic"))
  
  Epi_pdui_ccle <- pdui_ccle[-nrow(pdui_ccle),na.omit(which(pdui_ccle[nrow(pdui_ccle),] == "Normoxic"))]
  Mes_pdui_ccle <- pdui_ccle[-nrow(pdui_ccle),na.omit(which(pdui_ccle[nrow(pdui_ccle),] == "Hypoxia"))]
  
  ids <- na.omit(match(epi_up[,1],rownames(pdui_ccle)))
  Epi_pdui_ccle <- Epi_pdui_ccle[ids,]
  Mes_pdui_ccle <- Mes_pdui_ccle[ids,]
  pdui_ccle <- as.data.frame(t(cbind(Epi_pdui_ccle,Mes_pdui_ccle)))
  pdui_ccle$TYPE <- c(rep(0,ncol(Epi_pdui_ccle)),rep(1,ncol(Mes_pdui_ccle)))
  pdui_ccle$TYPE <- factor(pdui_ccle$TYPE)
  
  for (i in 1:(ncol(pdui_ccle)-1)) {
    ids_i <- na.omit(which(is.na(pdui_ccle[,i])))
    pdui_ccle[ids_i,i] <- median(as.numeric(pdui_ccle[,i]),na.rm = T)
  }
  
  cor(apply(pdui_ccle[,1:(ncol(pdui_ccle)-1)], 2, as.numeric))%>% corrplot(method="color")
  pdui_ccle[,1:(ncol(pdui_ccle)-1)] <- apply(pdui_ccle[,1:(ncol(pdui_ccle)-1)], 2, as.numeric)
  predicti2_ccle <- predict.glm(prei,type='response',newdata=pdui_ccle)
  modelroc <- roc(pdui_ccle$TYPE,predicti2_ccle)
  
  plot(modelroc,                         
       print.auc=TRUE,                   
       print.auc.x=0.5,print.auc.y=0.5,  
       auc.polygon=TRUE,                 
       auc.polygon.col="skyblue",        
       grid= FALSE,                      
       legacy.axes=TRUE)                 
  aucValues[k] <- as.numeric(modelroc$auc)
}
aucValues
paste0(mean(aucValues),"±",sd(aucValues))
