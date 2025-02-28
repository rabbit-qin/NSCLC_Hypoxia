rm(list = ls())

expession_data <- read.csv(file = "./rawData/expression_TCGA.csv",sep = ",",row.names = 1,header = T)
hypoxia_geneset <- read.csv(file = "./NMFScore/NMF_geneset.csv",sep = ",",header = T)

gene <- intersect(hypoxia_geneset$NMF.marker.genes,rownames(expession_data))
hypoxiaExp <- expession_data[gene,]
dim(hypoxiaExp)
hypoxiaExp <- hypoxiaExp[-which(apply(hypoxiaExp,1,mean)<1),]
sum(is.na(hypoxiaExp))

library(NMF)
log2fpkm.hypoxia <- log2(hypoxiaExp+1)
ranks <- 2:10
estim <- nmf(log2fpkm.hypoxia,ranks, nrun=25)
duplicated(colnames(log2fpkm.hypoxia))
png(paste0("./NMFScore/figure/clusterInformationTcga.png"),width=5,height=5,unit="in",res=900)
plot(estim)
dev.off()

seed = 7
nmf.rank4 <- nmf(log2fpkm.hypoxia, 
                 rank = 2, 
                 nrun=25,
                 seed = seed, 
                 method = "brunet")
jco <- c("#2874C5","#EABF00","#C6524A","#868686")
index <- extractFeatures(nmf.rank4,"max") 
sig.order <- unlist(index)
NMF.Exp.rank4 <- log2fpkm.hypoxia[sig.order,]
NMF.Exp.rank4 <- na.omit(NMF.Exp.rank4)
group <- predict(nmf.rank4)
table(group)

png(paste0("./NMFScore/figure/cluster2types.png"),width=5,height=5,unit="in",res=900)
consensusmap(nmf.rank4,
             labRow = NA,
             labCol = NA,
             annCol = data.frame("cluster"=group[colnames(NMF.Exp.rank4)]),
             annColors = list(cluster=c("1"=jco[1],"2"=jco[2],"3"=jco[3],"4"=jco[4])))
dev.off()

library(GSVA)
library(dplyr)
library(GSEABase)
group <- data.frame(group)

gene_exp <- apply(expession_data, 2, as.numeric)
rownames(gene_exp) <- rownames(expession_data)
other_set=getGmt("./MsigDB/HALLMARK_HYPOXIA.v2022.1.Hs.gmt")
Es=gsva(expr=as.matrix(gene_exp),gset.idx.list=other_set,kcdf="Gaussian",parallel.sz=4,method="ssgsea")  #ssgsea
Es <- data.frame(t(Es))
Es$type <- "NA"
ids <- na.omit(match(rownames(Es),rownames(group)))
Es[,2] <- group[ids,1]

library(ggplot2)
library(ggpubr)
mycon <- list(c('2','1'))
picture2 <- ggviolin(Es,x = 'type',y = 'HALLMARK_HYPOXIA',color = 'type',fill = 'type',add.params = list(fill="white",size=1),
                     size = 2,add='boxplot',palette=c('red','blue'),alpha = 0.5) + 
            stat_compare_means(comparisons = mycon,method = "t.test",paired=F,bracket.size=0.6,step.increase=0.08,label.y = 2.3) + 
            ylim(1,2.5)+
            theme_bw() +
            theme(legend.position = 'none') +
            theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
            theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
ggsave("./NMFScore/figure/violin2.png",picture2,width = 3,height = 3.5,device = "png",dpi = 900,units = "in")

group$type <- "NA"
group$type[which(group$group == '1')] <- "hypoxia"
group$type[which(group$group == '2')] <- "normoxic"
table(group$type)
write.table(group,file = "./NMFScore/cluster2Tcga.csv",quote = FALSE,sep=",",row.names=T,col.names=TRUE)
