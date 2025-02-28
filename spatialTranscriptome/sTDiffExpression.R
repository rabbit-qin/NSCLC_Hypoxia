rm(list = ls())
library(Seurat)
library(SPATA2)
library(ggpubr)
library(png)
library(grid)
library(cowplot)

sample = "OTAR_LNGsp10206157"
gene = "SELENBP1"
Data = readRDS(paste0("./NSCLC_STData/CNVOutput/RDSData/sTRNA(",sample,").rds"))
cnv_results = readRDS(paste0("./NSCLC_STData/CNVOutput/sTcnvOutput/",sample,"/cnv.results.rds"))
sTcnv_CellFraction = paste0("./NSCLC_STData/cell2LocationOutput/outputCSV/",sample,".csv")
cell2location=read.table(sTcnv_CellFraction,sep=",",header=TRUE,check.names = F,row.names=1)

predictions.assay = CreateAssayObject(t(cell2location))
Data[["predictions"]] <- predictions.assay
Data=SCTransform(Data, assay = "Spatial", verbose = FALSE, variable.features.n = 3000,
                 vars.to.regress = "percent.mt",return.only.var.genes = TRUE,vst.flavor = "v2")
Data <- ScaleData(object = Data,do.scale = FALSE,do.center = FALSE)
maxcell=data.frame(cell=Data@assays$predictions@data["Cancer",]>6)
Data@meta.data$cell=maxcell[rownames(Data@meta.data),"cell"]
TumorData=subset(Data,subset=cell=="TRUE")

hc=hclust(dist(t(cnv_results$cnv_mtr),method="euclidean"),method="ward.D2")
plot(hc)
clusters=data.frame(cluster=cutree(hc,k=3))
TumorData@meta.data$cluster=as.factor(clusters[rownames(TumorData@meta.data),])

sTcnv_EMTFile = "./MSigDBData/HALLMARK_HYPOXIA.v2022.1.Hs.gmt"
genes=read.table(sTcnv_EMTFile,sep="\t",header=FALSE,check.names=FALSE)
genes=data.frame(t(genes[3:ncol(genes)]))
colnames(genes)[1]="genename"
genes=list("genes"=genes$genename[genes$genename %in% rownames(TumorData)])
EMTsetscore=AddModuleScore(object=TumorData,features=genes,name="geneset",assay = "SCT",slot = "data")
FigureData=data.frame(Group=EMTsetscore@meta.data$cluster,Value=EMTsetscore@meta.data$geneset1)

meta_cols <- c('Normoxia'='blue','Hypoxia'='red','Other'='#B2B2B2')
my_comparisons <- list(c("1", "2"),c("2", "3"),c("1", "3"))
p0 <- ggboxplot(FigureData, x = "Group", y = "Value",color = "Group", palette = "rainbow12",size=1.8,add='nejm')+
  stat_compare_means(comparisons=my_comparisons)+border()+
  scale_color_manual(values=meta_cols)+ylab("HYPOXIA score")+
  theme(legend.position="none",plot.title=element_text(size=0,hjust=0.5,face = "bold"),
        axis.text=element_text(face = "bold",hjust=0.5),axis.title.x=element_text(size=0),
        axis.title.y=element_text(face = "bold"))+
  theme_bw() +
  theme(legend.position = 'none') +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
p0
ggsave(filename=paste0("./NSCLC_STData/CNVOutput/DiffExprFigure/",sample,"/sample1.png"),p0,width = 8,height = 8,device = "png",dpi = 900,units = "in")

levels(FigureData$Group) <- c(levels(FigureData$Group), "Hypoxia","Other","Normoxia")
FigureData$Group[which(FigureData$Group == "1")]= "Hypoxia"
FigureData$Group[which(FigureData$Group == "2")]= "Normoxia"
FigureData$Group[which(FigureData$Group == "3")]= "Other"
my_comparisons <- list(c("Hypoxia", "Normoxia"),c("Hypoxia", "Other"),c("Other", "Normoxia"))
p1 <- ggboxplot(FigureData, x = "Group", y = "Value",color = "Group", palette = "rainbow12",size=1.8,add='nejm')+
  stat_compare_means(comparisons=my_comparisons)+border()+
  scale_color_manual(values=meta_cols)+ylab("HYPOXIA score")+
  theme(legend.position="none",plot.title=element_text(size=0,hjust=0.5,face = "bold"),
        axis.text=element_text(face = "bold",hjust=0.5),axis.title.x=element_text(size=0),
        axis.title.y=element_text(face = "bold"))+
  theme_bw() +
  theme(legend.position = 'none') +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
p1
ggsave(filename=paste0("./NSCLC_STData/CNVOutput/DiffExprFigure/",sample,"/sample2.png"),p1,width = 8,height = 8,device = "png",dpi = 900,units = "in")

TumorData@meta.data$cnvcluster=FigureData$Group

gene = "SELENBP1"
CARM1_expression <- GetAssayData(TumorData, assay = "SCT", slot = "counts")[gene, ]
cnvcluster <- TumorData@meta.data[["cnvcluster"]]
CARM1_data <- data.frame( CARM1_expression = CARM1_expression, cnvcluster = cnvcluster)
colnames(CARM1_data) <- c("Value","Group")
mean(CARM1_data$Value[which(CARM1_data$Group == "Hypoxia")])
t.test(CARM1_data$Value[which(CARM1_data$Group == "Hypoxia")],CARM1_data$Value[which(CARM1_data$Group == "Normoxia")])
table(CARM1_data$Group)
table(CARM1_data[which(CARM1_data$Group == "Hypoxia"),1])
table(CARM1_data[which(CARM1_data$Group == "Normoxia"),1])

meta_cols <- c('Normoxia'='blue','Hypoxia'='red')
my_comparisons <- list(c("Normoxia", "Hypoxia"))
p2 <- ggboxplot(CARM1_data, x = "Group", y = "Value",color = "Group", palette = "rainbow12",size=1.8,add='nejm')+
  stat_compare_means(comparisons=my_comparisons)+border()+
  scale_color_manual(values=meta_cols)+
  theme(legend.position="none",plot.title=element_text(size=0,hjust=0.5,face = "bold"),
        axis.text=element_text(face = "bold",hjust=0.5),axis.title.x=element_text(size=0),
        axis.title.y=element_text(face = "bold"))+
  theme_bw() +
  theme(legend.position = 'none') +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
p2
ggsave(filename=paste0("./NSCLC_STData/CNVOutput/DiffExprFigure/",sample,"/",gene,"/gene.png"),p2,width = 8,height = 8,device = "png",dpi = 900,units = "in")

TumorData_hpy = subset(TumorData,subset=cnvcluster=="Hypoxia")
TumorData_nor = subset(TumorData,subset=cnvcluster=="Normoxia")
range_hpy <- range(GetAssayData(TumorData_hpy, slot = "data")[gene, ])
range_nor <- range(GetAssayData(TumorData_nor, slot = "data")[gene, ])
common_range <- range(range_hpy, range_nor)
color_map <- c("#D8BFD8", "lightcoral", "red", "darkred")
p3 <- SpatialFeaturePlot(TumorData_hpy, features = gene, pt.size.factor = 1.5)
p3 <- p3 + scale_fill_gradientn(colors = color_map, limits = common_range)

p4 <- SpatialFeaturePlot(TumorData_nor, features = gene, pt.size.factor = 1.5)
p4 <- p4 + scale_fill_gradientn(colors = color_map, limits = common_range)
p3+p4
ggsave(filename=paste0("./NSCLC_STData/CNVOutput/DiffExprFigure/",sample,"/",gene,"/DimplotNormoxia.png"),p4,width = 8,height = 8,device = "png",dpi = 900,units = "in")
