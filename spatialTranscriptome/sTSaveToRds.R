rm(list = ls())
library(Seurat)
library(ggplot2)
library(dplyr)

sample = "OTAR_LNGsp9476039"
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
clusters=data.frame(cluster=cutree(hc,k=2))
TumorData@meta.data$cluster=as.factor(clusters[rownames(TumorData@meta.data),])

#Plot Figures
sTcnv_EMTFile = "./NSCLC_SingleCell/NSCLCData/MSigDBData/HALLMARK_HYPOXIA.v2022.1.Hs.gmt"
genes=read.table(sTcnv_EMTFile,sep="\t",header=FALSE,check.names=FALSE)
genes=data.frame(t(genes[3:ncol(genes)]))
colnames(genes)[1]="genename"
genes=list("genes"=genes$genename[genes$genename %in% rownames(TumorData)])
EMTsetscore=AddModuleScore(object=TumorData,features=genes,name="geneset",assay = "SCT",slot = "data")
FigureData=data.frame(Group=EMTsetscore@meta.data$cluster,Value=EMTsetscore@meta.data$geneset1)
medianvalue=c()
groups=unique(FigureData$Group)
for(i in groups){
  Each=subset(FigureData,Group==i)
  medianvalue=c(medianvalue,median(Each$Value))
}
FigureData$Group=ifelse(FigureData$Group==groups[order(medianvalue)][1],"Normoxia",
                        ifelse(FigureData$Group==groups[order(medianvalue)][length(medianvalue)],"Hypoxia","Other"))
TumorData@meta.data$cnvcluster=FigureData$Group

sTcnv_OutputFold = "./NSCLC_STData/CNVOutput/HNRDSData"
saveRDS(TumorData,paste0(sTcnv_OutputFold,"/HNsTRNA(",sample,").rds"))
