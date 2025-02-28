rm(list = ls())
library(Seurat)
library(SPATA2)
library(ggpubr)
library(png)
library(grid)
library(cowplot)

Read10X_HiresImage <- function(image.dir, filter.matrix = TRUE, ...) {
  library(png)
  library(jsonlite)
  image <- readPNG(source = file.path(image.dir, 'tissue_hires_image.png'))
  scale.factors <- fromJSON(txt = file.path(image.dir, 'scalefactors_json.json'))
  tissue.positions.path <- Sys.glob(paths = file.path(image.dir, 'tissue_positions*'))
  tissue.positions <- read.csv( file = tissue.positions.path,
                                col.names = c('barcodes', 'tissue', 'row', 'col', 'imagerow', 'imagecol'),
                                header = ifelse( test = basename(tissue.positions.path) == "tissue_positions.csv",
                                                 yes = TRUE,
                                                 no = FALSE
                                ),
                                as.is = TRUE,
                                row.names = 1)
  if (filter.matrix) {
    tissue.positions <- tissue.positions[which(x = tissue.positions$tissue == 1), , drop = FALSE]
  }
  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_hires_scalef
  spot.radius <-  unnormalized.radius / max(dim(x = image))
  return(new(
    Class = 'VisiumV1',
    image = image,
    scale.factors = scalefactors(
      spot = scale.factors$spot_diameter_fullres,
      fiducial = scale.factors$fiducial_diameter_fullres,
      hires = scale.factors$tissue_hires_scalef,
      scale.factors$tissue_lowres_scalef
    ),
    coordinates = tissue.positions,
    spot.radius = spot.radius
  ))
}

sample = "OTAR_LNGsp10206165"
SampleFold=paste0('./NSCLC_STData/spaceRangerOutput/',sample,'/outs/')
sTcnv_OutputFold = "./NSCLC_STData/CNVOutput/RDSData"
Image=Read10X_HiresImage(paste0(SampleFold,"spatial/"))
Data=Load10X_Spatial(SampleFold,image = Image)
Data@images$slice1@scale.factors$lowres=Data@images$slice1@scale.factors$hires
Data[["percent.mt"]] <- PercentageFeatureSet(Data,pattern = "^MT-")
Data@meta.data$percent.mt[which(Data@meta.data$percent.mt=="NaN")]=0
Data[["percent.rb"]] <- PercentageFeatureSet(Data, pattern = "^RP[SL]")
Data@meta.data$percent.rb[which(Data@meta.data$percent.rb=="NaN")]=0
Data[["percent.HB"]]<-PercentageFeatureSet(Data,features=CaseMatch(c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ"), 
                                                                   rownames(Data)))
Data@meta.data$percent.HB[which(Data@meta.data$percent.HB=="NaN")]=0
saveRDS(Data,paste0(sTcnv_OutputFold,"/sTRNA(",sample,").rds"))

sTcnv_CellFraction = paste0("./NSCLC_STData/cell2LocationOutput/outputCSV/",sample,".csv")
sTcnv_OutputFold = paste0("./NSCLC_STData/CNVOutput/sTcnvOutput/",sample,"/")

cell2location=read.table(sTcnv_CellFraction,sep=",",header=TRUE,check.names = F,row.names=1)
predictions.assay = CreateAssayObject(t(cell2location))
Data[["predictions"]] <- predictions.assay
Data=SCTransform(Data, assay = "Spatial", verbose = FALSE, variable.features.n = 3000,
                 vars.to.regress = "percent.mt",return.only.var.genes = TRUE,vst.flavor = "v2")
Data <- ScaleData(object = Data,do.scale = FALSE,do.center = FALSE)
maxcell=data.frame(cell=Data@assays$predictions@data["Cancer",]>6)
Data@meta.data$cell=maxcell[rownames(Data@meta.data),"cell"]
TumorData=subset(Data,subset=cell=="TRUE")
spata2Data=asSPATA2(TumorData,sample_name = Images(TumorData), spatial_method = "Visium",assay_name = "SCT",image_name = Images(TumorData))
spata2Data=runCnvAnalysis(object=spata2Data,directory_cnv_folder=sTcnv_OutputFold,cnv_prefix="Chr")
cnv_results=getCnvResults(spata2Data)
saveRDS(cnv_results,paste0(sTcnv_OutputFold,"/cnv.results.rds"))

hc=hclust(dist(t(cnv_results$cnv_mtr),method="euclidean"),method="ward.D2")
plot(hc)
clusters=data.frame(cluster=cutree(hc,k=2))
TumorData@meta.data$cluster=as.factor(clusters[rownames(TumorData@meta.data),])

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
FigureData$Group=factor(FigureData$Group,levels=c("Normoxia","Other","Hypoxia"))
num=1
InforUniGroup=unique(FigureData$Group)
for(j in 1:(length(InforUniGroup)-1)){
  for(k in (j+1):length(InforUniGroup)){
    if(num==1){
      my_comparisons=list(c(as.character(InforUniGroup[j]),as.character(InforUniGroup[k])))
    }else{
      my_comparisons[[num]]=c(as.character(InforUniGroup[j]),as.character(InforUniGroup[k]))
    }
    num=num+1
  }
}

meta_cols <- c('Normoxia'='blue','Hypoxia'='red','Other'='#B2B2B2')
p <- ggboxplot(FigureData, x = "Group", y = "Value",color = "Group", palette = "rainbow12",size=1.8,add='nejm')+
  stat_compare_means(comparisons=my_comparisons)+border()+
  scale_color_manual(values=meta_cols)+ylab("HYPOXIA score")+
  theme(legend.position="none",plot.title=element_text(size=0,hjust=0.5,face = "bold"),
        axis.text=element_text(face = "bold",hjust=0.5),axis.title.x=element_text(size=0),
        axis.title.y=element_text(face = "bold"))+
  theme_bw() +
  theme(legend.position = 'none') +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
p
ggsave(filename=paste0(sTcnv_OutputFold,"/figure/cnv.cluster.Hypoxia.png"),p,width = 8,height = 8,device = "png",dpi = 900,units = "in")

slice=names(TumorData@images)[1]
SliceLocation=TumorData@images[[slice]]@coordinates
SliceLocation$imagerow_scaled = SliceLocation$imagerow*TumorData@images[[slice]]@scale.factors$hires
SliceLocation$imagecol_scaled = SliceLocation$imagecol*TumorData@images[[slice]]@scale.factors$hires
SliceLocation$Cell=TumorData@assays$predictions@data["Cancer",rownames(SliceLocation)]
SliceLocation$Cell=SliceLocation$Cell/30
SliceLocation$cluster=TumorData@meta.data[rownames(SliceLocation),]$cluster
SliceLocation$updatedcluster=ifelse(SliceLocation$cluster==groups[order(medianvalue)][1],"Normoxia",
                                    ifelse(SliceLocation$cluster==groups[order(medianvalue)][length(medianvalue)],"Hypoxia","Other"))

sTcnv_ImagePath = paste0("./NSCLC_STData/spaceRangerOutput/",sample,"/outs/spatial/tissue_hires_image.png")
img=readPNG(sTcnv_ImagePath)
img_grob=rasterGrob(img,interpolate=FALSE,width=grid::unit(1,"npc"),height=grid::unit(1,"npc"))
p1=ggplot()+annotation_custom(grob=img_grob,xmin=0,xmax=ncol(img),ymin=0,ymax=-nrow(img))+
  geom_point(data=SliceLocation,aes(x=imagecol_scaled,y=imagerow_scaled,size=Cell,alpha=Cell),
             color="blue")+
  scale_y_reverse()+ylim(nrow(img),0)+xlim(0,ncol(img))+theme_half_open(11,rel_small=1)+theme_void()+
  coord_fixed(ratio=1,xlim=NULL,ylim=NULL,expand=TRUE,clip="on")+
  scale_size_continuous(range=c(0,1))+scale_alpha_continuous(range=c(0,1))+labs(size=SliceLocation$Cell)+guides(alpha="none")+
  theme(legend.key.size=grid::unit(50,'pt'),legend.title=element_text(size=0),legend.text=element_text(size=12),legend.position="bottom")
p1
ggsave(filename=paste0(sTcnv_OutputFold,"/figure/spatial.tumor01.png"),p1,width = 8,height = 8,device = "png",dpi = 900,units = "in")

p2=ggplot()+annotation_custom(grob=img_grob,xmin=0,xmax=ncol(img),ymin=0,ymax=-nrow(img))+
  geom_point(data=SliceLocation,aes(x=imagecol_scaled,y=imagerow_scaled,colour=updatedcluster),size=1.1,alpha=1)+
  scale_y_reverse()+ylim(nrow(img),0)+xlim(0,ncol(img))+theme_half_open(11,rel_small=1)+theme_void()+
  coord_fixed(ratio=1,xlim=NULL,ylim=NULL,expand=TRUE,clip="on")+
  scale_color_manual(values=meta_cols)+guides(color=guide_legend(override.aes=list(size=5)))+labs(color=SliceLocation$updatedcluster)+
  theme(legend.key.size=grid::unit(50,'pt'),legend.title=element_text(size=0),legend.text=element_text(size=12),legend.position="bottom")
p2
ggsave(filename=paste0(sTcnv_OutputFold,"/figure/spatial.tumor02.png"),p2,width = 8,height = 8,device = "png",dpi = 900,units = "in")
