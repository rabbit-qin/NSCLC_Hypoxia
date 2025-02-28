rm(list = ls())
library(mistyR)
library(Seurat)
library(dplyr)
library(distances)
library(tibble)
library(purrr)
library(ggplot2)
library(cowplot)

future::plan(future::multisession,workers=8)
sample="OTAR_LNGsp10391237"
sTmisty_RDSFile=paste0("./NSCLC_STData/CNVOutput/HNRDSData/HNsTRNA(",sample,").rds")
spData=readRDS(sTmisty_RDSFile)
TumorData=subset(spData,subset=cnvcluster!="NA")
spatial_coord=GetTissueCoordinates(spData,col=c("row","col"),scale=NULL)
index=intersect(rownames(TumorData@meta.data),rownames(spatial_coord))
spData=subset(spData,cells=index)
TumorData=subset(TumorData,cells=index)
EpiCell=rownames(subset(TumorData@meta.data,cnvcluster=="Normoxia"))
MesCell=rownames(subset(TumorData@meta.data,cnvcluster=="Hypoxia"))
for(i in 1:dim(TumorData@assays$predictions@data)[1]){
  if(sum(as.double(TumorData@assays$predictions@data[i,])>6)==dim(TumorData@assays$predictions@data)[2]){
    cell=rownames(TumorData@assays$predictions@data)[i]
    break
  }
}
predictions=GetAssayData(spData,assay='predictions') %>% as.data.frame()
predictions["Normoxia",]=0
predictions["Normoxia",EpiCell]=predictions[cell,EpiCell]
predictions["Hypoxia",]=0
predictions["Hypoxia",MesCell]=predictions[cell,MesCell]
predictions[cell,c(EpiCell,MesCell)]=0
if(! any(as.double(predictions[cell,])!=0)){
  index=which(rownames(predictions)==cell)
  predictions=predictions[-index,]
}
mistyR::clear_cache()
geometry=GetTissueCoordinates(spData,col=c("row","col"),scale=NULL)
view.data=list()
features=rownames(predictions)
view.data$intra=predictions %>% as.matrix() %>% t() %>% as_tibble(rownames=NA) %>% slice(match(rownames(.),rownames(geometry))) %>% 
  rownames_to_column() %>% select(rowname,all_of(features)) %>% rename_with(make.names) %>% column_to_rownames()
view.data$juxta=view.data$intra
view.data$para=view.data$intra
spot.ids=rownames(view.data$intra)

mistyR::clear_cache()
all.views=create_initial_view(view.data$intra %>% rownames_to_column() %>% filter(rowname %in% spot.ids) %>% select(-rowname))

param=2
views.tmp=create_initial_view(view.data$juxta) %>% add_juxtaview(positions=geometry,neighbor.thr=param)
views.juxta=create_view("juxta",views.tmp[[paste0("juxtaview.",param)]]$data %>% mutate(rowname=rownames(view.data$juxta)) %>% 
                          filter(rowname %in% spot.ids) %>% select(-rowname))
all.views$juxta=views.juxta$juxta

param=3
views.tmp=create_initial_view(view.data$para) %>% add_paraview(positions=geometry,l=param)
views.para=create_view("para",views.tmp[[paste0("paraview.",param)]]$data %>% mutate(rowname=rownames(view.data$para)) %>% 
                         filter(rowname %in% spot.ids) %>% select(-rowname))
all.views$para=views.para$para

sTmisty_OutputFold=paste0("./NSCLC_STData/sTmistyOutput/",sample)
run_misty(all.views,sTmisty_OutputFold,cached=FALSE)
mistry_results=collect_results(sTmisty_OutputFold)
OutputData=mistry_results$importances %>% dplyr::select(-sample) %>% 
  dplyr::filter((Predictor == "Normoxia" | Predictor == "Hypoxia" | Target == "Normoxia" | Target == "Hypoxia") & 
                  (!(Predictor == "Normoxia" & (Target == "Hypoxia" | Target == "Normoxia"))) &
                  (!(Predictor == "Hypoxia" & (Target == "Hypoxia" | Target == "Normoxia"))) & (!(Predictor == cell | Target == cell)))
write.table(OutputData,paste0(sTmisty_OutputFold,"/misty.txt"),sep="\t",quote=FALSE,,row.names=FALSE,col.names=TRUE)

PlotData=OutputData %>% dplyr::filter(Importance > 0.5)
sTmisty_ImagePath = paste0("./NSCLC_STData/spaceRangerOutput/",sample,"/outs/spatial/tissue_hires_image.png")
img=png::readPNG(sTmisty_ImagePath)
img_grob=grid::rasterGrob(img,interpolate=FALSE,width=grid::unit(1,"npc"),height=grid::unit(1,"npc"))
spatial_coord=data.frame(spData@images[[names(spData@images)]]@coordinates)
spatial_coord$imagerow_scaled=spatial_coord$imagerow*spData@images[[names(spData@images)]]@scale.factors$hires
spatial_coord$imagecol_scaled=spatial_coord$imagecol*spData@images[[names(spData@images)]]@scale.factors$hires
spatial_coord=spatial_coord[rownames(view.data$intra),]
location=spatial_coord[,c('imagerow','imagecol')]
nnmatrix_intra=RANN::nn2(location,k=1)$nn.idx
nnmatrix_juxta=RANN::nn2(location,k=7)$nn.idx
nnmatrix_para=RANN::nn2(location,k=13)$nn.idx
for(i in 1:dim(PlotData)[1]){
  LRpair=c(PlotData[i,]$Target,PlotData[i,]$Predictor)
  expr=view.data$intra[,LRpair]
  countsum=colSums(expr)
  expr=t(log(t(expr)/countsum*median(countsum)+1))
  ligand=expr[,LRpair[1]]
  receptor=expr[,LRpair[2]]
  LRexp=rbind(ligand,receptor)
  if(PlotData[i,]$view=="intra"){
    neighexp=LRexp
  }else if(PlotData[i,]$view=="juxta"){
    neighexp=apply(nnmatrix_juxta,1,function(x){apply(LRexp[,x[2:dim(nnmatrix_juxta)[2]]],1,max)})
  }else if(PlotData[i,]$view=="para"){
    neighexp=apply(nnmatrix_para,1,function(x){apply(LRexp[,x[2:dim(nnmatrix_para)[2]]],1,max)})
  }
  LRadd=pmax(LRexp[1,]*neighexp[2,],LRexp[2,]*neighexp[1,])
  tmpLRadd=data.frame(x=spatial_coord$imagecol_scaled,y=spatial_coord$imagerow_scaled,LR=LRadd)
  tmpLRadd=tmpLRadd[c(EpiCell,MesCell),]
  if(PlotData[i,]$Predictor == "Normoxia" || PlotData[i,]$Predictor == "Hypoxia"){
    filepath=paste0(sTmisty_OutputFold,"/DEcolocalizationPlot/",PlotData[i,]$view,".",PlotData[i,]$Predictor,".",
                    PlotData[i,]$Target,".mistyR.png")
  }else{
    filepath=paste0(sTmisty_OutputFold,"/DEcolocalizationPlot/",PlotData[i,]$view,".",PlotData[i,]$Target,".",
                    PlotData[i,]$Predictor,".mistyR.png")
  }
  if(!file.exists(filepath)){
    p=ggplot()+annotation_custom(grob=img_grob,xmin=0,xmax=ncol(img),ymin=0,ymax=-nrow(img))+
      geom_point(data=tmpLRadd,aes(x=x,y=y,col=LR),size=1.5,alpha=1)+scale_y_reverse()+ylim(nrow(img),0)+xlim(0,ncol(img))+
      theme_half_open(11,rel_small=1)+theme_void()+coord_fixed(ratio=1,xlim=NULL,ylim=NULL,expand=TRUE,clip="on")+
      scale_colour_gradient2(low="white",high="red")+xlab("")+ylab("")+
      theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank())+labs(color="colocalization")+theme_minimal()+
      theme(axis.text=element_blank(),axis.title=element_blank(),panel.grid=element_blank())
    ggsave(filename=filepath,p,width = 8,height = 8,device = "png",dpi = 900,units = "in")
  }
}
