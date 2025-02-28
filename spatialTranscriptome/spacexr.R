rm(list = ls())
library(arrow)
library(Seurat)
library(spacexr)

FlexOutPath <- "./NSCLCWorks/HDData/cellrangerOut" # Path to cellranger aggr output folder
ColonFlex.data <- -Read10X_h5(paste0(FlexOutPath,"/filtered_feature_bc_matrix.h5")) 

FlexRef<-Read10X_h5(paste0(FlexOutPath,"/filtered_feature_bc_matrix.h5"))
MetaData<-readRDS('./NSCLCWorks/HDData/cellrangerOut/metadataNew.rds') 

KpIdents<-names(which(table(MetaData$ident)>25))
MetaData<-MetaData[MetaData$ident%in%KpIdents,]
MetaData$Barcode <- rownames(MetaData)
FlexRef<-FlexRef[,MetaData$Barcode]

CTRef<-MetaData$ident
CTRef<-gsub("/","_",CTRef)
CTRef<-as.factor(CTRef)
names(CTRef)<-MetaData$Barcode

reference <- Reference(FlexRef[,names(CTRef)], CTRef , colSums(FlexRef))
counts<-Read10X_h5("./NSCLCWorks/HDData/HDData/outsHD3/binned_outputs/square_008um/filtered_feature_bc_matrix.h5")
coords<-read_parquet("./NSCLCWorks/HDData/HDData/outsHD3/binned_outputs/square_008um/spatial/tissue_positions.parquet",as_data_frame = TRUE)
rownames(coords)<-coords$barcode
coords<-coords[colnames(counts),]
rownames(coords)<-coords$barcode
coords_1<-coords
coords<-coords[,3:4]
rownames(coords)<-coords_1$barcode
nUMI <- colSums(counts)

puck <- SpatialRNA(coords, counts, nUMI)
barcodes <- colnames(puck@counts)
myRCTD <- create.RCTD(puck, reference, max_cores = 10)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')

saveRDS(myRCTD,file="./NSCLCWorks/HDData/DeconvolutionResults/LUADHD3_Deconvolution_HD_New.rds")