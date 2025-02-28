rm(list = ls())
library(monocle)
library(Seurat)
library(tidyverse)

scRNAsub <- readRDS("./NSCLC_SingleCell/NSCLCData/rdsData/cancerInfercnv.rds")
table(scRNAsub@meta.data$seurat_clusters)
scRNAsub@meta.data$ident <- scRNAsub@meta.data$seurat_clusters
table(scRNAsub@meta.data$ident)

data <- as(as.matrix(scRNAsub@assays$SCT@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
mycds <- detectGenes(mycds,min_expr = 0.1)
head(fData(mycds))
head(pData(mycds))
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=3, relative_expr = TRUE)
var.genes <- VariableFeatures(scRNAsub)
mycds <- setOrderingFilter(mycds, var.genes)
p1 <- plot_ordering_genes(mycds)

expressed_genes <- row.names(subset(fData(mycds),num_cells_expressed >= 10))
diff.genes <- differentialGeneTest(mycds[expressed_genes,],fullModelFormulaStr="~seurat_clusters",cores=3) 
head(diff.genes)
deg <- subset(diff.genes, qval < 0.01)
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)
ordergene <- rownames(deg) 
mycds <- setOrderingFilter(mycds, ordergene)  
plot_ordering_genes(mycds)

mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
mycds <- orderCells(mycds)
p1 <- plot1 <- plot_cell_trajectory(mycds, color_by = "State",cell_size = 1.5,cell_link_size = 1, cell_name_size = 3, state_number_size = 4)+
  scale_color_manual(values = c("green","#FF4500","#FF1493","red","purple","yellow","blue"))+
  theme() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
ggsave("./NSCLC_SingleCell/Figure/state2.png",p1,width = 5,height = 5,device = "png",dpi = 900,units = "in")

mycds@phenoData@data[["ident"]] <- as.factor(mycds@phenoData@data[["ident"]])
p1 <- plot_cell_trajectory(mycds, color_by = "ident",cell_size = 1.5,cell_link_size = 1, cell_name_size = 3, state_number_size = 4)+
      scale_color_manual(values = c("blue","red","gray"))+
      theme() +
      theme(legend.position = 'none') +
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
      theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
ggsave("./NSCLC_SingleCell/Figure/seurat_clusters2Ref.png",p1,width = 5,height = 5,device = "png",dpi = 900,units = "in")

mycds <- orderCells(mycds,root_state = 1,num_paths = 3,reverse = T)
pData(mycds)$Pseudotime <- max(pData(mycds)$Pseudotime) - pData(mycds)$Pseudotime
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")+
  scale_color_gradient(low = "#0000FF",high = "#FF0000")+
  theme() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
ggsave("./NSCLC_SingleCell/Figure/pseudotime1Ref.png",plot3,width = 6,height = 6,device = "png",dpi = 600,units = "in")
