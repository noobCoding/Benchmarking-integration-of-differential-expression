# library(Seurat)
# install.packages('remotes')
# remotes::install_version("SDMTools", "1.1-221")
# remotes::install_version("multtest")
# BiocManager::install('devtools')
library(ggplot2)
library(Seurat)

lsdir <- list.dirs('data', recursive=FALSE)

# pl<-lapply(lsdir,function(x){
  # read data counts and cellinfo
x<-lsdir[3]
  counts <- read.table(paste0(x,'/counts.txt'), header = TRUE, sep='\t', fill = TRUE)
  cellinfo <- read.table(paste0(x,'/cellinfo.txt'), head=T, sep='\t', fill = TRUE)
  
  seu <- CreateSeuratObject(counts)
  seu <- SCTransform(seu)
  
  # Add cell type annotation to metadata
  # seu <- AddMetaData(seu, as.character(cellinfo$Batch), col.name = "Batch")
  seu <- AddMetaData(seu, cellinfo$Batch, col.name = "Batch")
  seu <- AddMetaData(seu, cellinfo$Group, col.name = "Group")

  # PCA
  seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
  # TSNE
  seu <- RunTSNE(seu, dims = 1:30, seed.use = 7968)
  # UMAP
  seu <- RunUMAP(seu, dims = 1:30, seed.use = 7968, uwot.sgd = TRUE)
  # ElbowPlot(seu, ndims = 30)
  
  DimPlot(seu, reduction = "pca", 
          group.by = "Group", shape.by="Batch", pt.size = 2, label = F, repel = T) +
    scale_color_brewer(type = "qual", palette = "Set1") + scale_shape_identity() + 
    theme(plot.title = element_blank())
 
  DimPlot(seu, reduction = "tsne",
          group.by = "Group", shape.by="Batch", pt.size = 2,label = F, repel = TRUE) +
    scale_color_brewer(type = "qual", palette = "Set1")+ scale_shape_identity()

  DimPlot(seu, reduction = "umap",
          group.by = "Group", shape.by="Batch", pt.size = 2,label = F, repel = TRUE) +
    scale_color_brewer(type = "qual", palette = "Set1")+ scale_shape_identity()
# })

# library(gridExtra)
# library(grid)
# library(ggplot2)
# library(lattice)
# 
# png(filename = 'group_vs_batch_distribution.png', units = 'px', width = 1000, height = 600)
# 
# tg <- textGrob('Group vs Batch data distribution', gp = gpar(fontsize = 16, fontface = 'bold'))
# # sg <- textGrob('', gp = gpar(fontsize = 10))
# margin <- unit(0.5, "line")
# grided <- gridExtra::grid.arrange(grobs = pl, ncol = 3)
# gridExtra::grid.arrange(tg, grided,
#                         heights = unit.c(grobHeight(tg) + 1.2*margin,
#                                          # grobHeight(sg) + margin,
#                                          unit(1,"null")))
# 
# 
# dev.off()