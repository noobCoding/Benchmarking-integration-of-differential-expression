rm(list=ls())
# BiocManager::install('Seurat')
library(Seurat)  # Seurat v3 version
library(cowplot)
library(ggplot2)
# 

lsdir <- list.dirs('data', recursive=FALSE) 
dir.create('demo_seurat3')
sapply(lsdir,function(x){
  
  x2 <- gsub('data/','',x)

  dir.create(paste0('demo_seurat3/',x2), showWarnings = FALSE)
  
  selection <- c('all')
  sapply(selection, function(s){
    
    dir.create(paste0('demo_seurat3/',x2,'/',s), showWarnings = FALSE)
    
    # read data counts and cellinfo    
    counts <- read.table(paste0(x,'/counts.txt'), head=T, sep='\t', fill = T)    
    cellinfo <- read.table(paste0(x,'/cellinfo.txt'), head=T, sep='\t', fill = T)
    rownames(cellinfo) <- factor(colnames(counts))
    if (anyNA(counts)) print ("NAs")
    
    pbmc <- CreateSeuratObject(counts = counts, project = '', min.cells = 3, min.features = 200, meta.data = cellinfo)
    pbmc.list <- SplitObject(pbmc, split.by = "Batch")
    
    # Run Seurat V3 integration
    t1 = Sys.time()
    for (i in names(pbmc.list)) {
      pbmc.list[[i]] <- SCTransform(pbmc.list[[i]], verbose = FALSE)
    }
    pbmc.features <- SelectIntegrationFeatures(object.list = pbmc.list, nfeatures = nrow(counts))
    pbmc.list <- PrepSCTIntegration(object.list = pbmc.list, anchor.features = pbmc.features)
    pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, normalization.method = "SCT", anchor.features = pbmc.features, dims = 1:50)
    immune.combined <- IntegrateData(anchorset = pbmc.anchors, normalization.method = "SCT", dims = 1:50)
    t2 = Sys.time()
    print(t2-t1)
    
    # save the output
    save(immune.combined, file = paste0('demo_seurat3/',x2,'/',s,"/output.rda"))
    seuratv3_integrated <- immune.combined@assays$integrated@data
    write.table(seuratv3_integrated, file = paste0('demo_seurat3/',x2,'/',s,"/output.txt"), row.names = T, col.names = T, sep="\t")
  })
})


