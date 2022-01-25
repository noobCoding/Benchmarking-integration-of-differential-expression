
rm(list=ls())
# BiocManager::install("scran")
# BiocManager::install("scales")
# BiocManager::install("Rtsne")
# BiocManager::install(version='devel')
# BiocManager::install("batchelor")
library(scran)
library(scales)
require(Rtsne)
library(Seurat)
library(ggplot2)
library(cowplot)
library(batchelor)
# 

dir.create('demo_MNN')
lsdir <- list.dirs('data', recursive=FALSE) 

sapply(lsdir,function(x){
  
  x2 <- gsub('data/','',x)
  dir.create(paste0('demo_MNN/',x2), showWarnings = FALSE)
  selection <- c('all')
  sapply(selection, function(s){
    
    dir.create(paste0('demo_MNN/',x2,'/',s), showWarnings = FALSE)
    
    # read data counts and cellinfo
    
    counts <- read.table(paste0(x,'/counts.txt'), head=T, sep='\t', fill=T)
    
    cellinfo <- read.table(paste0(x,'/cellinfo.txt'), head=T, sep='\t', fill=T)
    rownames(cellinfo) <- factor(colnames(counts))
    
    pbmc <- CreateSeuratObject(counts = counts, meta.data = cellinfo)
    pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize")
    
    pbmc.list <- SplitObject(pbmc, split.by = "Batch")    
    myData = c()
    for (i in pbmc.list){
      tmp <- i@assays$RNA@data
      myData <- c(myData, tmp)
    }
    
    View(myData)
    
    # Run MNN
    t1 = Sys.time()
    out.mnn.total <- batchelor::mnnCorrect(myData, k=21, sigma=0.1, cos.norm.in=TRUE, cos.norm.out=TRUE, var.adj=TRUE,
                                           correct.all=TRUE, auto.merge=TRUE)
    t2 = Sys.time()
    print(t2-t1)
    
    # save the output
    save(out.mnn.total, file = paste0('demo_MNN/',x2,'/',s,"/output.rda"))
    corre.mnn <- out.mnn.total@assays@data$corrected # @assays[['corrected']]
    write.table(corre.mnn, file = paste0('demo_MNN/',x2,'/',s,"/output.txt"), row.names = T, col.names = T, sep="\t")
   
  })
  
})


