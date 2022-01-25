rm(list=ls())
# BiocManager::install("limma")
library(limma)
library(Seurat)

lsdir <- list.dirs('data', recursive=FALSE) 
dir.create('demo_limma')

sapply(lsdir,function(x){
  x2 <- gsub('data/','',x)
  dir.create(paste0('demo_limma/',x2), showWarnings = T)
  selection <- c('all')
  
  sapply(selection, function(s){
    dir.create(paste0('demo_limma/',x2,'/',s), showWarnings = FALSE)
    # read data counts and cellinfo   
    counts <- read.table(paste0(x,'/counts.txt'), head=T, sep='\t', fill = T)    
    cellinfo <- read.table(paste0(x,'/cellinfo.txt'), head=T, sep='\t')
    cellinfo <- subset(cellinfo,select=c('Batch','Group'))
    names(cellinfo) <- c('batch','cell_type')
    rownames(cellinfo) <- factor(colnames(counts))
    
    # input limma : normalized matrix
    pbmc <- CreateSeuratObject(counts = counts, meta.data = cellinfo, project = "Limma", min.cells = 0)
    pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 1e4)
    pbmc <- ScaleData(object = pbmc)
    output_df <- pbmc@assays$RNA@data
    
    # run limma
    lm_df <- limma::removeBatchEffect(as.matrix(output_df), factor(cellinfo$batch))

    # save the output
    limma_srt <- CreateSeuratObject(counts = lm_df, meta.data = cellinfo, project = "Limma_normalized")
    save(limma_srt, file = paste0('demo_limma/',x2,'/',s,"/output.rda"))
    write.table(lm_df, file = paste0('demo_limma/',x2,'/',s,"/output.txt"), row.names = T, col.names = T, sep="\t")
  })
})
