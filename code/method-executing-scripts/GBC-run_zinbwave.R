# BiocManager::install("apeglm")
# BiocManager::install("DESeq2")
# BiocManager::install("zinbwave")
# BiocManager::install("scRNAseq")
# BiocManager::install("matrixStats")
# BiocManager::install("magrittr")
# BiocManager::install("dplyr")

library(zinbwave)
library(scRNAseq)
library(matrixStats)
library(magrittr)
library(ggplot2)
library(scater)
library(biomaRt)
library(tictoc)
library(cowplot)
library(lubridate)
library(dplyr)

# Register BiocParallel Serial Execution
# BiocParallel::register(BiocParallel::SerialParam())
# library(DESeq2)

base_names <- list.dirs('data')
base_names <-base_names[grepl('epithelial', base_names)]

dir.create('demo_zinbwave', showWarnings = F)

for(base_name in base_names){
  
  x2 <- gsub('data/','', base_name)
  dir.create(paste0('demo_zinbwave/',x2), showWarnings = FALSE)
  dir.create(paste0('demo_zinbwave/',x2,'/all'), showWarnings = FALSE)
  
  counts<- read.table(file = paste0(base_name,"/counts.txt"), sep = "\t", fill=T)
  cellinfo<-read.table(file = paste0(base_name,"/cellinfo.txt"), sep = "\t", fill=T)
  geneinfo<-read.table(file = paste0(base_name,"/geneinfo.txt"), sep = "\t", fill=T)
  count_df<-counts
  count_df[is.na(count_df)] = 0.  
  
  sce <- SingleCellExperiment(assay=list(counts=as.matrix(count_df)), colData=cellinfo)
  rowData(sce)$feature_symbol <- rownames(sce)
  
  sce_zinb_norm <- zinbwave(sce, K=2, X='~Batch', epsilon=1000, normalizedValues=TRUE, residuals = TRUE)
  W <- reducedDim(sce_zinb_norm)
  
  write.table(W, file=paste0('demo_zinbwave/', x2, '/all/', "W_zinbwave_corrected.txt"))
  saveRDS(sce_zinb_norm, paste0('demo_zinbwave/', x2, '/all/output.rds'))
  write.table(sce_zinb_norm@assays@data$normalizedValues, file=paste0('demo_zinbwave/', x2, '/all/output.txt'), quote = F, row.names = T, col.names = T, sep="\t") 
  
}
