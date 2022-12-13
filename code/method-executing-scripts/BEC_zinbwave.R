#raw count as input count=rawcount
run_zinbwave<-function(count, cellinfo,former.meth=''){
  library(scran)
  library(scales)
  library(zinbwave)
  library(scRNAseq)
  library(matrixStats)
  library(scater)
  library(lubridate)
  rownames(cellinfo)=cellinfo$Cell
  cellinfo<-cellinfo[colnames(count),]
  cellinfo$Group%<>%factor()
  cellinfo$Batch%<>%factor()
  
  sce <- SingleCellExperiment(assay=list(counts=as.matrix(count)), colData=cellinfo)
  rowData(sce)$feature_symbol <- rownames(sce)
  
  sce_zinb <- zinbwave(sce, K=2, epsilon=1000, 
                       normalizedValues=TRUE, residuals = TRUE,
                       observationalWeights = TRUE)
  res<-sce_zinb
  processed<-res@assays$normalizedValues
  save(res,processed,file=paste0('./',ifelse(former.meth=='','',paste0(former.meth,'+')),'zinbwave.rda'))
}