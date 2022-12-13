#raw count as input count=rawcount
run_mnn<-function(count,cellinfo){
  library(scran)
  library(scales)
  require(Rtsne)
  library(Seurat)
  library(batchelor)
  rownames(cellinfo)=cellinfo$Cell
  cellinfo<-cellinfo[colnames(count),]
  cellinfo$Group%<>%factor()
  cellinfo$Batch%<>%factor()
  
  dt <- CreateSeuratObject(counts = count, meta.data = cellinfo)
  dt <- LogNormalize(dt)
  res <- batchelor::mnnCorrect(as.matrix(dt@assays$RNA@data),batch=cellinfo$Batch, k=20, sigma=0.1, cos.norm.in=TRUE, cos.norm.out=TRUE, var.adj=TRUE, correct.all = T, auto.merge = T)
  processed<-res@assays@data$corrected
  save(res, processed, cellinfo, file='./mnn.RData')
}