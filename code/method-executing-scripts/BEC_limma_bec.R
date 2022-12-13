#raw count as input count=rawcount
run_limma_bec<-function(count,cellinfo){
  library(limma)
  library(Seurat)
  rownames(cellinfo)=cellinfo$Cell
  cellinfo<-cellinfo[colnames(count),]
  cellinfo$Group%<>%factor()
  cellinfo$Batch%<>%factor()
  
  ct <- CreateSeuratObject(counts = count, meta.data = cellinfo, project = "Limma", min.cells = 0)
  ct <- NormalizeData(object = ct, normalization.method = "LogNormalize", scale.factor = 1e4)
  ct <- ScaleData(object = ct)
  norm.count <- ct@assays$RNA@data
  # run limma
  res <- limma::removeBatchEffect(as.matrix(lm_df), factor(cellinfo$Batch))
  processed<-res
  save(res, processed, cellinfo, file='./limma_bec.rda')
}