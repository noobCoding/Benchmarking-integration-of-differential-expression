#raw count as input count=rawcount
run_limma_bec<-function(count,cellinfo,former.meth=''){
  library(limma)
  library(Seurat)
  rownames(cellinfo)=cellinfo$Cell
  cellinfo<-cellinfo[colnames(count),]
  cellinfo$Group%<>%factor()
  cellinfo$Batch%<>%factor()
  
  ct <- CreateSeuratObject(counts = count, meta.data = cellinfo, project = "Limma", min.cells = 0)
  ct@assays$RNA@data=LogNormalize(ct@assays$RNA@counts)
  # ct <- NormalizeData(object = ct, normalization.method = "LogNormalize", scale.factor = 1e4)
  ct <- ScaleData(object = ct)
  norm.count <- ct@assays$RNA@data
  # run limma
  res <- limma::removeBatchEffect(as.matrix(norm.count), factor(cellinfo$Batch))
  processed<-res
  
  
  
  res_name<-paste0(ifelse(former.meth=='','',paste0(former.meth,'+')),'limma_bec')
  save(res, processed, cellinfo, file=paste0('./',res_name,'.rda'))
  return(res_name)
}