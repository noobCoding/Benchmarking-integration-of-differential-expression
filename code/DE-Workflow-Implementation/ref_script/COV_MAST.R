#raw count as input processed=rawcount
run_MAST<-function(processed,cellinfo, cov=T,former.meth=''){
  library(MAST)
  library(Seurat)
  rownames(cellinfo)=cellinfo$Cell
  cellinfo<-cellinfo[colnames(processed),]
  cellinfo$Group%<>%factor()
  cellinfo$Batch%<>%factor()
  
  cellinfo.cov<-cellinfo[,c('Group','Batch')]
  dt<-CreateSeuratObject(counts=processed,meta.data=cellinfo.cov,min.cells=0)
  dt@assays$RNA@data<-LogNormalize(processed)
  Idents(dt)=cellinfo[colnames(processed),'Group']
  if(cov){
    markers <- FindMarkers(object = dt, ident.1 = sort(levels(Idents(dt)))[2],ident.2 = sort(levels(Idents(dt)))[1],
                           test.use='MAST', logfc.threshold = 0,
                           latent.vars=setdiff(colnames(cellinfo.cov),c('Group')),
                           min.cells.feature = 0,
                           min.cells.group = 0,
                           min.pct = 0,
                           only.pos = F)
  }else{
    markers <- FindMarkers(object = dt, ident.1 = sort(levels(Idents(dt)))[2],ident.2 = sort(levels(Idents(dt)))[1],
                           test.use='MAST', logfc.threshold = 0,
                           min.cells.feature = 0,
                           min.cells.group = 0,
                           min.pct = 0,
                           only.pos = F)
  }
  res<-markers
  
  res_name<-paste0(ifelse(former.meth=='','',paste0(former.meth,'+')),'MAST',ifelse(cov,'_Cov',''))
  save(res, cellinfo,cellinfo, file=paste0('./',res_name,'.rda'))
  return(res_name)
}