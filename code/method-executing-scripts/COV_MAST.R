#raw count as input processed=rawcount
run_MAST<-function(processed,cellinfo, cov=T){
  library(MAST)
  library(Seurat)
  rownames(cellinfo)=cellinfo$Cell
  cellinfo<-cellinfo[colnames(processed),]
  cellinfo$Group%<>%factor()
  cellinfo$Batch%<>%factor()
  
  cellinfo.cov<-cellinfo[,c('Group','Batch')]
  dt<-CreateSeuratObject(counts=processed)
  dt@assays$RNA@data<-LogNormalize(processed)
  dt<-Seurat::ScaleData(dt,do.scale=F,do.center=T, scale.max=10)
  Idents(dt)=cellinfo[colnames(processed),'Group']
  markers<-FindMarkers(dt,ident.1 = sort(levels(Idents(dt)))[2],ident.2 = sort(levels(Idents(dt)))[1],min.pct = 0,logfc.threshold = 0,min.cells.feature = 0,min.cells.group = 0)
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
  save(res, cellinfo, file='./MAST.rda')
}