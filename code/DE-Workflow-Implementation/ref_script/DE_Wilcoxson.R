#raw count as input processed=rawcount with is.log=F
#Log transformed BEC count as input processed=processed from BEC output with is.log=T. Also take same processed output when combined with zinbwave
run_wilcox<-function(processed,cellinfo, is.log=T,former.meth=''){
  library(Seurat)
  rownames(cellinfo)=cellinfo$Cell
  cellinfo<-cellinfo[colnames(processed),]
  cellinfo$Group%<>%factor()
  cellinfo$Batch%<>%factor()
  if(is.log){
    dt<-CreateSeuratObject(counts=processed)
    Idents(dt)=cellinfo[colnames(processed),'Group']
    markers<-FindMarkers(dt,ident.1 = sort(levels(Idents(dt)))[2],ident.2 = sort(levels(Idents(dt)))[1],min.pct = 0,logfc.threshold = 0,min.cells.feature = 0,min.cells.group = 0)
  }else{
    dt<-CreateSeuratObject(counts=processed)
    dt@assays$RNA@data<-LogNormalize(processed)
    dt<-Seurat::ScaleData(dt,do.scale=F,do.center=T, scale.max=10)
    Idents(dt)=cellinfo[colnames(processed),'Group']
    markers<-FindMarkers(dt,ident.1 = sort(levels(Idents(dt)))[2],ident.2 = sort(levels(Idents(dt)))[1],min.pct = 0,logfc.threshold = 0,min.cells.feature = 0,min.cells.group = 0)
  }
  res<-markers
  
  res_name<-paste0(ifelse(former.meth=='','',paste0(former.meth,'+')),'wilcox')
  save(res, cellinfo, file=paste0('./',res_name,'.rda'))
  return(res_name)
}
