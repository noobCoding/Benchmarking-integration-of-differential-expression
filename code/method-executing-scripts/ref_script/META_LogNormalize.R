#raw count as input processed=rawcount
run_LogNormalize<-function(count,cellinfo, separate=T,former.meth=''){
  library(Seurat)
  rownames(cellinfo)=cellinfo$Cell
  cellinfo<-cellinfo[colnames(count),]
  cellinfo$Group%<>%factor()
  cellinfo$Batch%<>%factor()
  
  if(separate){
    processed<-list()
    for(b in unique(cellinfo$Batch)){
      processed[[b]]=count[,cellinfo$Cell[which(cellinfo$Batch==b)]]%>%as.matrix()%>%LogNormalize()%>%as.matrix()
    }
    res=processed
  }else{
    res=processed=Seurat::LogNormalize(as.matrix(count))
  }
  
  res_name<-paste0(ifelse(former.meth=='','',paste0(former.meth,'+')),'LogNormalize',ifelse(separate,'_sep',''))
  save(res,processed, cellinfo, file=paste0('./',res_name,'.rda'))
  return(res_name)
}