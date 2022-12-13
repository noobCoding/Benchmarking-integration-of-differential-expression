#Take DE result(pvalue) of DESeq2/edgeR/limmatrend to res.
#Take count data which was taken as input for DESeq2 and limmatrend to processed. Count data is for calculating sample weight.
run_wFisher<-function(res, processed, former.meth, cellinfo){
  
  weight=get.weight(data.weight=processed,weight.cell='sample',weight.gene = T)
  p=res[['abs']]
  fc=res[['logFC']]
  K <- ncol(p)
  res <- matrix(NA,nrow(p), 4)
  colnames(res)=c('stat','pval','FDR','direction')
  rownames(res)=rownames(p)
  res%<>%as.data.frame()
  
  temp<-get.wFisher.bk(p=p, weight = weight ,fc=fc)
  
  res$stat=temp$stat
  res$pval <- temp$pval
  res$FDR <- temp$FDR
  res$direction<-temp$direction

  save(res, cellinfo, file=paste0('./',former.meth,'wfisher.rda'))
}
get.weight<-function(data.weight, cell.weight='sample', weight.gene=F){
  gene.weight<-sapply(data.weight,FUN = function(x){ (rowSums(x$count!=0)) })
  weight.temp<-sapply(data.weight,FUN = function(x){ rowSums(!is.na(x$count)) })
  if(weight.gene==T){
    res_weight=gene.weight
  }else{
    res_weight=weight.temp
  }
  if(cell.weight=='sqrt'){
    res_weight<-sqrt(res_weight)/ifelse(rowSums(sqrt(res_weight))==0,1,rowSums(sqrt(res_weight)))
  }else if(cell.weight=='sample'){
    res_weight<-res_weight/ifelse(rowSums(res_weight)==0,1,rowSums(res_weight))
  }
  return(res_weight)
}