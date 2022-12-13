#Take DE result(effectsize) of DESeq2/ES to res.
run_REM<-function(res, cellinfo,former.meth){
  library(MetaDE)
  if(is.null(res$ES)){
    ES=res$logFC
  }else if(is.null(res$logFC)){
    ES=res$ES
  }else{
    stop('Effect size is necessary')
  }
  V=res$Var
  rm(res)
  res.temp<-MetaDE:::get.REM2(em=ES,vm=V)
  res<-data.frame(logFC=mu.hat,pval=res.temp$pval,FDR=res.temp$FDR,row.names=rownames(ES))
  save(res, cellinfo, file=paste0('./',former.meth,'REM.rda'))
}
