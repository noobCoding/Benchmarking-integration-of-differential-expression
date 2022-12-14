#Take DE result(effectsize) of DESeq2/ES to res.
run_REM<-function(res, cellinfo,former.meth=''){
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
  res<-data.frame(mu.hat=res.temp$mu.hat,pval=res.temp$pval,FDR=res.temp$FDR,row.names=rownames(ES))
  
  res_name<-paste0(ifelse(former.meth=='','',paste0(former.meth,'_sep','+')),'REM')
  save(res, cellinfo, file=paste0('./',res_name,'.rda'))
  return(res_name)
}
