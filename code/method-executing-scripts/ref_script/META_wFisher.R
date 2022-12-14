#Take DE result(pvalue) of DESeq2/edgeR/limmatrend to res.
#Take count data which was taken as input for DESeq2 and limmatrend to processed. Count data is for calculating sample weight.
run_wFisher<-function(res, processed, cellinfo,former.meth=''){
  
  weight=get.weight(data.weight=processed,cell.weight='sample',weight.gene = T)
  p=res[['abs']]
  fc=res[['logFC']]
  K <- ncol(p)
  res <- matrix(NA,nrow(p), 4)
  colnames(res)=c('stat','pval','FDR','direction')
  rownames(res)=rownames(p)
  res%<>%as.data.frame()
  
  temp<-get.wFisher(p=p, weight = weight ,fc=fc)
  
  res$stat=temp$stat
  res$pval <- temp$pval
  res$FDR <- temp$FDR
  res$direction<-temp$direction

  res_name<-paste0(ifelse(former.meth=='','',paste0(former.meth,'_sep','+')),'wfisher')
  save(res, cellinfo, file=paste0('./',res_name,'.rda'))
  return(res_name)
}
get.weight<-function(data.weight, cell.weight='sample', weight.gene=F){
  gene.weight<-sapply(data.weight,FUN = function(x){ (rowSums(x!=0)) })
  weight.temp<-sapply(data.weight,FUN = function(x){ rowSums(!is.na(x)) })
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

get.wFisher<-function(p,weight,fc){
  k <- ncol(p)
  direction<-pval <- stat <- rep(NA, nrow(p))
  rnum <- 1:nrow(p)
  pval[rnum]<-unlist(sapply(X=rnum, FUN=function(x){metapro::wFisher(p[x,],weight=ifelse(all(weight[x,]==0),list(rep(1,ncol(weight))),list(weight[x,]))[[1]],is.onetail = F,eff.sign = c(sign(fc)[x,]))$p}))
  direction[rnum]<-unlist(sapply(X=rnum, FUN=function(x){metapro::wFisher(p[x,],weight=ifelse(all(weight[x,]==0),list(rep(1,ncol(weight))),list(weight[x,]))[[1]],is.onetail = F,eff.sign = c(sign(fc)[x,]))$overall.eff.direction}))
  qval <- p.adjust(pval, method = "BH")
  res <- list(stat = pval, pval = pval, FDR = qval, direction=direction)
  
  names(res$direction)<-names(res$stat) <- names(res$pval) <- names(res$FDR) <- rownames(p)
  return(res)
}
