#raw count as input count=rawcount
run_combat<-function(count,cellinfo,former.meth=''){
  library(magrittr)
  library(dplyr)
  library(Seurat)
  library(sva)
  rownames(cellinfo)=cellinfo$Cell
  cellinfo<-cellinfo[colnames(count),]
  cellinfo$Group%<>%factor()
  cellinfo$Batch%<>%factor()
  
  col_sums = apply(count,2, sum)
  med_trans = median(col_sums)
  norm_counts = med_trans* scale(count, center=FALSE, scale=col_sums)
  norm_counts = log1p(norm_counts)
  count_df <- norm_counts
  
  # Run COMBAT
  combat_output = ComBat(dat=as.matrix(count_df), 
                         batch=cellinfo$Batch,
                         mod=NULL, par.prior=TRUE, prior.plots=FALSE, mean.only=FALSE)
  res<-combat_output
  processed<-combat_output
  
  res_name<-paste0(ifelse(former.meth=='','',paste0(former.meth,'+')),'combat')
  save(res,processed,cellinfo, file=paste0('./',res_name,'.rda'))
  return(res_name)
}