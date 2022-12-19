#List of Log transformed count as input indexed by each batch input processed=list(batch1=batch1logcount, batch2=batch2logcount, ...)
run_limmatrend_sep<-function(processed,cellinfo,former.meth=''){
  library(limma)
  library(edgeR)
  rownames(cellinfo)=cellinfo$Cell
  cellinfo$Group%<>%factor()
  cellinfo$Batch%<>%factor()
  
  p=matrix(nrow=nrow(processed[[1]]),ncol=length(processed))
  rownames(p)=rownames(processed[[1]])
  colnames(p)=names(processed)
  res<-list()
  res[['abs']]=res[['logFC']]=p
  for(b in unique(cellinfo$Batch)){
    count_df<-processed[[b]]
    cellinfo.temp<-cellinfo[colnames(processed[[b]]),c('Group','Batch')]
    design <- model.matrix(~Group, data=cellinfo.temp)
    lmfit <- lmFit(count_df, design)
    lmfit <- eBayes(lmfit, trend=TRUE, robust = TRUE)
    res_table <- topTable(lmfit, n = Inf, adjust.method = "BH", coef = 2)
    res_table<-res_table[rownames(p),]
    logFC<-res_table$logFC
    pval<-res_table$P.Value
    res[['logFC']][,b]<-logFC
    res[['abs']][,b]<-pval
  }
  
  res_name<-paste0(ifelse(former.meth=='','',paste0(former.meth,'_sep','+')),'limmatrend_sep')
  save(res, processed, cellinfo, file=paste0('./',res_name,'.rda'))
  return(res_name)
}