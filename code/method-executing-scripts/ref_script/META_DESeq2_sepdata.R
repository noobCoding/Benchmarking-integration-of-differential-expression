#List of raw count as input indexed by each batch input processed=list(batch1=batch1rawcount, batch2=batch2rawcount, ...)
run_DESeq2_sep<-function(processed,cellinfo,former.meth=''){
  library(DESeq2)
  rownames(cellinfo)=cellinfo$Cell
  cellinfo$Group%<>%factor()
  cellinfo$Batch%<>%factor()
  
  p=matrix(nrow=nrow(processed[[1]]),ncol=length(processed))
  rownames(p)=rownames(processed[[1]])
  colnames(p)=names(processed)
  res<-list()
  res[['abs']]=res[['logFC']]=res[['Var']]=p
  for(b in unique(cellinfo$Batch)){
    count_df<-processed[[b]]
    count_df <- round(count_df, 0) + 1 # pseudo count prevents every gene includes at least one 0
    group=factor(cellinfo[colnames(processed[[b]]),'Group'])
    dds <- DESeqDataSetFromMatrix(countData = count_df, colData = data.frame(group=factor(group)), design = ~group)
    dds <- DESeq2::DESeq(dds)
    res_table<-lfcShrink(dds,coef=2, type="apeglm", lfcThreshold=0)
    pval=res_table$pvalue
    logFC=res_table$log2FoldChange
    Var=(res_table$lfcSE)^2
    res[['logFC']][,b]<-logFC
    res[['abs']][,b]<-pval
    res[['Var']][,b]=Var
  }
  
  res_name<-paste0(ifelse(former.meth=='','',paste0(former.meth,'_sep','+')),'deseq2_sep')
  save(res, processed, cellinfo, file=paste0('./',res_name,'.rda'))
  return(res_name)
}