#List of raw count as input indexed by each batch input processed=list(batch1=batch1rawcount, batch2=batch2rawcount, ...)
run_edgeR_sep<-function(processed,cellinfo,former.meth){
  library(edgeR)
  rownames(cellinfo)=cellinfo$Cell
  cellinfo$Group%<>%factor()
  cellinfo$Batch%<>%factor()
  
  p=matrix(nrow=nrow(processed[[1]]),ncol=length(processed))
  rownames(p)=rownames(processed)
  colnames(p)=names(processed)
  res<-list()
  res[['abs']]=res[['logFC']]=p
  for(b in unique(cellinfo$Batch)){
    count_df<-processed[[b]]
    group=factor(cellinfo[colnames(processed[[b]]),'Group'])
    y <- DGEList(counts=count_df, group=group)
    y <- calcNormFactors(y)
    design <- model.matrix(~factor(group))
    rownames(design) <- colnames(y)
    y <- estimateDisp(y, design, robust=TRUE)
    fit <- glmQLFit(y, design, robust=TRUE, prior.df = 0)
    qlf <- glmQLFTest(fit)
    logFC<-qlf$table$logFC
    pval<-qlf$table$PValue
    res[['logFC']][,ind.study]<-logFC
    res[['abs']][,ind.study]<-pval
  }
  save(res, processed, cellinfo, file=paste0('./',former.meth,'edgeR_sep.rda'))
}