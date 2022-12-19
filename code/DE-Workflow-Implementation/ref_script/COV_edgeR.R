#raw count or pseudobulk data as input processed=rawcount
run_edgeR<-function(processed,cellinfo,cov=T,Det=F,former.meth=''){
  library(edgeR)
  count_df<-processed
  rownames(cellinfo)=cellinfo$Cell
  cellinfo<-cellinfo[colnames(processed),]
  cellinfo$Group%<>%factor()
  cellinfo$Batch%<>%factor()
  
  cellinfo.cov<-cellinfo[,c('Group','Batch')]
  
  y <- DGEList(counts=count_df, group=cellinfo.cov$Group)
  y <- calcNormFactors(y)
  cellGroup <- factor(cellinfo.cov$Group)
  cellBatch <- factor(cellinfo.cov$Batch)
  cdr <- scale(colMeans(count_df > 0))
  
  if(Det){
    if(cov){
      design<-model.matrix(~cellGroup+cdr+cellBatch)
    }else{
      design <- model.matrix(~cellGroup+cdr)
    }
  }else{
    if(cov){
      design<-model.matrix(~cellGroup+cellBatch)
    }else{
      design <- model.matrix(~cellGroup)
    }
  }
  rownames(design) <- colnames(y)
  y <- estimateDisp(y, design, robust=TRUE)
  fit <- glmQLFit(y, design, robust=TRUE, prior.df = 0)
  qlf <- glmQLFTest(fit, coef=2)
  FDR<-p.adjust(qlf$table$PValue,method = "BH")
  qlf$table$FDR <- FDR
  res <- data.frame('pvalue' = qlf$table$PValue, 'adjpvalue' = qlf$table$FDR, 'logFC' = qlf$table$logFC)
  rownames(res) <- rownames(qlf)
  
  res_name<-paste0(ifelse(former.meth=='','',paste0(former.meth,'+')),'edgeR',ifelse(Det,'_Detrate',''),ifelse(cov,'_Cov',''))
  save(res, cellinfo, file=paste0('./',res_name,'.rda'))
  return(res_name)
}
