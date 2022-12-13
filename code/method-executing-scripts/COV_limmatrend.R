#raw count as input processed=rawcount
run_limmatrend<-function(processed,cellinfo,cov=T){
  library(limma)
  library(edgeR)
  count_df<-processed
  rownames(cellinfo)=cellinfo$Cell
  cellinfo<-cellinfo[colnames(processed),]
  cellinfo$Group%<>%factor()
  cellinfo$Batch%<>%factor()
  
  cellinfo.cov<-cellinfo[,c('Group','Batch')]
  ## Convert to an edgeR object
  dgeObj <- DGEList(count_df)
  ## Perform TMM normalisation
  dgeObj <- calcNormFactors(dgeObj)
  logCPM <- cpm(dgeObj, log=TRUE, prior.count=3)
  if(cov){
    design<-model.matrix(formula(paste(c("~ Group", setdiff(colnames(cellinfo.cov),c('Group'))), collapse = '+')), data=cellinfo.cov)
  }else{
    design <- model.matrix(~Group, data=cellinfo.cov)
  }
  lmfit <- lmFit(logCPM, design)
  lmfit <- eBayes(lmfit, trend=TRUE, robust = TRUE)
  res_table <- topTable(lmfit, n = Inf, adjust.method = "BH", coef = 2)
  res <- data.frame('pvalue' = res_table$P.Value, 
                             'adjpvalue' = res_table$adj.P.Val, 
                             'logFC' = res_table$logFC,
                             row.names = rownames(res_table))

  save(res, cellinfo, file=paste0('./limmatrend',ifelse(cov,'_Cov',''),'.rda'))
}
