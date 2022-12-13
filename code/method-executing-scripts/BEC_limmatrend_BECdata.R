#Log transformed BEC count as input processed=processed from BEC output
run_limmatrend_BECdata<-function(processed,cellinfo,former.meth=''){
  library(limma)
  library(edgeR)
  count_df<-processed
  rownames(cellinfo)=cellinfo$Cell
  cellinfo<-cellinfo[colnames(processed),]
  cellinfo$Group%<>%factor()
  cellinfo$Batch%<>%factor()
  
  design <- model.matrix(~Group, data=cellinfo)
  lmfit <- lmFit(count_df, design)
  lmfit <- eBayes(lmfit, trend=TRUE, robust = TRUE)
  res_table <- topTable(lmfit, n = Inf, adjust.method = "BH", coef = 2)
  res <- data.frame('pvalue' = res_table$P.Value, 
                    'adjpvalue' = res_table$adj.P.Val, 
                    'logFC' = res_table$logFC,
                    row.names = rownames(res_table))
  
  save(res, cellinfo, file=paste0('./',ifelse(former.meth=='','',paste0(former.meth,'+')),'limmatrend.rda'))
}