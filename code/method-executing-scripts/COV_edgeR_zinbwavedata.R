#zinbwave output as input to the processed parameter
run_edgeR_zinbwavedata<-function(processed,cellinfo,cov=T,former.meth){
  library(edgeR)
  rownames(cellinfo)=cellinfo$Cell
  cellinfo<-cellinfo[colnames(processed@assays@data$counts),]
  processed$Group <- factor(processed$Group)
  processed$Batch <- factor(processed$Batch)
  cellinfo$Group%<>%factor()
  cellinfo$Batch%<>%factor()
  
  y <- DGEList(assay(processed))
  y <- calcNormFactors(y)
  y$weights <- assay(processed, "weights")
  if(cov){
    design <- model.matrix(~Group+Batch, data=colData(processed))
  }else{
    design <- model.matrix(~Group, data=colData(processed))
  }
  y <- estimateDisp(y, design, robust=TRUE)
  fit <- glmFit(y, design)
  lrt <- glmWeightedF(fit, coef = 2)
  lrt <- topTags(lrt, n = Inf)
  
  res <- data.frame('pvalue' = lrt$table$PValue, 'adjpvalue' = lrt$table$FDR, 'logFC' = lrt$table$logFC)
  rownames(res) <- rownames(lrt)
  save(res, cellinfo, file=paste0('./',former.meth,'edgeR',ifelse(cov,'_Cov',''),'.rda'))
}