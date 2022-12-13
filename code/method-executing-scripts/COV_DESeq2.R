#raw count or pseudobulk data as input processed=rawcount
run_DESeq2<-function(processed,cellinfo,cov=T,former.meth=''){
  library(DESeq2)
  count_df<-processed
  rownames(cellinfo)=cellinfo$Cell
  cellinfo<-cellinfo[colnames(processed),]
  cellinfo$Group%<>%factor()
  cellinfo$Batch%<>%factor()
  
  cellinfo.cov<-cellinfo[,c('Group','Batch')]
  
  if(cov){
    design<-model.matrix(formula(paste(c("~ Group", setdiff(colnames(cellinfo.cov),c('Group'))), collapse = '+')), data=cellinfo.cov)
    count_df <- round(count_df, 0) + 1
    dds <- DESeqDataSetFromMatrix(countData = count_df, colData = cellinfo.cov, design = design)
  }else{
    count_df <- round(count_matrix, 0) + 1
    dds <- DESeqDataSetFromMatrix(countData = count_df, colData = cellinfo.cov, design = ~Group)
  }
  dds <- DESeq2::DESeq(dds)
  res_table <- lfcShrink(dds, coef =2 , type="apeglm", lfcThreshold=0)
  res <- data.frame('pvalue' = res_table$pvalue, 'adjpvalue' = res_table$padj, 'logFC' = res_table$log2FoldChange)
  rownames(res) <- rownames(dds)
  save(res, cellinfo, file=paste0('./',ifelse(former.meth=='','',paste0(former.meth,'+')),'deseq2',ifelse(cov,'_Cov',''),'.rda'))
}
