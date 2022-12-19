#zinbwave output as input to the processed parameter
run_DESeq2_zinbwavedata<-function(processed,cellinfo,cov=T,former.meth='zinbwave'){
  library(DESeq2)
  library(zinbwave)
  rownames(cellinfo)=cellinfo$Cell
  cellinfo<-cellinfo[colnames(processed@assays@data$counts),]
  processed$Group <- factor(processed$Group)
  processed$Batch <- factor(processed$Batch)
  cellinfo$Group%<>%factor()
  cellinfo$Batch%<>%factor()
  
  
  processed@assays@data$counts<-processed@assays@data$counts+1
  if(cov){
    dds <- DESeq2::DESeqDataSet(processed, design = ~Group+Batch)
  }else{
    dds <- DESeq2::DESeqDataSet(processed, design = ~Group)
  }
  dds <- DESeq2::DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6)
  res_table <- lfcShrink(dds, type="apeglm", coef=2, lfcThreshold=0)
  res <- data.frame('pvalue' = res_table$pvalue, 'adjpvalue' = res_table$padj, 'logFC' = res_table$log2FoldChange)
  rownames(res) <- rownames(dds)
  
  res_name<-paste0(ifelse(former.meth=='','',paste0(former.meth,'+')),'deseq2',ifelse(cov,'_Cov',''))
  save(res, cellinfo, file=paste0('./',res_name,'.rda'))
  return(res_name)
}
