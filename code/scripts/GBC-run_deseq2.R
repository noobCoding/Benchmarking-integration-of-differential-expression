# BiocManager::install("apeglm")
# BiocManager::install("DESeq2")
library(DESeq2)

base_names <- list.dirs('data')
base_names <-base_names[grepl('epithelial', base_names)]

for(base_name in base_names){
  counts<- read.table(file = paste0(base_name,"/counts.txt"), sep = "\t", fill=T)
  cellinfo<-read.table(file = paste0(base_name,"/cellinfo.txt"), sep = "\t", fill=T)
  geneinfo<-read.table(file = paste0(base_name,"/geneinfo.txt"), sep = "\t", fill=T)
  count_df<-counts 
  
  cellinfo$group <- factor(cellinfo$Group)
  cellinfo$batch <- factor(cellinfo$Batch)
  
  # dds <- DESeqDataSetFromMatrix(countData = count_df, colData = cellinfo, design = ~batch+group)
  dds <- DESeqDataSetFromMatrix(countData = count_df, colData = cellinfo, design = ~group)
  dds <- DESeq2::DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6)
  # res <- lfcShrink(dds, coef="group_tLung_vs_nLung", type="apeglm", lfcThreshold=0)
  res <- lfcShrink(dds, coef=2, type="apeglm", lfcThreshold=0)
  
  result.table <- data.frame('pvalue' = res$pvalue, 'adjpvalue' = res$padj, 'logFC' = res$log2FoldChange)
  rownames(result.table) <- rownames(dds)
  
  write.table(result.table, file=paste0(base_name,"/all_deseq2_result_table.txt"), sep = "\t")
} 

