# BiocManager::install("apeglm")
# BiocManager::install("DESeq2")
# BiocManager::install("zinbwave")
# BiocManager::install("scRNAseq")
# BiocManager::install("matrixStats")
# BiocManager::install("magrittr")
# BiocManager::install("biomaRt")

library(zinbwave)
library(scRNAseq)
library(matrixStats)
library(magrittr)
library(ggplot2)
library(biomaRt)
library(BiocParallel)

# Register BiocParallel Serial Execution
BiocParallel::register(BiocParallel::SerialParam())
library(DESeq2)

base_names <- list.dirs('data')
base_names <-base_names[grepl('epithelial', base_names)]
# base_name <- base_names[3]
for(base_name in base_names){
  counts<- read.table(file = paste0(base_name,"/counts.txt"), sep = "\t", fill=T)
  cellinfo<-read.table(file = paste0(base_name,"/cellinfo.txt"), sep = "\t", fill=T)
  geneinfo<-read.table(file = paste0(base_name,"/geneinfo.txt"), sep = "\t", fill=T)
  count_df <- counts
  count_df[is.na(count_df)] = 0.
  
  sce <- SingleCellExperiment(assay=list(counts=as.matrix(count_df)), colData=cellinfo)
  rowData(sce)$feature_symbol <- rownames(sce)
  sce_zinb <- zinbwave(sce, K=2, epsilon=1000)
  
  sce_zinb$Group <- factor(sce_zinb$Group)
  sce_zinb$Batch <- factor(sce_zinb$Batch)
  sce_zinb@assays@data$counts <-  sce_zinb@assays@data$counts + 1

  dds <- DESeqDataSet(sce_zinb, design = ~Group)
  dds <- DESeq2::DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6)
  res <- lfcShrink(dds, coef="Group_tLung_vs_nLung", type="apeglm", lfcThreshold=0)
  
  result.table <- data.frame('pvalue' = res$pvalue, 'adjpvalue' = res$padj, 'logFC' = res$log2FoldChange)
  rownames(result.table) <- rownames(dds)
  write.table(result.table, file=paste0(base_name,"/all_zinbwave_deseq2_table.txt"), sep = "\t")
}
