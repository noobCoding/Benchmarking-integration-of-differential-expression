# BiocManager::install('edgeR')
# BiocManager::install("apeglm")
# BiocManager::install("DESeq2")
# BiocManager::install("zinbwave")
# BiocManager::install("scRNAseq")
# BiocManager::install("matrixStats")
# BiocManager::install("magrittr")
# BiocManager::install("biomaRt")
library(edgeR)
library(zinbwave)
# library(scRNAseq)
# library(matrixStats)
library(magrittr)
# library(ggplot2)
library(biomaRt)

base_names <- list.dirs('data')
base_names <-base_names[grepl('epithelial', base_names)]

for(base_name in base_names){
 
  counts<- read.table(file = paste0(base_name,"/counts.txt"), sep = "\t", fill=T)
  cellinfo<-read.table(file = paste0(base_name,"/cellinfo.txt"), sep = "\t", fill=T)
  geneinfo<-read.table(file = paste0(base_name,"/geneinfo.txt"), sep = "\t", fill=T)
  count_df<-counts
  count_df[is.na(count_df)] = 0.
  
  sce <- SingleCellExperiment(assay=list(counts=as.matrix(count_df)), colData=cellinfo)
  rowData(sce)$feature_symbol <- rownames(sce)
  sce_zinb <- zinbwave(sce, K=2, epsilon=1000, observationalWeights = TRUE)
  
  y <- DGEList(assay(sce_zinb))
  y <- calcNormFactors(y)
  
  y$weights <- assay(sce_zinb, "weights")
  
  sce_zinb$Group <- factor(sce_zinb$Group)
  design <- model.matrix(~Group, data=colData(sce_zinb))
  
  y <- estimateDisp(y, design, robust=TRUE)
  fit <- glmFit(y, design)
  lrt <- glmWeightedF(fit, coef = 2)
  
  # fit <- glmQLFit(y, design, robust=TRUE, prior.df = 0)
  # qlf <- glmQLFTest(fit)
  
  # FDR<-p.adjust(qlf$table$PValue,method = "BH")
  # qlf$table$FDR <- FDR
  
  lrt <- topTags(lrt, n = Inf)
  result.table <- data.frame('pvalue' = lrt$table$PValue, 'adjpvalue' = lrt$table$FDR, 'logFC' = lrt$table$logFC)
  rownames(result.table) <- rownames(lrt)
  
  write.table(result.table, file=paste0(base_name,"/all_zinbwave_edger_table.txt"), sep = "\t")
  
}
