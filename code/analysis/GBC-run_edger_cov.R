# BiocManager::install('edgeR')
library(edgeR)
library(Seurat)

base_names <- list.dirs('data')
base_names <-base_names[grepl('epithelial', base_names)]

for(base_name in base_names){
  counts<- read.table(file = paste0(base_name,"/counts.txt"), sep = "\t", fill=T)
  cellinfo<-read.table(file = paste0(base_name,"/cellinfo.txt"), sep = "\t", fill=T)
  geneinfo<-read.table(file = paste0(base_name,"/geneinfo.txt"), sep = "\t", fill=T)
  count_df<-counts

  # Normalization (Seurat method)  
  y <- DGEList(counts=count_df, group=cellinfo$Group)
  y <- calcNormFactors(y)
  
  cellGroup <- factor(cellinfo$Group)
  cellBatch <- factor(cellinfo$Batch)
  design <- model.matrix(~cellBatch+cellGroup)
  rownames(design) <- colnames(y)
  y <- estimateDisp(y, design, robust=TRUE)
  y$common.dispersion
  
  design <- model.matrix(~cellBatch+cellGroup)
  rownames(design) <- colnames(y)
  
  fit <- glmQLFit(y, design, robust=TRUE, prior.df = 0)
  qlf <- glmQLFTest(fit)
  
  FDR<-p.adjust(qlf$table$PValue,method = "BH")
  qlf$table$FDR <- FDR
  
  result.table <- data.frame('pvalue' = qlf$table$PValue, 'adjpvalue' = qlf$table$FDR, 'logFC' = qlf$table$logFC)
  rownames(result.table) <- rownames(qlf) 
  
  write.table(result.table, file=paste0(base_name,"/all_edger_cov_table.txt"), sep = "\t") 
}
