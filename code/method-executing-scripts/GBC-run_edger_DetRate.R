# BiocManager::install('edgeR')
library(edgeR)

base_names <- list.dirs('data')
base_names <-base_names[grepl('epithelial', base_names)]

for(base_name in base_names){
  counts<- read.table(file = paste0(base_name,"/counts.txt"), sep = "\t", fill=T)
  cellinfo<-read.table(file = paste0(base_name,"/cellinfo.txt"), sep = "\t", fill=T)
  geneinfo<-read.table(file = paste0(base_name,"/geneinfo.txt"), sep = "\t", fill=T)
  count_df<-counts
  
  # Normalization (Seurat method)
  myFilteredData <- count_df
  rv_genes<-which(apply(myFilteredData,1,var)==0)
  rv_genes_names<-rownames(myFilteredData)[rv_genes]
  count_df<-myFilteredData[!(rownames(myFilteredData) %in% rv_genes_names),]
  geneinfo <- geneinfo[!(rownames(geneinfo) %in% rv_genes_names),]

  y <- DGEList(counts=count_df, group=cellinfo$Group)
  y <- calcNormFactors(y)
  cdr <- scale(colMeans(count_df > 0))
  
  cellGroup <- factor(cellinfo$Group)
  cellBatch <- factor(cellinfo$Batch)
  design <- model.matrix(~cdr + cellGroup)
  rownames(design) <- colnames(y)

  y <- estimateDisp(y, design, robust=TRUE)
  # y$common.dispersion
  
  fit <- glmQLFit(y, design, robust=TRUE, prior.df = 0)
  qlf <- glmQLFTest(fit)
  
  FDR<-p.adjust(qlf$table$PValue,method = "BH")
  qlf$table$FDR <- FDR
  
  result.table <- data.frame('pvalue' = qlf$table$PValue, 'adjpvalue' = qlf$table$FDR, 'logFC' = qlf$table$logFC)
  rownames(result.table) <- rownames(qlf)
  
  write.table(result.table, file=paste0(base_name,"/all_edger_detrate_table.txt"), sep = "\t")
}
