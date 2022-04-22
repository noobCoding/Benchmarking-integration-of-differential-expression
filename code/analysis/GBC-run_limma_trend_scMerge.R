library(edgeR)
library(limma)
library(Glimma)
library(gplots)
# library(org.Mm.eg.db)

base_names <- list.dirs('data', recursive = F)
base_names <-base_names[grepl('epithelial', base_names)]

for(base_name in base_names){
  x2 <- gsub('data/','',base_name)
  counts<- read.table(file = paste0('demo_scmerge/',x2, '/all/output.txt'), sep = "\t", fill=T)
  cellinfo<-read.table(file = paste0(base_name,"/cellinfo.txt"), sep = "\t", fill=T)
  count_df<-counts
  
  # count_df[is.na(count_df)] = 0.
  # count_df <- count_df + abs(min(count_df))
  
  ## Convert to an edgeR object
  # dgeObj <- DGEList(count_df)
  
  ## Perform TMM normalisation
  # dgeObj <- calcNormFactors(dgeObj)

  design <- model.matrix(~Group+Batch, data=cellinfo)
  
  # limma-trend
  # logCPM <- cpm(dgeObj, log = TRUE, prior.count = 3)
  # lmfit <- lmFit(logCPM, design)
  
  lmfit <- lmFit(count_df, design)
  lmfit <- eBayes(lmfit, trend=TRUE, robust = TRUE)
  res <- topTable(lmfit, n = Inf, adjust.method = "BH", coef = 2)
  
  result.table <- data.frame('pvalue' = res$P.Value, 
                             'adjpvalue' = res$adj.P.Val, 
                             'logFC' = res$logFC,
                             row.names = rownames(res))
  
  write.table(result.table, file=paste0(base_name,"/limma_trend_scMerge_result.txt"), sep = "\t")
}
