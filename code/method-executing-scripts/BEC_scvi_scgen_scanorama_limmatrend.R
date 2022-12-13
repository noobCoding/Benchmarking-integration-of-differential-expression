library(edgeR)
library(limma)
library(Glimma)

# METHODS
# vect_method <- c('scvi', 'scgen', 'scanorama')
vect_method <- c('scgen')

# SIMULATIONS
vect_simu <- gsub('data/','',list.dirs('data', recursive=FALSE))

for (simu in vect_simu){
  for (method in vect_method){
      # import data, sample
      counts<- read.csv(file = paste0('data/',simu, '/', method , '_corrected_data.csv'),header=T,row.names=1,check.names = F)
      cellinfo <- read.table(paste0("data/",simu,"/cellinfo.txt"),sep="\t",header=T,row.names=1)
      rownames(cellinfo) <- colnames(counts)
      count_df <- counts
      
      #### limmatrend
      design <- model.matrix(~Group, data=cellinfo)
      
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
      
      write.table(result.table, file=paste0('data/',simu, '/', method ,"_limmatrend_result.txt"), sep = "\t")
  }
}