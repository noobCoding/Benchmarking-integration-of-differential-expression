rm(list=ls())
#### Result summary

vect_simu <- gsub('data/','',list.dirs('data', recursive=FALSE))
vect_method <- c('scvi_limmatrend', 'scgen_limmatrend', 'scanorama_limmatrend', 'scvi', 'scgen', 'scanorama')

for (simu in vect_simu){
  # simu = 'simul1_dropout_37_b1_300_b2_750'
  load(paste0('data/splatter.', simu, '.RData'))
  extra <- result.list
  rm(result.list)
  result.list<-readRDS(paste0("data/sp80.", simu, '_full.RData'))
  # print(names(extra) %in% names(result.list))
  
  for (method in vect_method){
    if (method=='scvi_limmatrend'){
      scVI_limmatrend <-read.table(paste0('data/',simu,'/','scvi_limmatrend_result.txt'), head=T,  sep='\t', fill=T)
      colnames(scVI_limmatrend) <- c('pval', 'adj.pval', 'log2FC')
      result.list$scVI_limmatrend <- scVI_limmatrend

    } else if (method=='scgen_limmatrend'){
      scGen_limmatrend <-read.table(paste0('data/',simu,'/','scgen_limmatrend_result.txt'), head=T,  sep='\t', fill=T)
      colnames(scGen_limmatrend) <- c('pval', 'adj.pval', 'log2FC')
      result.list$scGen_limmatrend <- scGen_limmatrend

    } else if (method=='scanorama_limmatrend'){
      Scanorama_limmatrend <-read.table(paste0('data/',simu,'/','scanorama_limmatrend_result.txt'), head=T,  sep='\t', fill=T)
      colnames(Scanorama_limmatrend) <- c('pval', 'adj.pval', 'log2FC')
      result.list$Scanorama_limmatrend <- Scanorama_limmatrend
    }
    else if (method=='scvi')
    {
      S3 <- read.table(paste0('Seurat_DEGs_auc/',simu,'/','S3_batch12/','after_',method,'_','all','/degs_batch12_seurat_auc_DEG.txt'), head=T, sep='\t', fill=T)
      scVI_Wilcox <- S3[c('avg_log2FC', 'p_val', 'p_val_adj')]
      colnames(scVI_Wilcox) <- c('log2FC', 'pval', 'adj.pval')
      rownames(scVI_Wilcox) <- S3$X
      result.list$scVI_Wilcox <- scVI_Wilcox
    }
    else if (method=='scgen')
    {
      S3 <- read.table(paste0('Seurat_DEGs_auc/',simu,'/','S3_batch12/','after_',method,'_','all','/degs_batch12_seurat_auc_DEG.txt'), head=T, sep='\t', fill=T)
      scGen_Wilcox <- S3[c('avg_log2FC', 'p_val', 'p_val_adj')]
      colnames(scGen_Wilcox) <- c('log2FC', 'pval', 'adj.pval')
      rownames(scGen_Wilcox) <- S3$X
      result.list$scGen_Wilcox <- scGen_Wilcox
    }
    else if (method=='scanorama')
    {
      S3 <- read.table(paste0('Seurat_DEGs_auc/',simu,'/','S3_batch12/','after_',method,'_','all','/degs_batch12_seurat_auc_DEG.txt'), head=T, sep='\t', fill=T)
      Scanorama_Wilcox <- S3[c('avg_log2FC', 'p_val', 'p_val_adj')]
      colnames(Scanorama_Wilcox) <- c('log2FC', 'pval', 'adj.pval')
      rownames(Scanorama_Wilcox) <- S3$X
      result.list$Scanorama_Wilcox <- Scanorama_Wilcox
    }
  }
  
  tmp <- names(result.list)
  tmp[ which(tmp=="Pseudobulk_limma_trend") ] <- 'Pseudobulk_limmatrend'
  tmp[ which(tmp=="RISC_limma_trend") ] <- 'RISC_limmatrend'
  names(result.list) <- tmp
  
  for (name in names(extra)){
    result.list[[name]] <- extra[[name]]
  }
  
  counts <- read.table(paste0('data/',simu,'/counts.txt'), head=T, fill=T)
  geneinfo <- read.table(paste0('data/',simu,'/geneinfo.txt'), head=T, fill=T)
  up_genes <- read.table(paste0('data/',simu,'/true_up_genes.txt'), head=T, fill = T)
  down_genes <- read.table(paste0('data/',simu,'/true_down_genes.txt'), head=T, fill = T)
  de_genes <- read.table(paste0('data/',simu,'/de__genes.txt'), head=T, fill = T)
  
  GeneInfo <- list()
  GeneInfo$All <- geneinfo$Gene
  GeneInfo$Up <- up_genes$x
  GeneInfo$Down <- down_genes$x
  
  result.list$GeneInfo <- GeneInfo
  saveRDS(result.list, paste0("data/sp80.", simu, '_full.RData'))
}
    