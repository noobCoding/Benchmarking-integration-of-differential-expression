rm(list=ls())
source('DE-Seurat_DEG_analysis_auc.R')

######## S3 batch12 after normalization 

# METHODS
vect_method <- c('zinbwave')
vect_HVG <- c('all')

# 
vect_simu <- gsub('data/','',list.dirs('data', recursive=FALSE))

sapply(vect_simu,function(simu){
  sapply(vect_method,function(method){
    sapply(vect_HVG,function(HVG){
  
      base_name <- paste0('Seurat_DEGs_auc/',simu,'/')
      dir.create(base_name, showWarnings = FALSE)
      base_name <- paste0(base_name,'S3_batch12/')
      dir.create(base_name, showWarnings = FALSE)
      base_name <- paste0(base_name,'after_',method,'_',HVG,'/')
      dir.create(base_name, showWarnings = FALSE)
      
      # import data, sample
      inside <- list.files(paste0('demo_',method,'/',simu,'/',HVG), recursive=FALSE)
      if(grepl('.txt',inside[grep('output\\.((txt)|(csv))',inside)])){
        output_batch12 <- read.table(file = paste0('demo_',method,'/',simu,'/',HVG,'/output.txt'),sep="",header=T,row.names=1,check.names = F)
      } else {
        output_batch12 <- read.csv(file = paste0('demo_',method,'/',simu,'/',HVG,'/output.csv'),sep=",",header=T,row.names=1,check.names = F)
        output_batch12 <- t(output_batch12)
      }
      
      sample_batch12 <- read.table(paste0("data/",simu,"/cellinfo.txt"),sep="\t",header=T,row.names=1)
      dim(sample_batch12)
      colnames(output_batch12) = rownames(sample_batch12)
      # rownames(sample_batch12) = colnames(output_batch12)
      
      #### DEGs Group1 vs Group2 (batch 1 + batch 2)
      seurat_analysis_deg(TPM=output_batch12,
                          sample=sample_batch12,
                          group_col="Group",
                          base_name=paste0(base_name,'degs_batch12'),  
                          group1='nLung',
                          group2='tLung',
                          test.use="wilcox",
                          logfc.threshold = 0)
      
    })
  })
})
