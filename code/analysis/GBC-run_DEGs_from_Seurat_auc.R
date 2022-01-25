rm(list=ls())
library('Seurat') 
source('GBC-Seurat_DEG_analysis_auc.R')

######## S3 batch12 after normalization 

# METHODS
vect_method <- c('Combat', 'limma', 'MNN', 'seurat3', 'scmerge')
vect_HVG <- c('all')

dir.create('Seurat_DEGs_auc')

# SIMULATIONS
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
        output_batch12 <- read.table(file = paste0('demo_',method,'/',simu,'/',HVG,'/output.txt'),sep="\t",header=T,row.names=1,check.names = F, fill=T)
      } else {
        output_batch12 <- read.table(file = paste0('demo_',method,'/',simu,'/',HVG,'/output.csv'),sep=",",header=T,row.names=1,check.names = F, fill=T)
        output_batch12 <- t(output_batch12)
      }

      sample_batch12 <- read.table(paste0("data/",simu,"/cellinfo.txt"),sep="\t",header=T,row.names=1, fill=T)
      rownames(sample_batch12) = colnames(output_batch12)

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


# Raw
vect_simu <- gsub('data/','',list.dirs('data', recursive=FALSE))

sapply(vect_simu,function(simu){
  sapply(vect_HVG,function(HVG){

    # dir.create('Seurat_DEGs_auc')
    base_name <- paste0('Seurat_DEGs_auc/',simu,'/')
    # dir.create(base_name, showWarnings = FALSE)
    base_name <- paste0(base_name,'S3_batch12/')
    # dir.create(base_name, showWarnings = FALSE)
    base_name <- paste0(base_name,'raw_data_',HVG,'/')
    dir.create(base_name, showWarnings = FALSE)

    sample <- read.table(paste0("data/",simu,"/cellinfo.txt"),sep="\t",header=T,row.names=1, fill = T)

    # import data, sample
    if(HVG=='HVG'){
      data <- read.table(file = paste0('data/',simu,'/counts_HVG.txt'),sep="\t",header=T,row.names=1,check.names = F, fill = T)
      data <- t(data)
      rownames(data) <- gsub('.', '-', rownames(data), fixed = TRUE)
    } else {
      data <- read.table(file = paste0('data/',simu,'/counts.txt'),sep="\t",header=T,row.names=1,check.names = F, fill = T)
    }
    rownames(sample) = colnames(data)

    #### DEGs Group1 vs Group2 (batches 1+2)
    seurat_analysis_deg2(TPM=data,
                         sample=sample,
                         group_col="Group",
                         base_name=paste0(base_name,'degs_batch12_'),
                         group1='nLung',
                         group2='tLung',
                         test.use="wilcox",
                         logfc.threshold = 0)
  })
})
