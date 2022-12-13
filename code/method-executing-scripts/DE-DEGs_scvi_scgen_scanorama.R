rm(list=ls())
source('GBC-Seurat_DEG_analysis_auc.R')
#BiocManager::install('limma')
######## S3 batch12 after normalization 

# METHODS
# vect_method <- c('scvi', 'scgen', 'scanorama')
vect_method <- c('scgen')
vect_HVG <- c('all')
# SIMULATIONS
dir.create('Seurat_DEGs_auc')

vect_simu <- gsub('data/','',list.dirs('data', recursive=FALSE))

for (simu in vect_simu){
  for (method in vect_method){
    for (HVG in vect_HVG){
      base_name <- paste0('Seurat_DEGs_auc/',simu,'/')
      dir.create(base_name, showWarnings = FALSE)
      base_name <- paste0(base_name,'S3_batch12/')
      dir.create(base_name, showWarnings = FALSE)
      base_name <- paste0(base_name,'after_',method,'_',HVG,'/')
      dir.create(base_name, showWarnings = FALSE)
      
      # import data, sample
      output_batch12 <- read.csv(file = paste0('data/',simu, '/', method , '_corrected_data.csv'),header=T,row.names=1,check.names = F)
      
      sample_batch12 <- read.table(paste0("data/",simu,"/cellinfo.txt"),sep="\t",header=T,row.names=1)
      #dim(sample_batch12)
      #colnames(output_batch12) = rownames(sample_batch12)
      
      #### DEGs Group1 vs Group2 (batch 1 + batch 2)
      seurat_analysis_deg(TPM=output_batch12,
                          sample=sample_batch12,
                          group_col="Group",
                          base_name=paste0(base_name,'degs_batch12'),  
                          group1='Group1',
                          group2='Group2',
                          test.use="wilcox",
                          logfc.threshold = 0)
    }
  }
}