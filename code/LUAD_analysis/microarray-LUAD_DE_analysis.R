{
  library(GEOquery)
  library(stringr)
  library(dplyr)
  library(magrittr)
  library(tidyverse)
  library(survival)
  library(RegParallel)
  library(biomaRt)
  base_dir='~/'
  source(paste0(base_dir,'script_bk/Bulk sample_analysis.sources.R'))
  surv_dir=paste0(base_dir,'data/for surv/')
  setwd(base_dir)
  scaled_z_weight<-c()
  Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)
  geo_acc_tot=c('GSE43458', "GSE31210")
}

for(i in 1:length(geo_acc_tot)){
  geo_acc=geo_acc_tot[i]
  print(geo_acc)
  gse <- getGEO(geo_acc,GSEMatrix=TRUE)
  gse<-gse[[1]]
  fset <- fData(gse)
  eset <- exprs(gse)
  if(is.null(fset[['Gene Symbol']])){
    if(!is.null(fset[['Symbol']])){
      fset[['Gene Symbol']]=fset$Symbol
    }
  }
  if(!is.null(fset[['Gene Symbol']])){
    eset<-eset[!is.na(fset[['Gene Symbol']]),]
  }
  if(any(grep('^AFFX', rownames(eset)))){
    message('affy out')
    eset <- eset[-grep('^AFFX', rownames(eset)),]
  }
  eset<-do.mapping.microarray.microarray(geo_acc, eset,fset)
  pset<-pData(gse)

  if(geo_acc == c("GSE43458")){
    idx <- match(c('Gender:ch1', 'Age at Diagnosis:ch1','histology:ch1', 'smoking status:ch1'),colnames(pset))
    metadata <- data.frame(pset[,idx[!is.na(idx)]],
                           row.names = rownames(pset))
    colnames(metadata) <- c('gender', 'age', 'histology',
                            'smoking')[!is.na(idx)]
    metadata <- metadata[!apply(metadata, 1, function(x) any( is.na(x) )),]
    
    metadata$histology[metadata$histology=="Lung adenocarcinoma"]='tumor'
    metadata$histology[metadata$histology=="Normal lung tissue"]='normal'
  }else if(geo_acc %in% c("GSE31210")){
    
    idx <- match(c('gender:ch1', 'age (years):ch1','tissue:ch1', 'smoking status:ch1'),colnames(pset))
    metadata <- data.frame(pset[,idx[!is.na(idx)]],
                           row.names = rownames(pset))
    
    colnames(metadata) <- c('gender', 'age', 'histology',
                            'smoking')[!is.na(idx)]
    # metadata$stage[is.na(metadata$stage)]='not reported'
    metadata <- metadata[!apply(metadata, 1, function(x) any( is.na(x) )),]
    metadata$histology[metadata$histology=="primary lung tumor"]='tumor'
    metadata$histology[metadata$histology=="normal lung"]='normal'
    
    eset=log(eset,2)
  }else if(geo_acc %in% c("GSE30219")){
    
    idx <- match(c('gender:ch1', 'age at surgery:ch1','histology:ch1', 'smoking status:ch1'),colnames(pset))
    metadata <- data.frame(pset[,idx[!is.na(idx)]],
                           row.names = rownames(pset))
    
    colnames(metadata) <- c('gender', 'age', 'histology',
                            'smoking')[!is.na(idx)]
    metadata <- metadata[!apply(metadata, 1, function(x) any( is.na(x) )),]
    metadata$histology[metadata$histology=="ADC"]='tumor'
    metadata$histology[metadata$histology=="NTL"]='normal'
  }
  
  metadata%<>%dplyr::filter(histology %in% c('tumor','normal'))
  eset<-eset[,which(colnames(eset)%in%rownames(metadata))]
  metadata<-metadata[match(colnames(eset),rownames(metadata)),]
  
  #All data are log-transformed
  result.table<-run_limma_cov.microarray(dat=eset,metadata=metadata,min.gene.filter=5, log.transformed=T)
  write.table(result.table,file=paste0(base_dir,'data/analyzed/LUAD-',geo_acc,'.microarray','.limma.cov.txt'))
}


