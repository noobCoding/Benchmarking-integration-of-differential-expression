base_dir='set base directory/'
setwd(base_dir)
{
  library(GEOquery)
  library(stringr)
  library(dplyr)
  library(magrittr)
  library(tidyverse)
  library(biomaRt)
  library(limma)
  # Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)
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
  eset<-do.mapping.microarray(geo_acc, eset,fset)
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
  result.table<-run_limma_cov.microarray(dat=eset,metadata=metadata, log.transformed=T)
  write.table(result.table,file=paste0('LUAD-',geo_acc,'.microarray','.limma.cov.txt'))
}

run_limma_cov.microarray<-function(dat, metadata, log.transformed=T){
  count_df<-dat
  
  library(edgeR)
  library(limma)
  
  cellinfo<-metadata
  cellinfo$sample=NULL
  for(i in colnames(cellinfo)){
    if(i=='age'){
      cellinfo[[i]]%<>%as.numeric()
    }else{
      cellinfo[[i]]%<>%factor()
    }
  }
  
  expset<-ExpressionSet(assayData=count_df%>%as.matrix())
  if(all(colnames(cellinfo)=="histology")){
    design<-model.matrix(formula("~ histology"), data=cellinfo)
  }else{
    design<-model.matrix(formula(paste(c("~ histology", setdiff(colnames(cellinfo),c('histology'))), collapse = '+')), data=cellinfo)
  }
  
  fit <- lmFit(expset, design)
  fitbayes <- limma::eBayes(fit)
  pvalues <- fitbayes$p.value[, 2]
  adjpvalues <- p.adjust(pvalues, method = 'BH')
  logFC <- fitbayes$coefficients[, 2]
  score <- 1 - pvalues
  result.table <- data.frame('pvalue' = pvalues, 'adjpvalue' = adjpvalues, 'logFC' = logFC, 'score' = score)
  rownames(result.table) <- rownames(count_df)
  return(result.table)
}
do.mapping.microarray<-function(geo_acc,eset,fset){
  mart <- useMart('ENSEMBL_MART_ENSEMBL',host = 'useast.ensembl.org')
  mart <- useDataset("hsapiens_gene_ensembl", mart)
  if(geo_acc=='GSE32863'){
    platform='illumina_humanwg_6_v3'
    filt=attr='external_gene_name'
    
  }else if(geo_acc=='GSE10072'){
    platform='affy_hg_u133a'
    attr=filt='affy_hg_u133a'
  }else if(geo_acc=='GSE43458'){
    platform='affy_hugene_1_0_st_v1'
    filt=attr='affy_hugene_1_0_st_v1'
  }else if(geo_acc=="GSE31210"){
    platform='affy_hg_u133_plus_2'
    attr=filt='affy_hg_u133_plus_2'
  }else{
    platform='affy_hg_u133_plus_2'
    filt=attr='external_gene_name'
  }
  mrna_attributes <- getBM(mart = mart,
                           # attributes = c(platform,
                           attributes = c(platform,
                                          'ensembl_gene_id',
                                          'gene_biotype',
                                          'external_gene_name'),
                           filter = filt,
                           values = rownames(eset),
                           uniqueRows = TRUE)
  mrna_attributes%<>%dplyr::filter(gene_biotype=='protein_coding')
  eset<-eset[rownames(eset)%in%mrna_attributes[[attr]],]
  rownames(eset)<-mrna_attributes$external_gene_name[match(rownames(eset),mrna_attributes[[attr]])]
  eset<-eset[str_detect(rownames(eset),pattern='[///]', negate = T)%>%which(),]
  {
    eset.1<-eset[!duplicated(rownames(eset)),]
    eset.2<-eset[duplicated(rownames(eset)),]
    probavg<-apply(eset.2, 1, mean)%>%as.vector()
    
    maxprob<-sapply(unique(rownames(eset.2)),FUN=function(x){
      m1<-probavg[x==rownames(eset.2)]
      return(which(x==rownames(eset.2))[which(m1==max(m1))])
    })%>%unlist()%>%as.vector()
    res_dt<-rbind(eset.1,eset.2[maxprob,])
    
    res_dt%<>%as.data.frame()
    res_dt<-res_dt[setdiff(unique(rownames(res_dt))%>%sort(),c("")),]
  }
  return(res_dt)
}