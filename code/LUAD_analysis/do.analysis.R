library(stringr)
library(data.table)
library(feather)
library(dplyr)
library(magrittr)
library(tidyverse)
base_dir='~/'
script_dir=paste0(base_dir,'script_bk/')
source(paste0(script_dir,'de.analysis.R'))

do.zinbwave<-function(ct,base_dir,filter.rate = 0.01, filt_dir, count.base,FP='none',seed=1234,vs='none'){
  library(scran)
  library(scales)
  require(Rtsne)
  library(zinbwave)
  library(scRNAseq)
  library(matrixStats)
  library(magrittr)
  library(scater)
  library(biomaRt)
  library(tictoc)
  library(lubridate)
  library(dplyr)
  dir.create(paste0(base_dir,'data/',filt_dir,'/zinbwave'),showWarnings = F, recursive = T)
  if(FP=='tumor'){
    annot<-fread(paste0(base_dir,'data/','/FP/original','/',count.base,ct,'_FP.',seed,'seed.cell_annotation+patient.txt'))
  }else if(FP=='normal'){
    annot<-fread(paste0(base_dir,'data/','/FP/original','/',count.base,ct,'_FP2.',seed,'seed.cell_annotation+patient.txt'))
  }else if(FP=='none'){
    annot<-fread(paste0(base_dir,'data/','/original','/GSE131907_Lung_Cancer_cell_annotation+patient.txt'))
  }
  count<-read_feather(paste0(base_dir,'data/',filt_dir,'/analysis','/',count.base,ct,ifelse(is.null(filter.rate),'',paste0('.',filter.rate,'mincell_filtered')),'.feather'))
  GeneName<-count$GeneName
  count<-as.matrix(count[,2:ncol(count)])
  rownames(count)<-GeneName
  if(ct=='T-NK cells'){
    annot%<>%dplyr::filter(Cell_type.refined==gsub(pattern='[-]',replacement='/',ct))%>%dplyr::filter(Index %in% colnames(count))
  }else{
    annot%<>%dplyr::filter(Cell_type==ct)%>%dplyr::filter(Index %in% colnames(count))
  }
  if(vs=='subct'){
    annot$Sample_Origin=annot$Cell_subtype
  }else if(vs=='ct'){
    annot$Sample_Origin=annot$Cell_type.refined
  }else if(vs=='none'){
  }
  sample_names<-count%>%colnames()
  ##annot$Patient contains patient number
  ##annot$Sample_Origin is consist of 'nLung' and 'tLung'
  cellinfo<-data.frame(batch=annot$Patient[match(sample_names,annot$Index)],
                       group=annot$Sample_Origin[match(sample_names,annot$Index)],
                       row.names = sample_names)
  count_df<-count
  
  myFilteredData <- count_df
  rv_genes<-which(apply(myFilteredData,1,var)==0)
  rv_genes_names<-rownames(myFilteredData)[rv_genes]
  count_df<-myFilteredData[!(rownames(myFilteredData) %in% rv_genes_names),]
  count_df[is.na(count_df)] = 0.
  
  sce <- SingleCellExperiment(assay=list(counts=as.matrix(count_df)), colData=cellinfo)
  rowData(sce)$feature_symbol <- rownames(sce)
  
  sce_zinb <- zinbwave(sce, K=2, epsilon=1000, 
                       normalizedValues=TRUE, residuals = TRUE,
                       observationalWeights = TRUE)
  W <- reducedDim(sce_zinb)
  write.table(W,paste0(base_dir,'data/',filt_dir,'/zinbwave/',ct,ifelse(is.null(filter.rate),'',paste0('.',filter.rate,'mincell_filtered')),'_',"W_zinbwave_corrected.txt"))
  saveRDS(sce_zinb, paste0(base_dir,'data/',filt_dir,'/zinbwave/',ct,ifelse(is.null(filter.rate),'',paste0('.',filter.rate,'mincell_filtered')),'_','LC_output.rds'))
  write.table(sce_zinb@assays@data$normalizedValues, file=paste0(base_dir,'data/',filt_dir,'/zinbwave/',ct,ifelse(is.null(filter.rate),'',paste0('.',filter.rate,'mincell_filtered')),'_','LC_zinbwave_normalized.txt'), 
              quote = F, row.names = T, col.names = T, sep="\t")
}

do.zinbwave.DEG<-function(ct,base_dir,filter.rate = 0.01, filt_dir, count.base,FP='none',seed=1234,DEG.meth='DESeq2'){
  library(scran)
  library(scales)
  require(Rtsne)
  library(zinbwave)
  library(scRNAseq)
  library(matrixStats)
  library(magrittr)
  library(scater)
  library(biomaRt)
  library(tictoc)
  library(lubridate)
  library(dplyr)
  sce_zinb <- readRDS(paste0(base_dir,'data/',filt_dir,'/zinbwave/',ct,ifelse(is.null(filter.rate),'',paste0('.',filter.rate,'mincell_filtered')),'_','LC_output.rds'))
  sce_zinb$group <- factor(sce_zinb$group)
  sce_zinb$batch <- factor(sce_zinb$batch)
  if(DEG.meth=='DESeq2_pseudo'){
    library(DESeq2)
    sce_zinb@assays@data$counts<-sce_zinb@assays@data$counts+1
    dds <- DESeq2::DESeqDataSet(sce_zinb, design = ~group)
    dds <- DESeq2::DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6)
    res <- lfcShrink(dds, type="apeglm", coef=2, lfcThreshold=0)
    
    result.table <- data.frame('pvalue' = res$pvalue, 'adjpvalue' = res$padj, 'logFC' = res$log2FoldChange)
    rownames(result.table) <- rownames(dds)
    
  }else if(DEG.meth=='DESeq2_pseudo_Cov'){
    library(DESeq2)
    sce_zinb@assays@data$counts<-sce_zinb@assays@data$counts+1
    dds <- DESeq2::DESeqDataSet(sce_zinb, design = ~group+batch)
    dds <- DESeq2::DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6)
    res <- lfcShrink(dds, type="apeglm", coef=2, lfcThreshold=0)
    
    result.table <- data.frame('pvalue' = res$pvalue, 'adjpvalue' = res$padj, 'logFC' = res$log2FoldChange)
    rownames(result.table) <- rownames(dds)
    
  }else if(DEG.meth=='edgeR'){
    library(edgeR)
    
    y <- DGEList(assay(sce_zinb))
    y <- calcNormFactors(y)
    y$weights <- assay(sce_zinb, "weights")
    
    design <- model.matrix(~group, data=colData(sce_zinb))
    
    y <- estimateDisp(y, design, robust=TRUE)
    fit <- glmFit(y, design)
    lrt <- glmWeightedF(fit, coef = 2)
    
    lrt <- topTags(lrt, n = Inf)
    result.table <- data.frame('pvalue' = lrt$table$PValue, 'adjpvalue' = lrt$table$FDR, 'logFC' = lrt$table$logFC)
    rownames(result.table) <- rownames(lrt)
    
  }else if(DEG.meth=='edgeR_Cov'){
    library(edgeR)
    
    y <- DGEList(assay(sce_zinb))
    y <- calcNormFactors(y)
    y$weights <- assay(sce_zinb, "weights")
    
    design <- model.matrix(~group+batch, data=colData(sce_zinb))
    
    y <- estimateDisp(y, design, robust=TRUE)
    fit <- glmFit(y, design)
    lrt <- glmWeightedF(fit, coef=2)
    
    lrt <- topTags(lrt, n = Inf)
    result.table <- data.frame('pvalue' = lrt$table$PValue, 'adjpvalue' = lrt$table$FDR, 'logFC' = lrt$table$logFC)
    rownames(result.table) <- rownames(lrt)
  }else if(DEG.meth=='findmarkers'){
    library(Seurat)
    dt<-CreateSeuratObject(counts=sce_zinb@assays@data$normalizedValues)
    
    Idents(dt)=sce_zinb$group
    
    result.table<-FindMarkers(dt,ident.1 = sort(unique(Idents(dt)))[2],ident.2 = sort(unique(Idents(dt)))[1],min.pct = 0,logfc.threshold = 0)
  }
  write.table(result.table, file=  paste0(base_dir,'data/',filt_dir,'/zinbwave/',ct,ifelse(is.null(filter.rate),'',paste0('.',filter.rate,'mincell_filtered')),'_zinbwave||',DEG.meth,'_third_processed.txt'), sep = "\t")
  third_processed<-result.table
  save(third_processed,file=  paste0(base_dir,'data/',filt_dir,'/third_processed/',ct,ifelse(is.null(filter.rate),'',paste0('.',filter.rate,'mincell_filtered')),'_zinbwave||',DEG.meth,'_third_processed.RData'))
}


do.filter<-function(ct,base_dir, filter.rate=0.01,filt_dir,count.base){
  dir.create(paste0(base_dir,'data/',filt_dir),showWarnings = F)
  dir.create(paste0(base_dir,'data/',filt_dir,'/analysis'),showWarnings = F)
  dir.create(paste0(base_dir,'data/',filt_dir,'/first_processed'),showWarnings = F)
  dir.create(paste0(base_dir,'data/',filt_dir,'/second_processed'),showWarnings = F)
  dir.create(paste0(base_dir,'data/',filt_dir,'/third_processed'),showWarnings = F)
  library(Seurat)
  
  annot<-fread(paste0(base_dir,'data/','/original','/GSE131907_Lung_Cancer_cell_annotation+patient.txt'))
  count<-read_feather(paste0(base_dir,'data/','/analysis','/',count.base,ct,'.feather'))
  
  GeneName<-count$GeneName
  count_matrix<-as.matrix(count[,2:ncol(count)])
  rownames(count_matrix)<-GeneName
  dt<-CreateSeuratObject(counts=count_matrix, min.cells = filter.rate*ncol(count_matrix))
  g1<-dt@assays$RNA@counts%>%rownames()
  write_feather(count[which(count$GeneName %in% g1),],path=paste0(base_dir,'data/',filt_dir,'data/analysis','/',count.base,ct,ifelse(is.null(filter.rate),'',paste0('.',filter.rate,'mincell_filtered')),'.feather'))
}

do.filter.fp<-function(ct,base_dir, filter.rate=0.01,filt_dir,count.base,FP='normal',seed=1234, do.equal=T){
  dir.create(paste0(base_dir,'data/FP/'),showWarnings = F,recursive = T)
  dir.create(paste0(base_dir,'data/FP/',filt_dir),showWarnings = F,recursive = T)
  dir.create(paste0(base_dir,'data/FP/original/'),showWarnings = F,recursive = T)
  dir.create(paste0(base_dir,'data/FP/analysis/'),showWarnings = F,recursive = T)
  dir.create(paste0(base_dir,'data/FP/',filt_dir,'/first_processed'),showWarnings = F)
  dir.create(paste0(base_dir,'data/FP/',filt_dir,'/second_processed'),showWarnings = F)
  dir.create(paste0(base_dir,'data/FP/',filt_dir,'/third_processed'),showWarnings = F)
  
  library(Seurat)
  
  annot<-fread(paste0(base_dir,'data/original','/GSE131907_Lung_Cancer_cell_annotation+patient.txt'))
  annot<-annot[,2:ncol(annot)]
  if(FP=='tumor'){
    split.sample='tLung'
  }else if(FP=='normal'){
    split.sample='nLung'
  }
  count<-read_feather(paste0(base_dir,'data/analysis','/',count.base,ct,'.feather'))
  
  GeneName<-count$GeneName
  count_matrix<-as.matrix(count[,2:ncol(count)])
  rownames(count_matrix)<-GeneName
  
  annot%<>%dplyr::filter(Sample_Origin==split.sample)
  count_matrix<-count_matrix[,colnames(count_matrix)%in%annot$Index]
  dt<-CreateSeuratObject(counts=count_matrix, min.cells = filter.rate*ncol(count_matrix))
  g1<-dt@assays$RNA@counts%>%rownames()
  count<-count[which(count$GeneName %in% g1),]
  annot%<>%dplyr::filter(Index %in% colnames(count))
  
  num_p<-table(annot$Patient)
  set.seed(seed)
  if(do.equal){
    random.fac=rep(1,length(unique(annot$Patient)))
    message('equal split')
  }else{
    #seeds=c(1234:1247)
    
    random.fac=rep(c(4,7/3,1.5,1,2/3,3/7,1/4)[seed%%7],length(unique(annot$Patient)))
  }
  names(random.fac)=unique(annot$Patient)
  for(p in unique(annot$Patient)){
    set.seed(seed)
    annot$Sample_Origin[sample(which(annot$Patient==p),round(num_p[p]/(1+random.fac[p])))]=setdiff(c('nLung','tLung'),split.sample)
  }
  
  if(FP=='tumor'){
    if(ct=="Epithelial cells"){
      stop('tumor Epi patient not supported')
      annot%<>%dplyr::filter(Patient %in% c("P0018", "P0019", "P0020", "P0030", "P0034"))
    }
    write.table(annot,paste0(base_dir,'data/FP/original','/',count.base,ct,'_FP.',seed,'seed.cell_annotation+patient.txt'))
  }else if(FP=='normal'){
    write.table(annot,paste0(base_dir,'data/FP/original','/',count.base,ct,'_FP2.',seed,'seed.cell_annotation+patient.txt'))
  }else if(FP=='none'){
    message('Must set sample groups to perform FP analysis.')
    break()
    write.table(annot,paste0(base_dir,'data/FP/original','/',count.base,ct,'.cell_annotation+patient.txt'))
  }
  count<-count[,c(1,which(colnames(count)%in%annot$Index))]
  write_feather(count,path=paste0(base_dir,'data/',filt_dir,'/FP/analysis','/',count.base,ct,ifelse(is.null(filter.rate),'',paste0('.',filter.rate,'mincell_filtered')),'.feather'))
}



do.first<-function(ct,meth.comb,base_dir, filter.rate=0.01, filt_dir,count.base,FP='none',seed=1234, vs='original'){
  ind.analysis=F
  first.meth=str_split(meth.comb,pattern='[+]')[[1]][1]
  second.meth=str_split(meth.comb,pattern='[+]')[[1]][2]
  message(second.meth)
  if(FP=='tumor'){
    annot<-fread(paste0(base_dir,'data/original','/',count.base,ct,'_FP.',seed,'seed.cell_annotation+patient.txt'))
  }else if(FP=='normal'){
    annot<-fread(paste0(base_dir,'data/original','/',count.base,ct,'_FP2.',seed,'seed.cell_annotation+patient.txt'))
  }else if(FP=='none'){
    annot<-fread(paste0(base_dir,'data/original','/GSE131907_Lung_Cancer_cell_annotation+patient.txt'))
  }
  
  count<-read_feather(paste0(base_dir,'data/',filt_dir,'/analysis/',count.base,ct,ifelse(is.null(filter.rate),'',paste0('.',filter.rate,'mincell_filtered')),'.feather'))
  
  GeneName<-count$GeneName
  count<-as.matrix(count[,2:ncol(count)])
  rownames(count)<-GeneName
  annot<-annot[,2:ncol(annot)]
  if(ct=='T-NK cells'){
    annot%<>%dplyr::filter(Cell_type.refined==gsub(pattern='[-]',replacement='/',ct))%>%dplyr::filter(Index %in% colnames(count))
  }else{
    annot%<>%dplyr::filter(Cell_type==ct)%>%dplyr::filter(Index %in% colnames(count))
  }
  if(str_detect(second.meth,pattern='ind')){
    message('ind true')
    ind.analysis=T
  }
  if(vs=='subct'){
    annot$Sample_Origin=annot$Cell_subtype
  }else if(vs=='ct'){
    annot$Sample_Origin=annot$Cell_type.refined
  }else if(vs=='original'){
  }
  
  if(first.meth%in%c('','raw')){
    first_processed<-get.raw(count, annot, ind.analysis)
  }else if(first.meth=='mnn'){
    first_processed<-get.mnn(count, annot, ind.analysis)
  }else if(first.meth=='mnn_hai'){
    first_processed<-get.mnn_hai(count, annot, ind.analysis)
  }else if(first.meth=='mnn_opt'){
    first_processed<-get.mnn_opt(count, annot, ind.analysis)
  }else if(first.meth=='LogNormalize'){
   first_processed<-get.LogNormalize(count, annot, ind.analysis)
  }else if(first.meth=='scMerge'){
    first_processed<-get.scMerge(count, annot, ind.analysis)
  }else if(first.meth=='combat'){
    first_processed<-get.combat(count, annot, ind.analysis)
  }else if(first.meth=='voom'){
    first_processed<-get.voom(count, annot, ind.analysis)
  }else if(first.meth=='limma_bec'){
    first_processed<-get.limma_bec(count, annot, ind.analysis)
  }else if(first.meth=='Seurat'){
    first_processed<-get.Seurat(count, annot, ind.analysis)
  }
  message(paste0(base_dir,filt_dir,'/first_processed/',ct,ifelse(is.null(filter.rate),'',paste0('.',filter.rate,'mincell_filtered')),'_',first.meth,ifelse(ind.analysis,'.ind',''),'_first_processed.RData'))
  save(first_processed,file=paste0(base_dir,filt_dir,'/first_processed/',ct,ifelse(is.null(filter.rate),'',paste0('.',filter.rate,'mincell_filtered')),'_',first.meth,ifelse(ind.analysis,'.ind',''),'_first_processed.RData'))
}
do.second<-function(ct,meth.comb,base_dir,filter.rate=0.01,filt_dir,count.base,FP='none',seed=1234){
  ind.analysis=F
  first.meth=str_split(meth.comb,pattern='[+]')[[1]][1]
  second.meth=str_split(meth.comb,pattern='[+]')[[1]][2]
  third.meth=str_split(meth.comb,pattern='[+]')[[1]][3]
  if(str_detect(second.meth,pattern='ind')){
    ind.analysis=T
  }
  load(paste0(base_dir,'data/',filt_dir,'/first_processed/',ct,ifelse(is.null(filter.rate),'',paste0('.',filter.rate,'mincell_filtered')),'_',first.meth,ifelse(ind.analysis,'.ind',''),'_first_processed.RData'))
  if(second.meth%in%c('')){
    second_processed<-get.raw.second(first_processed)
  }else if(second.meth=='ind.ES'){
    second_processed<-get.ES.ind(first_processed, using='metaDE')
  }else if(second.meth=='ind.DESeq2.ES'){
    second_processed<-get.ES.ind(first_processed, using='DESeq2')
  }else if(second.meth=='ind.DESeq2.pval'){
    second_processed<-get.pval.ind(first_processed, using='DESeq2')
  }else if(second.meth=='ind.edgeR.pval'){
    second_processed<-get.pval.ind(first_processed, using='edgeR')
  }else if(second.meth=='ind.findmarkers'){
    message('ind.findmarkers')
    second_processed<-get.pval.ind(first_processed, using='findmarkers')
  }else if(second.meth=='ind.limma_trend'){
    second_processed<-get.pval.ind(first_processed, using='limma_trend')
  }else if(second.meth=='ind.modt'){
    second_processed<-get.pval.ind(first_processed, using='modt')
  }
  save(second_processed,file=paste0(base_dir,'data/',filt_dir,'/second_processed/',ct,ifelse(is.null(filter.rate),'',paste0('.',filter.rate,'mincell_filtered')),'_',first.meth,'|',second.meth,'_second_processed.RData'))
}

do.third<-function(ct,meth.comb,base_dir,filter.rate=0.01,filt_dir,count.base,FP='none',seed=1234){
  ind.analysis=F
  first.meth=str_split(meth.comb,pattern='[+]')[[1]][1]
  second.meth=str_split(meth.comb,pattern='[+]')[[1]][2]
  third.meth=str_split(meth.comb,pattern='[+]')[[1]][3]
  if(str_detect(second.meth,pattern='ind')){
    ind.analysis=T
  }
  load(paste0(base_dir,'data/',filt_dir,'/second_processed/',ct,ifelse(is.null(filter.rate),'',paste0('.',filter.rate,'mincell_filtered')),'_',first.meth,'|',second.meth,'_second_processed.RData'))
  if(third.meth == 'findmarkers'){
    if(first.meth=='raw'){
      third_processed<-do.findmarker(second_processed)
    }else{
      third_processed<-do.findmarker.prenorm(second_processed)
    }
  }else if(third.meth == 'limma_trend'){
    if(first.meth=='raw'){
      third_processed<-do.limma_trend.raw(second_processed, meth='sep')
    }else{
      third_processed<-do.limma_trend(second_processed)
    }
  }else if(third.meth == 'limma_trend_Cov'){
    if(first.meth=='raw'){
      third_processed<-do.limma_trend.raw(second_processed, meth='cov')
    }
  }else if(third.meth == 'DESeq2'){
    third_processed<-do.DESeq2(second_processed, meth='sep')
  }else if(third.meth == 'DESeq2_Cov'){
    third_processed<-do.DESeq2(second_processed, meth='cov')
  }else if(third.meth == 'edgeR'){
    third_processed<-do.edgeR(second_processed, meth='sep')
  }else if(third.meth == 'edgeR_Cov'){
    third_processed<-do.edgeR(second_processed, meth='cov')
  }else if(third.meth == 'edgeR_DetRate'){
    third_processed<-do.edgeR(second_processed, meth='Det')
  }else if(third.meth == 'edgeR_DetRate_Cov'){
    third_processed<-do.edgeR(second_processed, meth='Det_cov')
  }else if(third.meth == 'MAST'){
    third_processed<-do.MAST(second_processed, meth='sep')
  }else if(third.meth == 'MAST_Cov'){
    third_processed<-do.MAST(second_processed, meth='cov')
  }else if(third.meth == 'limma'){
    third_processed<-do.limma(second_processed, meth='sep')
  }else if(third.meth == 'limma_Cov'){
    third_processed<-do.limma(second_processed, meth='cov')
  }else if(third.meth == 'REM'){
    third_processed<-do.Meta.ES(second_processed,meth='REM')
  }else if(third.meth == 'FEM'){
    third_processed<-do.Meta.ES(second_processed,meth='FEM')
  }else if(third.meth == 'wFisher'){
    load(paste0(base_dir,filt_dir,'/first_processed/',ct,ifelse(is.null(filter.rate),'',paste0('.',filter.rate,'mincell_filtered')),'_',first.meth,ifelse(ind.analysis,'.ind',''),'_first_processed.RData'))
    weight<-get.weight(data.QC.filtered = first_processed,cell.weight = 'sample',weight.gene=T)
    second_processed[['weight']]<-weight
    second_processed[['two_tailed']]<-two_tailed(second_processed)
    third_processed<-do.Meta.pval(second_processed,meth='wFisher')
  }else if(third.meth == 'Fisher'){
    third_processed<-do.Meta.pval(second_processed,meth='Fisher')
  }
  save(third_processed,file=paste0(base_dir,'data/',filt_dir,'/third_processed/',ct,ifelse(is.null(filter.rate),'',paste0('.',filter.rate,'mincell_filtered')),'_',first.meth,'|',second.meth,'|',third.meth,'_third_processed.RData'))
}
