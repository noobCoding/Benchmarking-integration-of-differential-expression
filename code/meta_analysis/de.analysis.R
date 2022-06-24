get.raw<-function(count, annot, ind.analysis){
  if(ind.analysis){
    raw_processed<-list()
    for(b in unique(annot$Patient)){
      raw_processed[[b]][['count']]=as.matrix(count[,annot$Index[which(annot$Patient==b)]])
      raw_processed[[b]][['group']]=annot$Sample_Origin[match(colnames(raw_processed[[b]][['count']]),annot$Index)]
    }
  }else{
    raw_processed<-list()
    raw_processed[['count']]=as.matrix(count)
    raw_processed[['group']]=annot$Sample_Origin[match(colnames(raw_processed[['count']]),annot$Index)]
    raw_processed[['batch']]=annot$Patient[match(colnames(raw_processed[['count']]),annot$Index)]
  }
  return(raw_processed)
}
get.scMerge<-function(count,annot,ind.analysis){
  library(SingleCellExperiment)
  library(scater)
  #BiocManager::install("scMerge")
  library(scMerge)
  library(batchelor)
  library(Seurat)
  
  sample_names<-count%>%colnames()
  cellinfo<-data.frame(Batch=factor(annot$Patient[match(sample_names,annot$Index)]),group=factor(annot$Sample_Origin[match(sample_names,annot$Index)]),row.names = sample_names)
  colnames(cellinfo) <- c('batch','cell_type')
  chunk <- function(x,n){
    vect <- c(1:x)
    num <- ceiling(x/n)
    split(vect,rep(1:num,each=n,len=x))
  }
  
  log_fun <- function(data){
    if(dim(data)[2]>200000){
      l <- lapply(chunk(ncol(data),200000), function(i){log2(data[,i]+1)})
      res <- do.call(cbind,l)
    } else {
      res <- log2(data+1)
    }
    return(res)
  }
  
  sob <- CreateSeuratObject(counts = count)
  sob <- NormalizeData(sob)
  sob <- FindVariableFeatures(sob, selection.method = "vst", nfeatures = 2000)
  hvgenes <- head(VariableFeatures(sob), 1000)
  rm(sob)
  
  sce <- SingleCellExperiment(assays = list(counts = count),colData = cellinfo)
  sce <- logNormCounts(sce)
  result = scSEGIndex(exprs_mat = count)
  seg <- rownames(count)[which(result$segIdx < .5)]
  scMerge_res <- scMerge(
    sce_combine = sce, 
    ctl = seg,
    kmeansK = rep(2,length(unique(cellinfo$Batch))),
    assay_name = "scMerge_res",
    cell_type = sce$cell_type,
    replicate_prop = 1, verbose=T,
    # marker = rownames(counts))
    marker = hvgenes)
  
  scmerge_norm <- scMerge_res@assays@data$scMerge_res 
  
  scmerge_processed<-list()
  scmerge_processed[['count']]=scmerge_norm
  scmerge_processed[['group']]=annot$Sample_Origin[match(scmerge_norm%>%colnames(),annot$Index)]
  scmerge_processed[['batch']]=annot$Patient[match(scmerge_norm%>%colnames(),annot$Index)]
  return(scmerge_processed)
}
get.Seurat<-function(count, annot, ind.analysis){
  library(Seurat)
  
  if(ind.analysis){
    stop('not implemented')
  }else{
    
    sample_names<-count%>%colnames()
    
    cellinfo<-data.frame(Batch=factor(annot$Patient[match(sample_names,annot$Index)]),group=factor(annot$Sample_Origin[match(sample_names,annot$Index)]),row.names = sample_names)
    colnames(cellinfo) <- c('batch','cell_type')
    
    
    seurat_integrated <- CreateSeuratObject(counts = count, project = '', min.cells = 0, min.features = 0, meta.data = cellinfo)
    seurat_integrated<-SCTransform(seurat_integrated,variable.features.n=nrow(seurat_integrated),ncells = ncol(seurat_integrated),return.only.var.genes = F,conserve.memory = F, batch_var='batch', n_genes=NULL,min_cells=0)
    
    seurat_integrated<-seurat_integrated@assays$SCT@data
  
    seurat_processed<-list()
    seurat_processed[['count']]=seurat_integrated
    seurat_processed[['group']]=annot$Sample_Origin[match(seurat_integrated%>%colnames(),annot$Index)]
    seurat_processed[['batch']]=annot$Patient[match(seurat_integrated%>%colnames(),annot$Index)]
    return(seurat_processed)
  }
}
get.voom<-function(count,annot,ind.analysis){
  if(ind.analysis){
    voom_processed=raw_processed=list()
    for(b in unique(annot$Patient)){
      raw_processed[[b]][['count']]=as.matrix(count[,annot$Index[which(annot$Patient==b)]])
      raw_processed[[b]][['group']]=factor(annot$Sample_Origin[match(colnames(raw_processed[[b]][['count']]),annot$Index)])
    }
    for(b in unique(annot$Patient)){
      voom_processed[[b]][['x']]=as.matrix(count[,annot$Index[which(annot$Patient==b)]])
      y_temp<-annot$Sample_Origin[match(colnames(raw_processed[[b]][['count']]),annot$Index)]
      voom_processed[[b]][['y']]=y_temp
      nf <- edgeR::calcNormFactors(voom_processed[[b]]$x, method = 'TMM')
      voom.data <- limma::voom(voom_processed[[b]]$x, design = model.matrix(~factor(voom_processed[[b]]$y)), lib.size = colSums(voom_processed[[b]]$x) * nf)
      voom.data$genes <- rownames(voom_processed[[b]]$x)
      voom_processed[[b]]$count<-voom_processed[[b]]$x
      voom_processed[[b]]$x<-voom.data$E
    }
  }
  return(voom_processed)
}
get.combat<-function(count,annot,ind.analysis){
  library(magrittr)
  library(dplyr)
  library(Seurat)
  library(sva)
  
  sample_names<-count%>%colnames()
  cellinfo<-data.frame(Batch=factor(annot$Patient[match(sample_names,annot$Index)]),group=factor(annot$Sample_Origin[match(sample_names,annot$Index)]),row.names = sample_names)
  col_sums = apply(count,2, sum)
  med_trans = median(col_sums)
  norm_counts = med_trans* scale(count, center=FALSE, scale=col_sums)
  myFilteredData = log(norm_counts + 1)
  
  
  count_df <- myFilteredData
  
  # Run COMBAT
  combat_output = ComBat(dat=as.matrix(count_df), 
                         batch=cellinfo$Batch, 
                         mod=NULL, par.prior=TRUE, prior.plots=FALSE, mean.only=FALSE)
  combat_processed<-list()
  combat_processed[['count']]=combat_output
  combat_processed[['group']]=annot$Sample_Origin[match(combat_output%>%colnames(),annot$Index)]
  combat_processed[['batch']]=annot$Patient[match(combat_output%>%colnames(),annot$Index)]
  return(combat_processed)
}
get.mnn<-function(count,annot,ind.analysis){
  library(scran)
  library(scales)
  require(Rtsne)
  library(Seurat)
  library(batchelor)
  
  sample_names<-count%>%colnames()
  cellinfo<-data.frame(Batch=annot$Patient[match(sample_names,annot$Index)],group=annot$Sample_Origin[match(sample_names,annot$Index)],row.names = sample_names)
  dt <- CreateSeuratObject(counts = count, meta.data = cellinfo)
  dt <- NormalizeData(dt, normalization.method = "LogNormalize")
  t1 = Sys.time()
  out.mnn.total <- batchelor::mnnCorrect(as.matrix(dt@assays$RNA@data),batch=cellinfo$Batch, k=20, sigma=1, cos.norm.in=TRUE, cos.norm.out=TRUE, var.adj=TRUE)
  t2 = Sys.time()
  print(t2-t1)
  
  corre.mnn <- out.mnn.total@assays@data$corrected # @assays[['corrected']]
  mnn_processed<-list()
  mnn_processed[['count']]=corre.mnn
  mnn_processed[['group']]=annot$Sample_Origin[match(corre.mnn%>%colnames(),annot$Index)]
  mnn_processed[['batch']]=annot$Patient[match(corre.mnn%>%colnames(),annot$Index)]
  return(mnn_processed)
}
get.mnn_hai<-function(count,annot,ind.analysis){
  library(scran)
  library(scales)
  require(Rtsne)
  library(Seurat)
  library(batchelor)
  
  sample_names<-count%>%colnames()
  cellinfo<-data.frame(Batch=annot$Patient[match(sample_names,annot$Index)],group=annot$Sample_Origin[match(sample_names,annot$Index)],row.names = sample_names)
  dt <- CreateSeuratObject(counts = count, meta.data = cellinfo)
  dt <- NormalizeData(dt, normalization.method = "LogNormalize")
  t1 = Sys.time()
  out.mnn.total <- batchelor::mnnCorrect(as.matrix(dt@assays$RNA@data),batch=cellinfo$Batch, k=20, sigma=0.1, cos.norm.in=TRUE, cos.norm.out=TRUE, var.adj=TRUE)
  t2 = Sys.time()
  print(t2-t1)
  
  corre.mnn <- out.mnn.total@assays@data$corrected # @assays[['corrected']]
  mnn_processed<-list()
  mnn_processed[['count']]=corre.mnn
  mnn_processed[['group']]=annot$Sample_Origin[match(corre.mnn%>%colnames(),annot$Index)]
  mnn_processed[['batch']]=annot$Patient[match(corre.mnn%>%colnames(),annot$Index)]
  return(mnn_processed)
}
get.mnn_opt<-function(count,annot,ind.analysis){
  library(scran)
  library(scales)
  require(Rtsne)
  library(Seurat)
  library(batchelor)
  
  sample_names<-count%>%colnames()
  cellinfo<-data.frame(Batch=annot$Patient[match(sample_names,annot$Index)],group=annot$Sample_Origin[match(sample_names,annot$Index)],row.names = sample_names)
  dt <- CreateSeuratObject(counts = count, meta.data = cellinfo)
  dt <- NormalizeData(dt, normalization.method = "LogNormalize")
  t1 = Sys.time()
  out.mnn.total <- batchelor::mnnCorrect(as.matrix(dt@assays$RNA@data),batch=cellinfo$Batch, k=20, sigma=0.1, cos.norm.in=TRUE, cos.norm.out=TRUE, var.adj=TRUE, correct.all = T, auto.merge = T)
  t2 = Sys.time()
  print(t2-t1)
  
  corre.mnn <- out.mnn.total@assays@data$corrected # @assays[['corrected']]
  mnn_processed<-list()
  mnn_processed[['count']]=corre.mnn
  mnn_processed[['group']]=annot$Sample_Origin[match(corre.mnn%>%colnames(),annot$Index)]
  mnn_processed[['batch']]=annot$Patient[match(corre.mnn%>%colnames(),annot$Index)]
  return(mnn_processed)
  
  
  
}
get.limma_bec<-function(count,annot,ind.analysis){
  library(SingleCellExperiment)
  library(limma)
  library(Seurat)
  library(dplyr)
  library(magrittr)
  limma_bec_processed<-list()
  sample_names<-count%>%colnames()
  cellinfo<-data.frame(Batch=annot$Patient[match(sample_names,annot$Index)],group=annot$Sample_Origin[match(sample_names,annot$Index)],row.names = sample_names)
  
  ct <- CreateSeuratObject(counts = count, meta.data = cellinfo, project = "Limma", min.cells = 0)
  ct <- NormalizeData(object = ct, normalization.method = "LogNormalize", scale.factor = 1e4)
  ct <- ScaleData(object = ct)
  lm_df <- ct@assays$RNA@data
  
  # run limma
  lm_df <- limma::removeBatchEffect(as.matrix(lm_df), factor(cellinfo$Batch))
  limma_bec_processed[['count']]=lm_df
  limma_bec_processed[['group']]=annot$Sample_Origin[match(lm_df%>%colnames(),annot$Index)]
  limma_bec_processed[['batch']]=annot$Patient[match(lm_df%>%colnames(),annot$Index)]
  
  return(limma_bec_processed)
}
get.LogNormalize<-function(count, annot, ind.analysis){
  library(Seurat)
  
  if(ind.analysis){
    message('LogNormalize ind run')
    LogNormalize_processed<-list()
    for(b in unique(annot$Patient)){
      LogNormalize_processed[[b]][['count']]=count[,annot$Index[which(annot$Patient==b)]]%>%as.matrix()%>%NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000)%>%as.matrix()
      LogNormalize_processed[[b]][['group']]=annot$Sample_Origin[match(LogNormalize_processed[[b]][['count']]%>%colnames(),annot$Index)]
      
    }
  }else{
    message('LogNormalize run')
    LogNormalize_processed<-list()
    LogNormalize_processed[['count']]=count%>%as.matrix()%>%NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000)%>%as.matrix()
    LogNormalize_processed[['group']]=annot$Sample_Origin[match(LogNormalize_processed[['count']]%>%colnames(),annot$Index)]
    LogNormalize_processed[['batch']]=annot$Patient[match(LogNormalize_processed[['count']]%>%colnames(),annot$Index)]
  }
  return(LogNormalize_processed)
}
one_tailed<-function(x,ind.tail,logFC){
  if(ind.tail=='high'&logFC>0 || ind.tail=='low'&logFC<=0){
    # one_tailed_p<-(1-0.5*x)
    one_tailed_p<-0.5*x
  }else{
    # one_tailed_p<-0.5*x
    one_tailed_p<-(1-0.5*x)
  }
  return(one_tailed_p)
}

get.raw.second<-function(first_processed){
  return(first_processed)
}
get.pval.ind<-function(first_processed, using){
  library(DESeq2)
  library(edgeR)
  library(Seurat)
  library(dplyr)
  library(magrittr)
  ind.Res<-list()
  weight=rep(1,length(first_processed))
  for (wi in 1:length(weight)) {
    weight[wi] = sqrt(length(first_processed[[wi]]$group))
  }
  weight<-weight/sum(weight)
  ind.Res[['weight']]=weight
  tails<-c('high','low')
  p=matrix(nrow=nrow(first_processed[[1]]$count),ncol=length(first_processed))
  rownames(p)=rownames(first_processed[[1]]$count)
  colnames(p)=names(first_processed)
  ind.Res[[tails[1]]]=ind.Res[[tails[2]]]=ind.Res[['logFC']]=p
  if(using=='DESeq2'){
    for(ind.study in names(first_processed)){
      count_df<-first_processed[[ind.study]]$count
      count_df <- round(count_df, 0) + 1 # pseudo count prevents every gene includes at least one 0
      # group=first_processed[[ind.study]]$group
      group=factor(first_processed[[ind.study]]$group)
      dds <- DESeqDataSetFromMatrix(countData = count_df, colData = data.frame(group=factor(group)), design = ~group)
      dds <- DESeq2::DESeq(dds) #, fitType ='mean')
      # res<-lfcShrink(dds,coef='group_tLung_vs_nLung', type="apeglm", lfcThreshold=0)
      res<-lfcShrink(dds,coef=2, type="apeglm", lfcThreshold=0)
      pval=res$pvalue
      logFC=res$log2FoldChange
      
      ind.Res[['logFC']][,ind.study]<-logFC
      for(ind.tail in c('high', 'low')){
        tailed=switch (ind.tail,
                       "low" = "left",
                       "high" = "right"
        )
        ind.Res[[ind.tail]][,ind.study]<-unlist(mapply(function(x,y) one_tailed(x=x,logFC=y,ind.tail=ind.tail),pval,logFC))
      }
    }
  }else if(using=='edgeR'){
    for(ind.study in names(first_processed)){
      count_df<-first_processed[[ind.study]]$count
      group=first_processed[[ind.study]]$group
      y <- DGEList(counts=count_df, group=factor(group))
      y <- calcNormFactors(y)
      design <- model.matrix(~factor(group))
      rownames(design) <- colnames(y)
      y <- estimateDisp(y, design, robust=TRUE)
      fit <- glmQLFit(y, design, robust=TRUE, prior.df = 0)
      qlf <- glmQLFTest(fit)
      logFC<-qlf$table$logFC
      pval<-qlf$table$PValue
      ind.Res[['logFC']][,ind.study]<-logFC
      for(ind.tail in c('high', 'low')){
        tailed=switch (ind.tail,
                       "low" = "left",
                       "high" = "right"
        )
        ind.Res[[ind.tail]][,ind.study]<-unlist(mapply(function(x,y) one_tailed(x=x,logFC=y,ind.tail=ind.tail),pval,logFC))
      }
    }
  }else if(using=='findmarkers'){
    for(ind.study in names(first_processed)){
      processed.ind<-first_processed[[ind.study]]
      markers<-do.findmarker(processed.ind)
      logFC<-markers$avg_log2FC[match(rownames(p),rownames(markers))]
      pval<-markers$p_val[match(rownames(p),rownames(markers))]
      qval<-markers$p_val_adj[match(rownames(p),rownames(markers))]
      ind.Res[['logFC']][,ind.study]<-logFC
      for(ind.tail in c('high', 'low')){
        tailed=switch (ind.tail,
                       "low" = "left",
                       "high" = "right"
        )
        #save pval for ind analysis
        #We'll apply bonferroni(findmarker default correction) later.
        ind.Res[[ind.tail]][,ind.study]<-pval
      }
    }
  }else if(using=='modt'){
    MetaDE.Res<-list()
    MetaDE.Res.ind<-list()
    ndata=length(first_processed)
    MetaDE.test<-first_processed
    for(i in 1:length(MetaDE.test)){
      control<-sort(first_processed[[1]]$y%>%unique())[1]
      MetaDE.test[[i]]$y=ifelse(MetaDE.test[[i]]$y==control,0,1)
    }
    ind.Res1<-MetaDE::ind.analysis(MetaDE.test,ind.method=rep("modt",ndata),nperm=300,tail='abs')
    ind.Res2<-MetaDE::ind.analysis(MetaDE.test,ind.method=rep("modt",ndata),nperm=300,tail='high')
    ind.Res3<-MetaDE::ind.analysis(MetaDE.test,ind.method=rep("modt",ndata),nperm=300,tail='low')
    for(ind.study in names(first_processed)){
      #assign only signs for modt
      logFC<-ifelse(ind.Res2$p[,ind.study]<ind.Res3$p[,ind.study],1,-1)
      pval<-ind.Res1$p[,ind.study]
      ind.Res[['logFC']][,ind.study]<-logFC
      for(ind.tail in c('high', 'low')){
        tailed=switch (ind.tail,
                       "low" = "left",
                       "high" = "right"
        )
        ind.Res[[ind.tail]][,ind.study]<-unlist(mapply(function(x,y) one_tailed(x=x,logFC=y,ind.tail=ind.tail),pval,logFC))
      }
    }
  }else if(using=='limma_trend'){
    for(ind.study in names(first_processed)){
      library(limma)
      library(edgeR)
      count_df<-first_processed[[ind.study]]$count
      count_df[is.na(count_df)] = 0.
      cellinfo<-data.frame(group=first_processed[[ind.study]]$group%>%factor(), 
                           row.names = colnames(count_df))
      design <- model.matrix(~group, data=cellinfo)
      
      lmfit <- lmFit(count_df, design)
      
      lmfit <- eBayes(lmfit, trend=TRUE, robust = TRUE)
      res <- topTable(lmfit, n = Inf, adjust.method = "BH", coef = 2)
      logFC<-res$logFC
      pval<-res$P.Value
      ind.Res[['logFC']][,ind.study]<-logFC
      for(ind.tail in c('high', 'low')){
        tailed=switch (ind.tail,
                       "low" = "left",
                       "high" = "right"
        )
        ind.Res[[ind.tail]][,ind.study]<-unlist(mapply(function(x,y) one_tailed(x=x,logFC=y,ind.tail=ind.tail),pval,logFC))
      }
    }
  }
  return(ind.Res)
}

get.ES.ind<-function(first_processed, using){
  library(DESeq2)
  library(MetaDE)
  if(using=='DESeq2'){
    ind.Res<-list()
    temp_M<-matrix(nrow=nrow(first_processed[[1]]$count),ncol=length(first_processed))
    rownames(temp_M)=rownames(first_processed[[1]]$count)
    colnames(temp_M)=names(first_processed)
    ind.Res$ES=ind.Res$SE=temp_M
    for(ind.study in names(first_processed)){
      count_df<-first_processed[[ind.study]]$count
      count_df <- round(count_df, 0) + 1 # pseudo count prevents every gene includes at least one 0
      group=first_processed[[ind.study]]$group
      dds <- DESeqDataSetFromMatrix(countData = count_df, colData = data.frame(group=factor(group)), design = ~group)
      dds <- DESeq2::DESeq(dds) #, fitType ='mean')
      res<-lfcShrink(dds,coef=2, type="apeglm", lfcThreshold=0)
      ind.Res$ES[,ind.study]=res$log2FoldChange
      ind.Res$SE[,ind.study]=res$lfcSE
    }
    return(ind.Res)
  }else if(using=='metaDE'){
    if(is.null(first_processed[[1]][['x']])){
      first_processed_pseudo<-lapply(first_processed,FUN=function(x){y<-list()
      y$count<-x$count
      y$group<-x$group
      return(y)})
    }else{
      first_processed_pseudo<-lapply(first_processed,FUN=function(x){y<-list()
      y$count<-x$x
      y$group<-x$y
      return(y)})
    }
    ind.Res<-ind.cal.ES.bk(first_processed_pseudo, vs=first_processed_pseudo[[1]]$group%>%unique()%>%sort(decreasing = T))
    
    ind.Res$SE=sqrt(ind.Res$Var)
    colnames(ind.Res$ES)=colnames(ind.Res$Var)=colnames(ind.Res$SE)=names(first_processed_pseudo)
  }
  return(ind.Res)
}
do.findmarker<-function(second_processed){
  library(Seurat)
  count_matrix<-second_processed$count
  dt<-CreateSeuratObject(counts=count_matrix)
  dt[["percent.mt"]] <- PercentageFeatureSet(dt, pattern = "^MT-")
  dt<-NormalizeData(dt, normalization.method = "LogNormalize", scale.factor = 10000)
  dt<-Seurat::ScaleData(dt,do.scale=F,do.center=T, scale.max=10)
  
  Idents(dt)=second_processed$group
  markers<-FindMarkers(dt,ident.1 = sort(unique(Idents(dt)))[2],ident.2 = sort(unique(Idents(dt)))[1],min.pct = 0,logfc.threshold = 0,min.cells.feature = 0,min.cells.group = 0)
  return(markers)
}
do.findmarker.prenorm<-function(second_processed){
  library(Seurat)
  count_matrix<-second_processed$count
  dt<-CreateSeuratObject(counts=count_matrix)
  Idents(dt)=second_processed$group
  # markers<-FindMarkers(dt,ident.1 = "tLung",ident.2 = 'nLung',min.pct = 0,logfc.threshold = 0)
  markers<-FindMarkers(dt,ident.1 = sort(unique(Idents(dt)))[2],ident.2 = sort(unique(Idents(dt)))[1],min.pct = 0,logfc.threshold = 0,min.cells.feature = 0,min.cells.group = 0)
  return(markers)
}
do.MAST<-function(second_processed, meth){
  library(MAST)
  library(Seurat)
  count_matrix<-second_processed$count
  cellinfo<-data.frame(group=second_processed$group%>%factor(), 
                       batch=second_processed$batch%>%factor(), 
                       row.names = colnames(count_matrix))
  dt<-CreateSeuratObject(counts=count_matrix, meta.data=cellinfo, project='MAST', min.cells=0)
  dt<-NormalizeData(dt, normalization.method = "LogNormalize", scale.factor = 10000)
  Idents(dt)=second_processed$group
  if(meth=='sep'){
    markers <- FindMarkers(object = dt, ident.1 = sort(unique(Idents(dt)))[2],ident.2 = sort(unique(Idents(dt)))[1],
                       test.use='MAST', logfc.threshold = 0,
                       min.cells.feature = 0,
                       min.cells.group = 0,
                       min.pct = 0,
                       only.pos = F)
  }else if(meth=='cov'){
    markers <- FindMarkers(object = dt, ident.1 = sort(unique(Idents(dt)))[2],ident.2 = sort(unique(Idents(dt)))[1],
                       test.use='MAST', logfc.threshold = 0,
                       latent.vars=c('batch'),
                       min.cells.feature = 0,
                       min.cells.group = 0,
                       min.pct = 0,
                       only.pos = F)
  }
  return(markers)
}
do.limma<-function(second_processed,meth){
  library(edgeR)
  library(limma)
  count_df<-second_processed$count
  cellinfo<-data.frame(group=second_processed$group%>%factor(), 
                       batch=second_processed$batch%>%factor(), 
                       row.names = colnames(count_df))
  
  nf <- edgeR::calcNormFactors(count_df, method = 'TMM')
  if(meth=='sep'){
    mod <- model.matrix(~group, data=cellinfo)
  }else if(meth=='cov'){
    mod <- model.matrix(~group+batch, data=cellinfo)
  }
  
  
  voom.data <- limma::voom(count_df, design = mod, lib.size = colSums(count_df) * nf)
  voom.data$genes <- rownames(count_df)
  voom.fitlimma <- limma::lmFit(voom.data, design = mod)
  voom.fitbayes <- limma::eBayes(voom.fitlimma)
  voom.pvalues <- voom.fitbayes$p.value[, 2]
  voom.adjpvalues <- p.adjust(voom.pvalues, method = 'BH')
  voom.logFC <- voom.fitbayes$coefficients[, 2]
  voom.score <- 1 - voom.pvalues
  result.table <- data.frame('pvalue' = voom.pvalues, 'adjpvalue' = voom.adjpvalues, 'logFC' = voom.logFC, 'score' = voom.score)
  rownames(result.table) <- rownames(count_df)
  return(result.table)
}
do.DESeq2<-function(second_processed, meth){
  
  library(Seurat)
  library(DESeq2)
  count_matrix<-second_processed$count
  
  cellinfo<-data.frame(group=second_processed$group%>%factor(), 
                       batch=second_processed$batch%>%factor(), 
                       row.names = colnames(count_matrix))
  
  
  
  if(meth=='cov'){
    count_matrix <- round(count_matrix, 0) + 1
    dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = cellinfo, design = ~group+batch)
    dds <- DESeq2::DESeq(dds) #, fitType ='mean')
    res <- lfcShrink(dds, coef =2 , type="apeglm", lfcThreshold=0)
    result.table <- data.frame('pvalue' = res$pvalue, 'adjpvalue' = res$padj, 'logFC' = res$log2FoldChange)
    rownames(result.table) <- rownames(dds)
  }else if(meth=='sep'){
    count_matrix <- round(count_matrix, 0) + 1
    dds_sep <- DESeqDataSetFromMatrix(countData = count_matrix, colData = cellinfo, design = ~group)
    dds_sep <- DESeq2::DESeq(dds_sep) #, fitType ='mean')
    res_sep <- lfcShrink(dds_sep, coef=2, type="apeglm", lfcThreshold=0)
    result.table <- data.frame('pvalue' = res_sep$pvalue, 'adjpvalue' = res_sep$padj, 'logFC' = res_sep$log2FoldChange)
    rownames(result.table) <- rownames(dds_sep)
  }
  return(result.table)
}
do.edgeR<-function(second_processed,meth){
  library(Seurat)
  library(edgeR)
  count_matrix<-second_processed$count
  
  cellinfo<-data.frame(group=second_processed$group%>%factor(), 
                       batch=second_processed$batch%>%factor(), 
                       row.names = colnames(count_matrix))
  
  
  y <- DGEList(counts=count_matrix, group=cellinfo$group)
  y <- calcNormFactors(y)
  cellGroup <- factor(cellinfo$group)
  cellBatch <- factor(cellinfo$batch)
  cdr <- scale(colMeans(count_matrix > 0))
  if(meth=='cov'){
    # design <- model.matrix(~cellBatch+cellGroup)
    design <- model.matrix(~cellGroup+cellBatch)
    rownames(design) <- colnames(y)
    y <- estimateDisp(y, design, robust=TRUE)
    
    fit <- glmQLFit(y, design, robust=TRUE, prior.df = 0)
    # qlf <- glmQLFTest(fit)
    qlf <- glmQLFTest(fit, coef=2)
    
    FDR<-p.adjust(qlf$table$PValue,method = "BH")
    qlf$table$FDR <- FDR
    
    result.table <- data.frame('pvalue' = qlf$table$PValue, 'adjpvalue' = qlf$table$FDR, 'logFC' = qlf$table$logFC)
    rownames(result.table) <- rownames(qlf)
  }else if(meth=='sep'){
    design <- model.matrix(~cellGroup)
    rownames(design) <- colnames(y)
    y <- estimateDisp(y, design, robust=TRUE)
    
    fit <- glmQLFit(y, design, robust=TRUE, prior.df = 0)
    # qlf <- glmQLFTest(fit)
    qlf <- glmQLFTest(fit, coef=2)
    
    FDR<-p.adjust(qlf$table$PValue,method = "BH")
    qlf$table$FDR <- FDR
    
    result.table <- data.frame('pvalue' = qlf$table$PValue, 'adjpvalue' = qlf$table$FDR, 'logFC' = qlf$table$logFC)
    rownames(result.table) <- rownames(qlf)
  }else if(meth=='Det'){
    # design <- model.matrix(~cdr + cellGroup)
    design <- model.matrix(~cellGroup+cdr)
    rownames(design) <- colnames(y)
    y <- estimateDisp(y, design, robust=TRUE)
    
    fit <- glmQLFit(y, design, robust=TRUE, prior.df = 0)
    # qlf <- glmQLFTest(fit)
    qlf <- glmQLFTest(fit, coef=2)
    
    FDR<-p.adjust(qlf$table$PValue,method = "BH")
    qlf$table$FDR <- FDR
    
    result.table <- data.frame('pvalue' = qlf$table$PValue, 'adjpvalue' = qlf$table$FDR, 'logFC' = qlf$table$logFC)
    rownames(result.table) <- rownames(qlf)
  }else if(meth=='Det_cov'){
    # design <- model.matrix(~cdr + cellBatch + cellGroup)
    design <- model.matrix(~cellGroup+cdr + cellBatch)
    rownames(design) <- colnames(y)
    
    y <- estimateDisp(y, design, robust=TRUE)
    # y$common.dispersion
    
    fit <- glmQLFit(y, design, robust=TRUE, prior.df = 0)
    # qlf <- glmQLFTest(fit)
    qlf <- glmQLFTest(fit, coef=2)
    
    FDR<-p.adjust(qlf$table$PValue,method = "BH")
    qlf$table$FDR <- FDR
    
    
    result.table <- data.frame('pvalue' = qlf$table$PValue, 'adjpvalue' = qlf$table$FDR, 'logFC' = qlf$table$logFC)
    rownames(result.table) <- rownames(qlf)
  }
  return(result.table)
}
do.Meta.ES<-function(second_processed,meth){
  Meta.Res<-list()
  for(ind.tail in c('high', 'low','abs')){
    if(meth=='FEM'){
      Meta.Res[[ind.tail]]<-run_FEM(effect_size=second_processed$ES,standard_error=second_processed$SE,tail=ind.tail)
    }else if(meth=='REM'){
      Meta.Res[[ind.tail]]<-run_REM(effect_size=second_processed$ES,standard_error=second_processed$SE,tail=ind.tail)
    }
  }
  return(Meta.Res)
}
do.Meta.pval<-function(second_processed,meth){
  library(MetaDE)
  library(metapro)
  Meta.Res<-list()
  Meta.Res[['logFC_table']]=second_processed[['logFC']]
  if(meth=='wFisher'){
    
    Meta.Res[["logFC"]]<-rowSums(second_processed[['logFC']]*second_processed[['weight']])
  }else if(meth=='Fisher'){
    # Meta.Res[["logFC"]]<-apply(second_processed[['logFC']],1,mean)
  }
  
  if(meth=='wFisher'){
    x<-list()
    x$p<-second_processed[["two_tailed"]]
    x$fc<-second_processed[["logFC"]]
    Meta.Res[['two_tailed']]<-run_wFisher(x,weight=second_processed[["weight"]], weight.cell='sqrt.sample',weight.gene=T)
  }else{
    for(ind.tail in c('high', 'low')){
      tailed=switch (ind.tail,
                     "low" = "left",
                     "high" = "right"
      )
      x<-list()
      
      
      if(meth=='Fisher'){
        x$p<-second_processed[[ind.tail]]
        Meta.Res[[ind.tail]]<-MetaDE.pvalue(x,meta.method=c("Fisher"),asymptotic = T)
      }
    }
    #assign signs to logFC of Fisher
    logFC=Meta.Res[["logFC"]]=ifelse(Meta.Res[['high']][['meta.analysis']][['pval']]<0.5,1,-1)
    
    two_tailed_Fisher=high=Meta.Res[['high']][['meta.analysis']][['pval']]
    low=Meta.Res[['low']][['meta.analysis']][['pval']]
    
    two_tailed_Fisher[which(logFC>0)]=high[which(logFC>0)]*2
    two_tailed_Fisher[which(logFC<=0)]=low[which(logFC<=0)]*2
    # high[which(logFC>0)]=high[which(logFC>0)]*2
    # high[which(logFC<0)]=(1-high[which(logFC<0)])*2
    # Meta.Res[['two_tailed']][['pval']]=high
    Meta.Res[['two_tailed']][['pval']]=two_tailed_Fisher
    Meta.Res[['two_tailed']][['FDR']]=p.adjust(two_tailed_Fisher,method='BH')
  }
  
  
  return(Meta.Res)
}

do.limma_trend.raw<-function(second_processed, meth='sep'){
  library(limma)
  library(edgeR)
  count_df<- second_processed$count
  count_df[is.na(count_df)] = 0.
  cellinfo<-data.frame(group=second_processed$group%>%factor(), 
                       batch=second_processed$batch%>%factor(), 
                       row.names = colnames(count_df))
  
  
  # count_df <- count_df + abs(min(count_df))
  
  ## Convert to an edgeR object
  dgeObj <- DGEList(count_df)
  
  ## Perform TMM normalisation
  dgeObj <- calcNormFactors(dgeObj)
  
  # design <- model.matrix(~Group+Batch, data=cellinfo)
  
  logCPM <- cpm(dgeObj, log=TRUE, prior.count=3)
  if(meth=='sep'){
    design <- model.matrix(~group, data=cellinfo)
    lmfit <- lmFit(logCPM, design)
    
    lmfit <- eBayes(lmfit, trend=TRUE, robust = TRUE)
    res <- topTable(lmfit, n = Inf, adjust.method = "BH", coef = 2)
  }else if(meth=='cov'){
    design <- model.matrix(~group+batch, data=cellinfo)
    lmfit <- lmFit(logCPM, design)
    
    lmfit <- eBayes(lmfit, trend=TRUE, robust = TRUE)
    res <- topTable(lmfit, n = Inf, adjust.method = "BH", coef = 2)
  }
  
  
  result.table <- data.frame('pvalue' = res$P.Value, 
                             'adjpvalue' = res$adj.P.Val, 
                             'logFC' = res$logFC,
                             row.names = rownames(res))
  return(result.table)
}
do.limma_trend<-function(second_processed){
  library(limma)
  count_df<- second_processed$count
  count_df[is.na(count_df)] = 0.
  cellinfo<-data.frame(group=second_processed$group%>%factor(), 
                       batch=second_processed$batch%>%factor(), 
                       row.names = colnames(count_df))
  
  
  # count_df <- count_df + abs(min(count_df))
  
  ## Convert to an edgeR object
  # dgeObj <- DGEList(count_df)
  
  ## Perform TMM normalisation
  # dgeObj <- calcNormFactors(dgeObj)
  
  # design <- model.matrix(~Group+Batch, data=cellinfo)
  design <- model.matrix(~group, data=cellinfo)
  
  
  lmfit <- lmFit(count_df, design)
  
  lmfit <- eBayes(lmfit, trend=TRUE, robust = TRUE)
  res <- topTable(lmfit, n = Inf, adjust.method = "BH", coef = 2)
  
  result.table <- data.frame('pvalue' = res$P.Value, 
                             'adjpvalue' = res$adj.P.Val, 
                             'logFC' = res$logFC,
                             row.names = rownames(res))
  
  return(result.table)
}
run_FEM<-function(effect_size, standard_error, tail='abs'){
  weight<-1/standard_error^2
  Avg<-rowSums(weight*effect_size)/rowSums(weight)
  Var<-1/rowSums(weight)
  Z<-Avg/sqrt(Var)
  if (tail=="low"){
    z.p<-pnorm(Z)
  }else if(tail=="high"){
    z.p<-(1-pnorm(Z))
  }else if(tail=="abs"){
    z.p<-2*(1-pnorm(abs(Z)))
  }else{
    stop("Tail should be low, high or abs")
  }
  qval<-p.adjust(z.p,method="BH")
  res<-list(mu.hat=Avg,mu.var=Var,zval=Z,pval=z.p,FDR=qval)
  return(res)
}

run_REM<-function(effect_size, standard_error, tail='abs'){
  weight<-1/standard_error^2
  k<-ncol(effect_size)
  Q.val<-get.Q(effect_size,standard_error)
  tau.sq<-get.tau.sq(Q.val,standard_error,k)
  weight.rem<-1/((1/weight)+tau.sq)	
  Avg<-rowSums(weight.rem*effect_size)/rowSums(weight.rem)
  Var<-1/rowSums(weight.rem)
  Qpval <- pchisq(Q.val, df = k - 1, lower.tail = FALSE)
  Z<-Avg/sqrt(Var)
  if (tail=="low"){
    z.p<-pnorm(Z)
  }else if(tail=="high"){
    z.p<-(1-pnorm(Z))
  }else if(tail=="abs"){
    z.p<-2*(1-pnorm(abs(Z)))
  }else{
    stop("Tail should be low, high or abs")
  }
  qval<-p.adjust(z.p,method="BH")
  res<-list(mu.hat=Avg,mu.var=Var,Qval=Q.val,Qpval=Qpval,tau2=tau.sq,zval=Z,pval=z.p,FDR=qval)
  return(res)
}

get.Q<-function(effect_size,standard_error){
  weight <- 1/standard_error^2
  Avg<-rowSums(weight*effect_size)/rowSums(weight)
  Q <- rowSums(weight*(effect_size-Avg)^2)
  return(Q)
}
get.tau.sq<-function(Q,standard_error,k){
  weight<-1/standard_error^2
  tau.sq<- (Q-(k-1))/(rowSums(weight)-(rowSums(weight^2)/rowSums(weight)))
  tau.sq<-sapply(tau.sq,FUN=function(x)max(x,0))
  return(tau.sq)
}

ind.cal.ES.bk<-function(x, vs){
  K <- length(x)
  res <- get.ES.bk(x, vs)
  if (is.null(names(x))) {
    colnames(res$ES) <- colnames(res$Var) <- paste("dataset", 
                                                   1:K, sep = "")
  }else {
    colnames(res$ES) <- colnames(res$Var) <- names(x)
  }
  result <- list(ES = res$ES, Var = res$Var)
  return(result)  
}
get.ES.bk<-function(x, vs){
  K <- length(x)
  ES.m <- Var.m <- N <- n <- NULL
  for (k in 1:K) {
    y <- x[[k]][[1]]
    l <- x[[k]][[2]]
    temp <- cal.ES.bk(y, l, vs=vs)
    ES.m <- cbind(ES.m, temp[, "dprime"])
    Var.m <- cbind(Var.m, temp[, "vardprime"])
    N <- c(N, length(l))
    n <- c(n, table(l))
  }
  rownames(ES.m) <- rownames(y)
  rownames(Var.m) <- rownames(y)
  colnames(ES.m) <- names(x)
  colnames(Var.m) <- names(x)
  res <- list(ES = ES.m, Var = Var.m)
  return(res)
}
cal.ES.bk<-function (y, l, vs){
  library(gmp)
  l <- unclass(factor(l))
  n <- table(factor(l))
  ind <- diag(rep(1, length(n)))[l, ]
  ym <- y %*% ind %*% diag(1/n)
  ntilde <- 1/sum(1/n)
  m = sum(n) - 2
  # cm = gamma(m/2)/(sqrt(m/2) * gamma((m - 1)/2))
  cm = as.numeric(factorialZ((m/2)-1)/(sqrt(m/2) * factorialZ(((m - 1)/2)-1)))
  ind <- diag(rep(1, length(n)))[l, ]
  ym <- y %*% ind %*% diag(1/n)
  
  s <- sqrt((1/(sum(n)- 2)) * ((((y-ym[,l])^2)%*%ind)%*%diag(1/(n - 1))%*%(n - 1)))
  
  # s <- sqrt((1/(sum(n) - 2) * ((y^2 %*% ind) %*% diag(1/(n - 1)) - ym^2 %*% diag(n/(n - 1))) %*% (n - 1)))
  # d <- (ym[, 2] - ym[, 1])/s
  d <- (ym[, match(vs,levels(l))[1]] - ym[, match(vs,levels(l))[2]])/s
  
  d[which(s==0)]=Inf
  d[which((ym[, match(vs,levels(l))[1]] - ym[, match(vs,levels(l))[2]])==0)]=0
  
  dprime = d - 3 * d/(4 * (sum(n) - 2) - 1)
  terme1 = 1/ntilde
  vard = terme1 + d^2 * (terme1 * ntilde - 1/cm^2)
  vardprime = sum(1/n) + dprime^2/(2 * sum(n))
  result = cbind(dprime, vardprime)
  colnames(result) = c("dprime", "vardprime")
  rownames(result) <- rownames(y)
  return(result)
}
run_wFisher<-function(x, weight, weight.cell='sample',weight.gene=T){
  weight.method=paste0(weight.cell,ifelse(weight.gene, '', '+gene'),'_weights')
  meth.name=paste0(weight.method,'+wFisher')
  K <- ncol(x$p)
  meta.res <- list(stat = NA, pval = NA, FDR = NA, direction=NA, AW.weight = NA)
  meta.res$direction <- meta.res$stat <- meta.res$pval <- meta.res$FDR <- matrix(NA,nrow(x$p), 1)
  
  temp<-get.wFisher.bk(p=x$p, weight = weight ,fc=x$fc)
  
  meta.res$stat[,1]=temp$stat
  meta.res$pval[, 1] <- temp$pval
  meta.res$FDR[, 1] <- temp$FDR
  meta.res$direction[,1]<-temp$direction
  colnames(meta.res$direction)<-colnames(meta.res$stat) <- colnames(meta.res$pval) <- colnames(meta.res$FDR) <- meth.name
  rownames(meta.res$direction)<-rownames(meta.res$stat) <- rownames(meta.res$pval) <- rownames(meta.res$FDR) <- rownames(x$p)
  return(meta.res)
}
# {
#   metapro::wFisher(second_processed$high[3,],weight=weight[3,])$p*2
#   metapro::wFisher(second_processed$low[3,],weight=weight[3,])
#   metapro::wFisher(second_processed$twotailed[3,],weight=weight[3,],is.onetail=F,eff.sign=c(sign(second_processed$logFC)[3,]))$overall.eff.direction
# }
get.wFisher.bk<-function(p,weight,fc){
  k <- ncol(p)
  direction<-pval <- stat <- rep(NA, nrow(p))
  rnum <- 1:nrow(p)
  # pval[rnum] <- unlist(apply(p[rnum, ], 1, metapro::wFisher, weight = weight[rnum,], is.onetail = T))
  # pval[rnum] <- unlist(apply(p[rnum, ], 1, metapro::wFisher, weight = weight[1,], is.onetail = T))
  # pval[rnum]<-unlist(sapply(X=rnum, FUN=function(x){metapro::wFisher(p[x,],weight=weight[x,],is.onetail = T)}))
  # for(x in rnum){
  #   print(x)
  #   metapro::wFisher(p[x,],weight=weight[x,],is.onetail = F,eff.sign = c(sign(fc)[x,]))$p
  # }
  pval[rnum]<-unlist(sapply(X=rnum, FUN=function(x){metapro::wFisher(p[x,],weight=ifelse(all(weight[x,]==0),list(rep(1,ncol(weight))),list(weight[x,]))[[1]],is.onetail = F,eff.sign = c(sign(fc)[x,]))$p}))
  direction[rnum]<-unlist(sapply(X=rnum, FUN=function(x){metapro::wFisher(p[x,],weight=ifelse(all(weight[x,]==0),list(rep(1,ncol(weight))),list(weight[x,]))[[1]],is.onetail = F,eff.sign = c(sign(fc)[x,]))$overall.eff.direction}))
  qval <- p.adjust(pval, method = "BH")
  res <- list(stat = pval, pval = pval, FDR = qval, direction=direction)
  
  names(res$direction)<-names(res$stat) <- names(res$pval) <- names(res$FDR) <- rownames(p)
  return(res)
}
two_tailed<-function(second_processed){
  logFC<-second_processed$logFC
  high<-second_processed$high
  # low<-second_processed$low
  high[which(logFC>0,arr.ind = T)]=high[which(logFC>0,arr.ind = T)]*2
  high[which(logFC<0,arr.ind = T)]=(1-high[which(logFC<0,arr.ind = T)])*2
  # low[which(logFC>0,arr.ind = T)]=(1-low[which(logFC>0,arr.ind = T)])*2
  # low[which(logFC<0,arr.ind = T)]=low[which(logFC<0,arr.ind = T)]*2
  return(high)
}
# two_tailed_for_ef.meta<-function(third_processed){
#   high<-third_processed$high$pval
#   mu.hat<-third_processed$high$mu.hat
#   high[which(mu.hat>0)]=high[which(mu.hat>0)]*2
#   high[which(mu.hat<0)]=(1-high[which(mu.hat<0)])*2
#   return(high)
# }
get.weight<-function(data.weight, cell.weight='sqrt', weight.gene=F){
  gene.weight<-sapply(data.weight,FUN = function(x){ (rowSums(x$count!=0)) })
  weight.temp<-sapply(data.weight,FUN = function(x){ rowSums(!is.na(x$count)) })
  
  
  if(weight.gene==T){
    res_weight=gene.weight
  }else{
    res_weight=weight.temp
  }
  if(cell.weight=='sqrt'){
    # res_weight<-sqrt(res_weight)/rowSums(sqrt(res_weight))
    res_weight<-sqrt(res_weight)/ifelse(rowSums(sqrt(res_weight))==0,1,rowSums(sqrt(res_weight)))
  }else if(cell.weight=='sample'){
    # res_weight<-res_weight/rowSums(res_weight)
    res_weight<-res_weight/ifelse(rowSums(res_weight)==0,1,rowSums(res_weight))
  }
  return(res_weight)
}

