
rm(list=ls())

library(SingleCellExperiment)
library(scater)
#BiocManager::install("scMerge")
library(scMerge)
library(batchelor)
library(Seurat)
# 

lsdir <- list.dirs('data', recursive=FALSE) 
dir.create('demo_scmerge')

sapply(lsdir,function(x){

  x2 <- gsub('data/','',x)
  dir.create(paste0('demo_scmerge/',x2), showWarnings = FALSE)
  
  # selection <- c('all', 'HVG')
  selection <- c('all')
  sapply(selection, function(s){
  # s<- 'all'  
    
    dir.create(paste0('demo_scmerge/',x2,'/',s), showWarnings = FALSE)
    
    # read data counts and cellinfo
    counts <- read.table(paste0(x,'/counts.txt'), head=T, sep='\t', fill=T)
    cellinfo <- read.table(paste0(x,'/cellinfo.txt'), head=T, sep='\t', fill=T)
    cellinfo <- subset(cellinfo,select=c('Batch','Group'))
    colnames(cellinfo) <- c('batch','cell_type')
    rownames(cellinfo) <- factor(colnames(counts))
    
    counts <- counts[rowSums(counts)>1,]
    
    # remove genes with variance equals to 0
    rv_genes <- which(apply(counts, 1, var)==0) 
    rv_genes_names <- rownames(counts)[rv_genes]
    count_df <- counts[!(rownames(counts) %in% rv_genes_names),]
    count_df <- count_df[,rownames(cellinfo)]
    
    sob <- CreateSeuratObject(counts = count_df, min.cells = 3, min.features = 200)
    sob <- NormalizeData(sob)
    sob <- FindVariableFeatures(sob, selection.method = "vst", nfeatures = 2000)
    hvgenes <- head(VariableFeatures(sob), 1000)
    
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
    
    sce <- SingleCellExperiment(assays = list(counts = count_df),colData = cellinfo)
    sce <- logNormCounts(sce)
    
    result = scSEGIndex(exprs_mat = count_df)
    seg <- rownames(count_df)[which(result$segIdx < .5)]
    print (length(seg))
    print (min(result$segIdx))
    print (max(result$segIdx))
    
    # Run ScMerge
    t1 = Sys.time()
    scMerge_res <- scMerge(
      sce_combine = sce, 
      ctl = seg,
      kmeansK = c(2, 2, 2, 2, 2, 2, 2),
      assay_name = "scMerge_res",
      cell_type = sce$cell_type,
      replicate_prop = 1, verbose=T,
      # marker = rownames(counts))
      marker = hvgenes)
    t2 = Sys.time()
    print(t2-t1)
    
    scMerge_res = runPCA(scMerge_res,
                                exprs_values = "scMerge_res")
    
    scater::plotPCA(
      scMerge_res,
      colour_by = "cell_type",
      shape_by = "batch")
    
    # save the output
    save(scMerge_res, file = paste0('demo_scmerge/',x2,'/',s,"/output.rda"))
    scmerge_norm <- scMerge_res@assays@data$scMerge_res       
    write.table(scmerge_norm, file = paste0('demo_scmerge/',x2,'/',s,"/output.txt"), row.names = T, col.names = T, sep="\t")
  
  })
  
})


