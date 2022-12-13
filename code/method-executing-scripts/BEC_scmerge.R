scSEGIndex_used<-function (exprs_mat, cell_type = NULL, BPPARAM = BiocParallel::SerialParam(progressbar = TRUE), 
                         return_all = FALSE){
  library(BiocParallel)
  library(scMerge)
  if (is.null(exprs_mat)) {
    stop("exprs_mat is NULL.")
  }
  if (is.null(cell_type)) {
    message("Calculating scSEG index without F-statistics \n")
  }else {
    if (length(cell_type) != ncol(exprs_mat)) {
      stop("length of cell type information is not equal to the number of column of exprs_mat.")
    }
    cell_type <- as.factor(cell_type)
    message("Calculating scSEG index with F-statistics \n")
  }
  z.all <- rowMeans(exprs_mat == 0)
  del <- which(z.all > 0.8)
  if (length(del) > 0) {
    
    if(length(which(z.all<=0.8))>=100){
      message("0.8 zero rate cutoff for gene filtering")
      del <- which(z.all > 0.8)
    }else if(length(which(z.all<=0.9))>=100){
      message("0.9 zero rate cutoff for gene filtering")
      del <- which(z.all > 0.9)
    }else{
      message(paste0(nrow(exprs_mat_filt),' genes left'))
      stop("Not enough gene pass the QC!")
    }
    exprs_mat_filt <- exprs_mat[-del, ]
    message(paste0(nrow(exprs_mat_filt),' genes left'))
  }else {
    exprs_mat_filt <- exprs_mat
  }
  message("Fitting the mixture model... \n")
  paraMat <- scMerge:::make_para_gn_parallel(as.matrix(exprs_mat_filt), 
                                             BPPARAM = BPPARAM)
  r <- paraMat$rho
  s <- paraMat$sigma
  m <- paraMat$mu
  m.scaled <- (m - min(m))/(max(m) - min(m))
  names(r) <- names(s) <- names(m) <- names(m.scaled) <- rownames(exprs_mat_filt)
  genes <- rownames(exprs_mat_filt)
  z <- z.all[genes] * m.scaled[genes]
  if (is.null(cell_type)) {
    x1 = rank(r)/(length(r) + 1)
    x2 = 1 - rank(s)/(length(s) + 1)
    x3 = 1 - rank(z)/(length(z) + 1)
    segIdx <- base::rowMeans(cbind(x1, x2, x3))
    resMat <- data.frame(segIdx = segIdx, rho = r, sigma = s, 
                         mu = m, mu.scaled = m.scaled, zero = z)
  }else {
    message("Fitting ANOVA model for F-stats... \n")
    aovStats <- apply(exprs_mat_filt, 1, function(x) {
      tryCatch(stats::aov(as.numeric(x) ~ cell_type), error = function(e) {
        NULL
      })
    })
    f <- log2(as.numeric(lapply(aovStats, FUN = function(x) {
      tryCatch(summary(x)[[1]]$`F value`[1], error = function(e) {
        NA
      })
    })))
    x1 = rank(r)/(length(r) + 1)
    x2 = 1 - rank(s)/(length(s) + 1)
    x3 = 1 - rank(z)/(length(z) + 1)
    x4 = 1 - rank(f)/(length(f) + 1)
    segIdx <- base::rowMeans(cbind(x1, x2, x3, x4))
    resMat <- data.frame(segIdx = segIdx, rho = r, sigma = s, 
                         mu = m, mu.scaled = m.scaled, zero = z, f_stats = f)
  }
  if (return_all) {
    resMat2 = cbind(gene = as.character(rownames(resMat)), 
                    resMat)
    all_genes = data.frame(gene = as.character(rownames(exprs_mat)), 
                           zero.all = rowMeans(exprs_mat == 0))
    result = merge(x = all_genes, y = resMat2, by = "gene", 
                   all.x = TRUE)
    rownames(result) = rownames(exprs_mat)
  }else {
    result = resMat
  }
  return(result)
}



#raw count as input count=rawcount
run_scMerge<-function(count,cellinfo){
  library(SingleCellExperiment)
  library(scater)
  library(scMerge)
  library(batchelor)
  library(Seurat)
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
  rownames(cellinfo)=cellinfo$Cell
  cellinfo<-cellinfo[colnames(count),]
  cellinfo$Group%<>%factor()
  cellinfo$Batch%<>%factor()
  sob <- CreateSeuratObject(counts = count)
  sob <- NormalizeData(sob)
  sob <- FindVariableFeatures(sob, selection.method = "vst", nfeatures = 2000)
  hvgenes <- head(VariableFeatures(sob), 1000)
  rm(sob)
  
  sce <- SingleCellExperiment(assays = list(counts = count),colData = cellinfo)
  sce <- logNormCounts(sce)
  result = scSEGIndex_used(exprs_mat = count)
  seg <- rownames(count)[which(result$segIdx < .5)]
  scMerge_res <- scMerge(
    sce_combine = sce, 
    ctl = seg,
    kmeansK = rep(2,length(unique(cellinfo$Batch))),
    assay_name = "scMerge_res",
    cell_type = sce$Group,
    replicate_prop = 1, verbose=T,
    marker = hvgenes)
  
  corrected <- scMerge_res@assays@data$scMerge_res 
  save(scMerge_res,corrected,cellinfo, file='./scmerge.rda')
}
