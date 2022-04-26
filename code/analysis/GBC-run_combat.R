# BiocManager::install('sva')
library(Seurat)
library(sva)

dir.create('demo_Combat')
lsdir <- list.dirs('data', recursive=FALSE) 

sapply(lsdir,function(x){

  x2 <- gsub('data/','',x)
  dir.create(paste0('demo_Combat/',x2), showWarnings = FALSE)
  
  selection <- c('all')
  sapply(selection, function(s){
    
    dir.create(paste0('demo_Combat/',x2,'/',s), showWarnings = FALSE)
    
    # read data counts and cellinfo
    counts <- read.table(paste0(x,'/counts.txt'), head=T, sep='\t')
    cellinfo <- read.table(paste0(x,'/cellinfo.txt'), head=T, sep='\t')
    rownames(cellinfo) <- factor(colnames(counts))
    
    print("Median normalizing counts and log-transforming")
    col_sums = apply(counts,2, sum)
    med_trans = median(col_sums)
    norm_counts = med_trans* scale(counts, center=FALSE, scale=col_sums)
    myFilteredData = log(norm_counts + 1)
    
    # remove genes with variance equals to 0
    rv_genes <- which(apply(myFilteredData, 1, var)==0) # apply on normalized data
    rv_genes_names <- rownames(myFilteredData)[rv_genes]
    count_df <- myFilteredData[!(rownames(myFilteredData) %in% rv_genes_names),]
    
    # Run COMBAT
    t1 = Sys.time()
    combat_output = ComBat(dat=as.matrix(count_df), 
                           batch=cellinfo$Batch, 
                           mod=NULL, par.prior=TRUE, prior.plots=FALSE, mean.only=FALSE)
    t2 = Sys.time()
    
    # save the output
    save(combat_output,file=paste0('demo_Combat/',x2,'/',s,"/output.rda"))
    write.table(combat_output, file = paste0('demo_Combat/',x2,'/',s,"/output.txt"), quote=FALSE, row.names = T, col.names = T, sep="\t")
    
  })
  
})
