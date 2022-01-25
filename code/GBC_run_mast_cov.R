# BiocManager::install('MAST')
library(MAST)
library(Seurat)

base_names <- list.dirs('data')
base_names <-base_names[grepl('epithelial', base_names)]

for(base_name in base_names){
  
  counts<- read.table(file = paste0(base_name,"/counts.txt"), sep = "\t", fill=T)
  cellinfo<-read.table(file = paste0(base_name,"/cellinfo.txt"), sep = "\t", fill=T)
  geneinfo<-read.table(file = paste0(base_name,"/geneinfo.txt"), sep = "\t", fill=T)
  count_df<-counts
  
  myFilteredData <- count_df
  rv_genes<-which(apply(myFilteredData,1,var)==0)
  rv_genes_names<-rownames(myFilteredData)[rv_genes]
  count_df<-myFilteredData[!(rownames(myFilteredData) %in% rv_genes_names),]
  geneinfo <- geneinfo[!(rownames(geneinfo) %in% rv_genes_names),]
  
  cellinfo$Group <- factor(cellinfo$Group)
  cellinfo$Batch <- factor(cellinfo$Batch)
  rownames(cellinfo) <- colnames(counts)
  
  seo <- CreateSeuratObject(counts = count_df, meta.data = cellinfo , project = "MAST", min.cells = 0,)
  seo <- NormalizeData(object = seo, normalization.method = "LogNormalize", scale.factor = 1e4)
  # seo <- SCTransform(seo)
  
  Idents(seo) <- cellinfo$Batch
  Idents(seo) <- cellinfo$Group

  res <- FindMarkers(object = seo, ident.1 = 'nLung', ident.2 = 'tLung',test.use='MAST', logfc.threshold = 0,
                     latent.vars = c('Group', 'Batch'))
  
  write.table(res, file=paste0(base_name,"/all_MAST_Cov_table.txt"), sep = "\t")
}
