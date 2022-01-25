seurat_analysis_deg <- function(TPM=TPM,
                                sample=sample,
                                group_col="Site",
                                base_name="GC_cluster7", 
                                mt="^MT-", 
                                scale=F, 
                                scale_factor=1e4,
                                group1=NULL,
                                group2=NULL,
                                test.use="wilcox",
                                col.low = "#FF00FF",
                                col.mid = "#000000", col.high = "#FFFF00",
                                logfc.threshold = 0.0,
                                plotDEGhm=T,topn=NULL){
  
  library(Seurat)
  library(dplyr)
  library(Matrix)
  
  pbmc.data <- TPM
  
  common_cell <- intersect(row.names(sample),colnames(pbmc.data))
  pbmc.data <- pbmc.data[,common_cell]
  if(ncol(sample)==1) sample <- cbind(sample,sample)
  sample <- sample[common_cell,]
  
  if((!is.null(group1)) && (!is.null(group2))){
    group12 <- row.names(sample)[as.character(sample[,group_col])==as.character(group1) | as.character(sample[,group_col])==as.character(group2)]
    pbmc.data <- pbmc.data[,group12]
  }
  
  if(!is.null(sample)){
    colnames(pbmc.data) <- paste(sample[colnames(pbmc.data),group_col],colnames(pbmc.data),sep="_")
  }
  
  pbmc <- CreateSeuratObject(counts = pbmc.data, min.cells = 0, project = base_name, min.features = 0)
  
  pbmc.markers <- FindMarkers(object = pbmc, ident.1 = group1, ident.2 = group2,test.use=test.use,
                              logfc.threshold = logfc.threshold, 
                              min.cells.feature = 0,
                              min.cells.group = 0,
                              min.pct = 0,
                              only.pos = F)

  ####################
  
  # deg <- pbmc.markers[pbmc.markers$p_val_adj<cutoff,]
  deg <- pbmc.markers
  deg <- deg[order(deg$avg_log2FC,decreasing = T),]
  genelist <- rownames(deg)
  write.table(deg,paste(base_name,"_seurat_auc_DEG.txt",sep=""),sep="\t",col.names=NA, quote=F)
}


seurat_analysis_deg2 <- function(TPM=TPM,
                                 sample=sample,
                                 group_col="Site",
                                 base_name="GC_cluster7", 
                                 mt="^MT-", 
                                 scale=F, 
                                 scale_factor=1e4,
                                 group1=NULL,
                                 group2=NULL,
                                 test.use="wilcox",
                                 col.low = "#FF00FF",
                                 col.mid = "#000000", col.high = "#FFFF00",
                                 logfc.threshold = 0.0,
                                 plotDEGhm=T,topn=NULL){
  
  library(Seurat)
  library(dplyr)
  library(Matrix)
  
  pbmc.data <- TPM
  
  common_cell <- intersect(row.names(sample),colnames(pbmc.data))
  pbmc.data <- pbmc.data[,common_cell]
  
  if(ncol(sample)==1) sample <- cbind(sample,sample)
  sample <- sample[common_cell,]
  
  if((!is.null(group1)) & (!is.null(group2))){
    group12 <- row.names(sample)[as.character(sample[,group_col])==as.character(group1) | as.character(sample[,group_col])==as.character(group2)]
    pbmc.data <- pbmc.data[,group12]
  }
  
  if(!is.null(sample)){
    colnames(pbmc.data) <- paste(sample[colnames(pbmc.data),group_col],colnames(pbmc.data),sep="_")
  }
  
  pbmc <- CreateSeuratObject(counts = pbmc.data, min.cells = 0, project = base_name)
  
  pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = scale_factor)
  
  pbmc.markers <- FindMarkers(object = pbmc, ident.1 = group1, ident.2 = group2,test.use=test.use,
                              logfc.threshold = logfc.threshold, 
                              min.cells.feature = 0,
                              min.cells.group = 0,
                              min.pct = 0,
                              only.pos = F)
  
  ####################  
  deg <- pbmc.markers
  deg <- deg[order(deg$avg_log2FC,decreasing = T),]
  genelist <- rownames(deg)
  
  write.table(deg, paste(base_name,"_seurat_auc_DEG.txt",sep=""),sep="\t",col.names=NA, quote=F)
}
