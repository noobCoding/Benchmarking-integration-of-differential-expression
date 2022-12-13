rm(list=ls())

library(splatter)
library(scater)
library(Seurat)

sparsity_level = .95
# dropout = rep(c(0.01, 0.05),times=3) # sparsity 40%
# dropout = rep(c(3.7, 3.9),times=3)   # sparsity 80%

b1 <- rep(c(300,300,300),each = 2)
b2 <- rep(c(750,750,750),each = 2)
combi <- data.frame(d=dropout,b1=b1,b2=b2)
group_ratio <-  list(c(0.8, 0.2), 
                     c(0.7, 0.3),
                     c(0.6, 0.4),
                     c(0.4, 0.6),
                     c(0.3, 0.7),
                     c(0.2, 0.8))

simulate <- function(id=0, nGroups=2, nGenes=5000, dropout=0.5, seed = 17,
                     batchCells = c(300, 600), group.prob = c(0.3, 0.7),
                     de.prob = c(0.1, 0.1),
                     de.downProb = c(0.4, 0.4),
                     de.facLoc = 0.1,
                     de.facScale = 0.1,
                     batch.facLoc=0.4,
                     batch.facScale=0.4
){
  
  method <- 'groups'
  dropout.type <- 'experiment'
  
  sim <- splatSimulate(group.prob=group.prob, nGenes=nGenes, batchCells=batchCells,method=method,
                       de.prob=de.prob, de.downProb=de.downProb,
                       de.facLoc=de.facLoc, de.facScale=de.facScale,
                       batch.facLoc=batch.facLoc, batch.facScale=batch.facScale
                       ,lib.loc = 10.5,  lib.scale = 0.5
  )

  tmpcount = counts(sim)
  colsum = colSums(tmpcount == 0)
  simsparsity = sum(colsum)/(nGenes*(batchCells[1] + batchCells[2]))
  print("sparsity")
  print(simsparsity)
  
  print("avg. depth")
  print(mean(tmpcount[tmpcount > 0])/3.83)
  
  counts <- as.data.frame(counts(sim))
  truecounts <- as.data.frame(assays(sim)$TrueCounts)
  cellinfo <- as.data.frame(colData(sim))
  geneinfo <- as.data.frame(rowData(sim))
  
  rowsum = rowSums(tmpcount==0)/length(colnames(tmpcount))
  non_sparse_genes = rownames(tmpcount)[which(rowsum < sparsity_level)]
  
  ##################################################
  print("bach1")
  mb1=as.matrix(counts[, cellinfo$Batch=='Batch1'])
  
  rowsum = rowSums(mb1==0)/length(colnames(mb1))
  non_sparse_mb1 = rownames(mb1)[which(rowsum < sparsity_level)]
  
  colsum = colSums(mb1 == 0)
  simsparsity = sum(colsum)/(length(rownames(mb1))*length(colnames(mb1)))
  print(simsparsity)
  
  print("batch2")
  mb2=as.matrix(counts[, cellinfo$Batch=='Batch2'])
  
  rowsum = rowSums(mb2==0)/length(colnames(mb2))
  non_sparse_mb2 = rownames(mb2)[which(rowsum < sparsity_level)]
  
  colsum = colSums(mb2 == 0)
  simsparsity = sum(colsum)/(length(rownames(mb2))*length(colnames(mb2)))
  print(simsparsity)

  # aggerate non-sparse genes  
  non_sparse_genes = intersect(non_sparse_mb1, non_sparse_genes)
  non_sparse_genes = intersect(non_sparse_mb2, non_sparse_genes)
  non_sparse_genes = unique(non_sparse_genes)
  
  truecounts <- truecounts[rownames(truecounts) %in% non_sparse_genes,]
  geneinfo <- geneinfo[rownames(geneinfo) %in% non_sparse_genes,]
  counts <- counts[rownames(counts) %in% non_sparse_genes, ]
  counts[is.na(counts)] = 0
  
  library(Seurat)
  
  seu <- CreateSeuratObject(counts)
  seu <- SCTransform(seu)
  
  # Add cell type annotation to metadata
  seu <- AddMetaData(seu, sim$Batch, col.name = "Batch")
  seu <- AddMetaData(seu, cellinfo$Group, col.name = "Group")
  
  # PCA
  seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
  # TSNE
  seu <- RunTSNE(seu, dims = 1:30, seed.use = 7968)
  # UMAP
  # seu <- RunUMAP(seu, dims = 1:30, seed.use = 7968, uwot.sgd = TRUE)
  
    pca<-DimPlot(seu, reduction = "pca",
          group.by = "Group", shape.by="Batch", pt.size = 1, label = F, repel = T) +
    scale_color_brewer(type = "qual", palette = "Set1") +
    theme(plot.title = element_blank())
  
  tsne<-DimPlot(seu, reduction = "tsne",
          group.by = "Group", shape.by="Batch", pt.size = 1, label = F, repel = TRUE) +
    scale_color_brewer(type = "qual", palette = "Set1")
  
  # DimPlot(seu, reduction = "umap",
  #         group.by = "Group", shape.by="Batch", pt.size = 2,label = F, repel = TRUE) +
  #   scale_color_brewer(type = "qual", palette = "Set1")
  ggsave(filename=paste0(id,'_tsne.png'), plot=tsne, width = 4, height = 3)
  
  return(list(simu=sim,counts=counts,cellinfo=cellinfo,geneinfo=geneinfo))
}

id = 0
for (x in 1:dim(combi)[1]){
  id = id + 1
  base_name <- paste0('simul',x,'_dropout_',gsub('\\.','',combi$d[[x]]),'_b1_',combi$b1[[x]],'_b2_',combi$b2[[x]],'/')
  dir.create(base_name, showWarnings = FALSE)
  
  sim <- simulate(id=id, dropout=combi$d[[x]], batchCells = c(combi$b1[[x]],combi$b2[[x]]),
                group.prob = group_ratio[[x]])
  
  counts <- ceiling(sim$counts/3.83)
  geneinfo <- sim$geneinfo
  cellinfo <- sim$cellinfo
  
  de_genes_ls <- rownames(geneinfo[(geneinfo$DEFacGroup1+geneinfo$DEFacGroup2)!=2,])
  de_genes_df <- geneinfo[de_genes_ls,]
  down_genes <- rownames(de_genes_df[de_genes_df$DEFacGroup1<de_genes_df$DEFacGroup2,])
  up_genes <- de_genes_ls[!de_genes_ls %in% down_genes]
  genes_deg <- c(up_genes, down_genes)
  
  write.table(counts, file = paste0(base_name,"/counts.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(geneinfo, file = paste0(base_name,"/geneinfo.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(cellinfo, file = paste0(base_name,"/cellinfo.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(down_genes, file = paste0(base_name,"/true_down_genes.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(up_genes, file = paste0(base_name,"/true_up_genes.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(de_genes_ls, file = paste0(base_name,"/de__genes.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
}
