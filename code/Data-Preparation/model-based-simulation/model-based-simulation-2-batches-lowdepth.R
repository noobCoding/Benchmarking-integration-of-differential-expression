rm(list=ls())

library(splatter)
library(scater)

sparsity_level = .95
dropout = rep(c(0.001, 0.005),times=3)   # sparsity 80% with lowdepth

b1 <- rep(c(500, 500, 500),each = 2)
b2 <- rep(c(500, 500, 500),each = 2)

combi <- data.frame(d=dropout,b1=b1,b2=b2)
group_ratio <-  list(c(0.8, 0.2), 
                     c(0.7, 0.3),
                     c(0.6, 0.4),
                     c(0.4, 0.6),
                     c(0.3, 0.7),
                     c(0.2, 0.8))

simulate <- function(id=0, nGroups=2, nGenes=5000,seed = 17,  dropout.mid=0.5, 
                     batchCells = c(300, 600), group.prob = c(0.3, 0.7),
                     de.prob = c(0.1, 0.1),
                     de.downProb = c(0.4, 0.4),
                     de.facLoc = 0.15,
                     de.facScale = 0.15,
                     batch.facLoc=0.04,
                     batch.facScale=0.04
){
  
  sim <- splatSimulate(group.prob=group.prob, nGenes=nGenes, batchCells=batchCells,
                       dropout.mid=dropout,  method='groups', 
                       de.prob=de.prob, de.downProb=de.downProb,
                       de.facLoc=de.facLoc, de.facScale=de.facScale,
                       batch.facLoc=batch.facLoc, batch.facScale=batch.facScale,
                       lib.loc=8.2,  lib.scale = 0.5   # avg. 4
                       
  )
  tmpcount = counts(sim)
  colsum = colSums(tmpcount == 0)
  simsparsity = sum(colsum)/(nGenes*(batchCells[1] + batchCells[2]))
  print("sparsity")
  print(simsparsity)
  
  print("depth")
  print(max(tmpcount))
  print(mean(tmpcount[tmpcount>0]))
  
  
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
  s
  colsum = colSums(mb2 == 0)
  simsparsity = sum(colsum)/(length(rownames(mb2))*length(colnames(mb2)))
  print(simsparsity)
  
  non_sparse_genes = intersect(non_sparse_mb1, non_sparse_genes)
  non_sparse_genes = intersect(non_sparse_mb2, non_sparse_genes)
  non_sparse_genes = unique(non_sparse_genes)
  
  truecounts <- truecounts[rownames(truecounts) %in% non_sparse_genes,]
  geneinfo <- geneinfo[rownames(geneinfo) %in% non_sparse_genes,]
  counts <- counts[rownames(counts) %in% non_sparse_genes, ]
  counts[is.na(counts)] = 0
  
  return(list(simu=sim,counts=counts,cellinfo=cellinfo,geneinfo=geneinfo,truecounts=truecounts))
}

id = 0
for (x in 1:dim(combi)[1]){
  id = id + 1
  base_name <- paste0('simul',x,'_dropout_',gsub('\\.','',combi$d[[x]]),'_b1_',combi$b1[[x]],'_b2_',combi$b2[[x]],'/')
  dir.create(base_name, showWarnings = FALSE)
  
  sim <- simulate(id=id, batchCells = c(combi$b1[[x]],combi$b2[[x]]), dropout.mid = combi$d[[x]],
                  group.prob = group_ratio[[x]])
  
  counts <- sim$counts
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
