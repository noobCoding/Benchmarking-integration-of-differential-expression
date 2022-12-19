rm(list=ls())

library(splatter)
library(scater)

sparsity_level = .95
# dropout = rep(c(0.01, 0.05),times=3) # sparsity 40%
dropout = rep(c(3.7, 3.9),times=3)   # sparsity 80%

s <- list()
gid <- c()
for (i in 1:6){
  s <- rbind(s, sample(c(100, 150, 200, 300, 400, 500, 750), 7))
  gid <- c(gid, sample(1:6, 1, replace=TRUE))
}
print(s)
print(gid)

combi <- data.frame(as.matrix(s))
combi$d <- dropout
combi$gid <- gid

group_ratio <-  list(c(0.8, 0.2),
                     c(0.7, 0.3),
                     c(0.6, 0.4),
                     c(0.4, 0.6),
                     c(0.3, 0.7),
                     c(0.2, 0.8))

simulate <- function(id=0, gid=0, nGroups=2, nGenes=5000, dropout=0.5, seed = 17,
                     batchCells = c(100, 150, 200, 300, 400, 500, 750), group.prob = c(0.3, 0.7),
                     de.prob = c(0.1, 0.1),
                     de.downProb = c(0.4, 0.4),
                     de.facLoc = c(0.15, 0.15),
                     de.facScale = c(0.15, 0.15),
                     batch.facLoc=0.2,
                     batch.facScale=0.2
){
  
  method <- 'groups'
  dropout.type <- 'experiment'
  
  sim <- splatSimulate(group.prob=group.prob, nGenes=nGenes, batchCells=batchCells, group.prob=group_ratio[gid],
                       dropout.type=dropout.type, method=method,seed=seed, dropout.shape=-1, 
                       dropout.mid=dropout, 
                       de.prob=de.prob, de.downProb=de.downProb,
                       de.facLoc=de.facLoc, de.facScale=de.facScale,
                       batch.facLoc=batch.facLoc, batch.facScale=batch.facScale,
  )
  
  tmpcount = counts(sim)
  colsum = colSums(tmpcount == 0)
  simsparsity = sum(colsum)/(nGenes*(batchCells[1] + batchCells[2]  + batchCells[3]  + batchCells[4]
                                     + batchCells[5]  + batchCells[6]  + batchCells[7]))
  print("sparsity")
  print(simsparsity)
  
  counts <- as.data.frame(counts(sim))
  truecounts <- as.data.frame(assays(sim)$TrueCounts)
  cellinfo <- as.data.frame(colData(sim))
  geneinfo <- as.data.frame(rowData(sim))
  
  rowsum = rowSums(tmpcount==0)/length(colnames(tmpcount))
  non_sparse_genes = rownames(tmpcount)[which(rowsum < sparsity_level)]
  
  counts <- counts[rownames(counts) %in% non_sparse_genes, ]
  counts[is.na(counts)] = 0
  
  ##################################################
  mb1=as.matrix(counts[, cellinfo$Batch=='Batch1'])
  rowsum = rowSums(mb1==0)/length(colnames(mb1))
  non_sparse_mb1 = rownames(mb1)[which(rowsum < sparsity_level)]
  
  print("batch2")
  mb2=as.matrix(counts[, cellinfo$Batch=='Batch2'])
  rowsum = rowSums(mb2==0)/length(colnames(mb2))
  non_sparse_mb2 = rownames(mb2)[which(rowsum < sparsity_level)]
  
  print("batch3")
  mb3=as.matrix(counts[, cellinfo$Batch=='Batch3'])
  rowsum = rowSums(mb3==0)/length(colnames(mb3))
  non_sparse_mb3 = rownames(mb3)[which(rowsum < sparsity_level)]
  
  print("batch4")
  mb4=as.matrix(counts[, cellinfo$Batch=='Batch4'])
  rowsum = rowSums(mb4==0)/length(colnames(mb4))
  non_sparse_mb4 = rownames(mb4)[which(rowsum < sparsity_level)]
  
  print("batch5")
  mb5=as.matrix(counts[, cellinfo$Batch=='Batch5'])
  rowsum = rowSums(mb5==0)/length(colnames(mb5))
  non_sparse_mb5 = rownames(mb5)[which(rowsum < sparsity_level)]
  
  print("batch6")
  mb6=as.matrix(counts[, cellinfo$Batch=='Batch6'])
  rowsum = rowSums(mb6==0)/length(colnames(mb6))
  non_sparse_mb6 = rownames(mb6)[which(rowsum < sparsity_level)]
  
  print("batch7")
  mb7=as.matrix(counts[, cellinfo$Batch=='Batch7'])
  rowsum = rowSums(mb7==0)/length(colnames(mb7))
  non_sparse_mb7 = rownames(mb3)[which(rowsum < sparsity_level)]

  # aggerate non-sparse genes  
  non_sparse_genes = intersect(non_sparse_mb1, non_sparse_genes)
  non_sparse_genes = intersect(non_sparse_mb2, non_sparse_genes)
  non_sparse_genes = intersect(non_sparse_mb3, non_sparse_genes)
  non_sparse_genes = intersect(non_sparse_mb4, non_sparse_genes)
  non_sparse_genes = intersect(non_sparse_mb5, non_sparse_genes)
  non_sparse_genes = intersect(non_sparse_mb6, non_sparse_genes)
  non_sparse_genes = intersect(non_sparse_mb7, non_sparse_genes)
  non_sparse_genes = unique(non_sparse_genes)
  
  counts <- counts[rownames(counts) %in% non_sparse_genes,]
  geneinfo <- geneinfo[rownames(geneinfo) %in% non_sparse_genes,]
  
  return(list(simu=sim,counts=counts,cellinfo=cellinfo,geneinfo=geneinfo))
}

id = 0
for (x in 1:dim(combi)[1]){
  id = id + 1
  base_name <- paste0('simul',x,'_dropout_',gsub('\\.','',combi$d[[x]]),'_gid_',combi$gid[[x]])
  dir.create(base_name, showWarnings = FALSE)
  
  cmb<-unlist(t(combi)[,x])[1:7]
  sim <- simulate(id=id, dropout=combi$d[[x]], gid=combi$gid[[x]], batchCells = cmb)
  
  counts <- sim$counts
  geneinfo <- sim$geneinfo
  cellinfo <- sim$cellinfo
  print(nrow(geneinfo))
  print(nrow(counts))
  
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