# Tutorial: https://bioconductor.org/packages/devel/bioc/vignettes/splatter/inst/doc/splat_params.R
rm(list=ls())

library(splatter)
library(scater)
library(Seurat)

sparsity_level = .95
# dropout = rep(c(0.01, 0.05),times=3) # sparsity 40%
# dropout = rep(c(1.1, 1.3),times=3)     # sparsity 60%
dropout = rep(c(3.7, 3.9),times=3)

b1 <- rep(c(300,300,300),each = 2)
b2 <- rep(c(400,400,400),each = 2)
b3 <- rep(c(400,400,400),each = 2)
b4 <- rep(c(750,750,750),each = 2)
combi <- data.frame(d=dropout,b1=b1,b2=b2,b3=b3, b4=b4)
group_ratio <-  list(c(0.8, 0.2), 
                     c(0.7, 0.3),
                     c(0.6, 0.4),
                     c(0.4, 0.6),
                     c(0.3, 0.7),
                     c(0.2, 0.8))

simulate <- function(nGroups=2, nGenes=5000, dropout=0.5, seed = 17,
                     batchCells = c(300, 400, 400, 750), group.prob = c(0.3, 0.7),
                     de.prob = c(0.1, 0.1),
                     de.downProb = c(1/3, 1/3),
                     de.facLoc = 0.2,
                     de.facScale = 0.2,                     
                     batch.facLoc=0.4,
                     batch.facScale=0.4
){
  
  method <- 'groups'
  dropout.type <- 'experiment'
  sim <- splatSimulate(group.prob=group.prob, nGenes=nGenes, batchCells=batchCells,
                       dropout.type=dropout.type, method=method,seed=seed, dropout.shape=-1, 
                       dropout.mid=dropout, 
                       # lib.loc = lib.loc, lib.scale= lib.scale,
                       de.prob=de.prob, de.downProb=de.downProb,
                       de.facLoc=de.facLoc, de.facScale=de.facScale,
                       batch.facLoc=batch.facLoc, batch.facScale=batch.facScale
  )
  # print(class(sim))
  # sim <- splatSimulateGroups(myParams)
  
  tmpcount = counts(sim)
  colsum = colSums(tmpcount == 0)
  simsparsity = sum(colsum)/(nGenes*(batchCells[1] + batchCells[2] + batchCells[3] + batchCells[4]))
  print("sparsity")
  print(simsparsity)
  
  counts <- as.data.frame(counts(sim))
  truecounts <- as.data.frame(assays(sim)$TrueCounts)
  cellinfo <- as.data.frame(colData(sim))
  geneinfo <- as.data.frame(rowData(sim))
  
  rowsum = rowSums(tmpcount==0)/length(colnames(tmpcount))
  non_sparse_genes = rownames(tmpcount)[which(rowsum < sparsity_level)]
  
  counts <- counts[rownames(counts) %in% non_sparse_genes, ]
  
  rowsum = rowSums(counts == 0)
  simsparsity = sum(rowsum)/(length(rownames(counts)) * length(colnames(counts)) )
  print("sparsity after")
  print(simsparsity)
  
  truecounts <- truecounts[rownames(truecounts) %in% non_sparse_genes,]
  geneinfo <- geneinfo[rownames(geneinfo) %in% non_sparse_genes,]
  
  sim.normalized <- logNormCounts(sim)
  sim.pca <- runPCA(sim.normalized)
  pca <- plotPCA(sim.pca, colour_by = "Group", shape_by="Batch", theme_size=10, size_by="Batch") + 
    # scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'))+
    scale_size_manual(values=c(2, 2, 2, 2)) +
    guides(color = guide_legend(override.aes = list(size=2))) + labs(x = 'PC2', y = 'PC1') + coord_flip()
  
  tmp_int <- sample(1:10000000, size=1, replace = F, prob = NULL)
  ggsave(filename=paste0(tmp_int, '.png'), plot=pca, width = 4, height = 3)
  
  return(list(simu=sim,counts=counts,cellinfo=cellinfo,geneinfo=geneinfo,truecounts=truecounts, id=tmp_int))
}

sapply(1:dim(combi)[1],function(x){
  
  base_name <- paste0('simul',x,'_dropout_',gsub('\\.','',combi$d[[x]]),'_b1_',combi$b1[[x]],'_b2_',combi$b2[[x]],'/')
  dir.create(base_name, showWarnings = FALSE)
  
  sim <- simulate(dropout=combi$d[[x]], batchCells = c(combi$b1[[x]],combi$b2[[x]], combi$b3[[x]], combi$b4[[x]]),
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
})
