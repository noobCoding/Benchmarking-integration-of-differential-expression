rm(list=ls())

library(scater)
library(SeuratData)
library(Seurat)

base_name='pancreas'
# SeuratData::InstallData("panc8")
panc8 <- panc8.SeuratData::panc8
# split the object by dataset
pancreas.list <- SplitObject(panc8, split.by = "tech")

indrop_cells <- c(which(pancreas.list$indrop@meta.data$celltype=='beta'))

indrop_human1 <- c(which(pancreas.list$indrop$orig.ident == 'human1'))
indrop_cells_human1 <- intersect(indrop_cells, indrop_human1)
length(indrop_cells_human1) #  241 alpha

indrop_human2 <- c(which(pancreas.list$indrop$orig.ident == 'human2'))  
indrop_cells_human2 <- intersect(indrop_cells, indrop_human2)
length(indrop_cells_human2) #  659 alpha

tmp_batch1 = pancreas.list$indrop@assays$RNA@counts[,indrop_cells_human1]
tmp_batch2 = pancreas.list$indrop@assays$RNA@counts[,indrop_cells_human2]

b1 <- CreateSeuratObject(tmp_batch1)
b1 <- SCTransform(b1)

# PCA
b1 <- RunPCA(b1, npcs = 75, verbose = FALSE)
# TSNE
b1 <- RunTSNE(b1, dims = 1:75, seed.use = 7968)

b1 <- FindNeighbors(b1, verbose = FALSE, dims = 1:75)
b1 <- FindClusters(b1, algorithm = 3, random.seed = 7968, resolution = 0.5)
DimPlot(b1, reduction = "tsne", group.by = "seurat_clusters", label = TRUE)


b2 <- CreateSeuratObject(tmp_batch2)
b2 <- SCTransform(b2)
# PCA
b2 <- RunPCA(b2, npcs = 75, verbose = FALSE)
# TSNE
b2 <- RunTSNE(b2, dims = 1:75, seed.use = 7968)

b2 <- FindNeighbors(b2, verbose = FALSE, dims = 1:75)
b2 <- FindClusters(b2, algorithm = 3, random.seed = 7968, resolution = 0.5)
DimPlot(b2, reduction = "tsne", group.by = "seurat_clusters", label = TRUE)

library(dplyr)
# find markers for every cluster compared to all remaining cells, report only the positive ones
b2.markers <- FindAllMarkers(b2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
b2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top10_2 <- b2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(b2, features = top10_2$gene)

# find markers for every cluster compared to all remaining cells, report only the positive ones
b1.markers <- FindAllMarkers(b1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
b1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top10_1 <- b1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(b1, features = top10_1$gene)


b1c0 <- b1.markers$gene[b1.markers$cluster==0 & b1.markers$avg_log2FC > 0.5]
b1c1 <- b1.markers$gene[b1.markers$cluster==1 & b1.markers$avg_log2FC > 0.5]
b1c2 <- b1.markers$gene[b1.markers$cluster==2 & b1.markers$avg_log2FC > 0.5]
b1c3 <- b1.markers$gene[b1.markers$cluster==3 & b1.markers$avg_log2FC > 0.5]

b2c0 <- b2.markers$gene[b2.markers$cluster==0 & b2.markers$avg_log2FC > 0.5]
b2c1 <- b2.markers$gene[b2.markers$cluster==1 & b2.markers$avg_log2FC > 0.5]
b2c2 <- b2.markers$gene[b2.markers$cluster==2 & b2.markers$avg_log2FC > 0.5]

length(intersect(b1c0, b2c0)) # 0
length(intersect(b1c0, b2c1)) # 11

length(intersect(b1c1, b2c0)) # 1
length(intersect(b1c1, b2c1)) # 40

length(intersect(b1c2, b2c0)) # 0
length(intersect(b1c2, b2c1)) # 0

length(intersect(b1c3, b2c0)) # 0
length(intersect(b1c3, b2c1)) # 0

length(which(b1$seurat_clusters == 1))
length(which(b2$seurat_clusters == 1))

# K_intervals = c(2, 3, 4, 6, 7, 8)
K_intervals = 1:10
cutoff = 0.95
a = 2
b = 2
pp = 10

for (K in K_intervals)
{   
  # clustered matching
  batch1 <- tmp_batch1[ ,which(b1$seurat_clusters == 1)]
  batch2 <- tmp_batch2[ ,which(b2$seurat_clusters == 1)]
  batch <- c(rep(1, length(colnames(batch1))), rep(2, length(colnames(batch2))))
  
  # Shuffling replicates    
  ori_cell_order <- colnames(batch1)
  sp_cells <- colnames(batch1)
  sp_cells <- sample(sp_cells, length(sp_cells), replace = FALSE)
  batch1 <- batch1[, sp_cells]
  colnames(batch1) <- ori_cell_order
  
  ori_cell_order <- colnames(batch2)
  sp_cells <- colnames(batch2)
  sp_cells <- sample(sp_cells, length(sp_cells), replace = FALSE)
  batch2 <- batch2[, sp_cells]
  colnames(batch2) <- ori_cell_order
  
  group <- batch
  for (i in unique(group)){
    gid = which(group==i)
    group[gid[1:(length(gid)%/%10*2)]] = paste0(i, "_A")
    group[gid[((length(gid)%/%10*2) + 1):length(gid) ]] = paste0(i, "_B")
  }
  

  
  rowsum = rowSums(as.matrix(batch1)==0)/length(colnames(batch1))
  rowsum2 = rowSums(as.matrix(batch2)==0)/length(colnames(batch2))
  gb1 = rownames(batch1)[which(rowsum < cutoff)]
  print(length(gb1))
  gb2 = rownames(batch2)[which(rowsum2 < cutoff)]
  print(length(gb2))
  common_genes = intersect(gb1, gb2)
  
  # common_genes = sample(common_genes, 5000, replace = F)
  print(length(common_genes))
  
  batch1 = batch1[common_genes, ]  
  batch2 = batch2[common_genes, ]  
  colnames(batch1)<-(group[which(batch==1)])  
  colnames(batch2)<-(group[which(batch==2)])
  
  N = round(length(common_genes) *pp /50)  #
  if (N%%2==1) N = N + 1
  print(N)
  
  DEG_groundthruth = sample(common_genes, N, replace=FALSE)
  second_N_genes  = sample(DEG_groundthruth, N%/%2, replace = FALSE)
  first_N_genes= DEG_groundthruth[!DEG_groundthruth %in% second_N_genes]
  
  #Down sampling
  batch1_ori = batch1
  batch2_ori = batch2

  tmp <- batch1[first_N_genes, which(colnames(batch1)==unique(colnames(batch1))[1])]
  for (idx in 1:length(tmp@x)){
    prob = rbeta(1, a, b)
    tmp@x[idx] = rbinom(1, tmp@x[idx], prob=prob)
  }
  batch1[first_N_genes, which(colnames(batch1)==unique(colnames(batch1))[1])] <- tmp

  tmp <- batch1[second_N_genes, which(colnames(batch1)==unique(colnames(batch1))[2])]
  for (idx in 1:length(tmp@x)){
    prob = rbeta(1, a, b)
    tmp@x[idx] = rbinom(1, tmp@x[idx], prob=prob)
  }
  batch1[second_N_genes, which(colnames(batch1)==unique(colnames(batch1))[2])] <- tmp

  ############################# batch 2
  tmp <- batch2[first_N_genes, which(colnames(batch2)==unique(colnames(batch2))[1])]
  for (idx in 1:length(tmp@x)){
    prob = rbeta(1, a, b)
    tmp@x[idx] = rbinom(1, tmp@x[idx], prob=prob)
  }
  batch2[first_N_genes, which(colnames(batch2)==unique(colnames(batch2))[1])] <- tmp

  tmp <- batch2[second_N_genes, which(colnames(batch2)==unique(colnames(batch2))[2])]
  for (idx in 1:length(tmp@x)){
    prob = rbeta(1, a, b)
    tmp@x[idx] = rbinom(1, tmp@x[idx], prob=prob)
  }
  batch2[second_N_genes, which(colnames(batch2)==unique(colnames(batch2))[2])] <- tmp
  
  ##################################################
  print("bach1")
  mb1=as.matrix(batch1)
  # print(mean(mb1))
  print(mean(mb1[mb1 > 0]))
  print(max(mb1))
  
  colsum = colSums(as.matrix(mb1) == 0)
  simsparsity = sum(colsum)/(length(rownames(mb1))*length(colnames(mb1)))
  print(simsparsity)
  
  print("batch2")
  mb2=as.matrix(batch2)
  # print(mean(mb2))
  print(mean(mb2[mb2 > 0]))
  print(max(mb2))
  
  colsum = colSums(as.matrix(mb2) == 0)
  simsparsity = sum(colsum)/(length(rownames(mb2))*length(colnames(mb2)))
  print(simsparsity)
  
  print("all batches")
  
  ##################################################
  
  batch1 <- batch1[, which(colnames(batch1)=='1_B')]
  colnames(batch1) <- rep('1_B', length(colnames(batch1)))
  first_half <- length(colnames(batch1)) %/% 2
  colnames(batch1)[1:first_half] <- rep('1_A', first_half)
  
  batch2 <- batch2[, which(colnames(batch2)=='2_B')]
  colnames(batch2) <- rep('2_B', length(colnames(batch2)))
  first_half <- length(colnames(batch2)) %/% 2
  colnames(batch2)[1:first_half] <- rep('2_A', first_half)
  
  tbatch <- c(rep(1, length(colnames(batch1))), rep(2, length(colnames(batch2))) )
  ##################################################
  newmat <- cbind.DataFrame(batch1, batch2)
  
  colsum = colSums(as.matrix(newmat) == 0)
  simsparsity = sum(colsum)/(length(rownames(newmat))*length(colnames(newmat)))
  print(simsparsity)
  
  newgroup<-colnames(newmat)
  newgroup[newgroup=='1_A'] = 'A'
  newgroup[newgroup=='1_B'] = 'B'
  newgroup[newgroup=='2_A'] = 'A'
  newgroup[newgroup=='2_B'] = 'B'
  
  newbatch<-tbatch
  geneinfo <- common_genes
  cellinfo <- as.data.frame(newgroup)
  cellinfo$batch <- factor(newbatch)
  
  colnames(cellinfo) <- c("Group", "Batch")
  up_genes <- second_N_genes
  down_genes <- first_N_genes
  de_genes <- DEG_groundthruth
  
  ds = paste0("pan_hm12_alpha_", cutoff*100,"_ab_", a, b,"_", pp,"_", K)
  dir.create(ds)
  write.table(newmat, file = paste0(ds,"/counts.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(geneinfo, file = paste0(ds,"/geneinfo.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(cellinfo, file = paste0(ds,"/cellinfo.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(down_genes, file = paste0(ds,"/true_down_genes.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(up_genes, file = paste0(ds,"/true_up_genes.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(de_genes, file = paste0(ds,"/de_genes.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
}
