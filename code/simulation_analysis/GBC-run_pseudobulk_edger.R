# BiocManager::install('edgeR')
# BiocManager::install("apeglm")
# BiocManager::install("DESeq2")
# BiocManager::install("zinbwave")
# BiocManager::install("scRNAseq")
# BiocManager::install("matrixStats")
# BiocManager::install("magrittr")
# BiocManager::install("biomaRt")
library(edgeR)
library(zinbwave)
# library(scRNAseq)
# library(matrixStats)
library(magrittr)
# library(ggplot2)
library(biomaRt)
library(BiocParallel)

BiocParallel::register(BiocParallel::SerialParam())

base_names <- list.dirs('data')
base_names <-base_names[grepl('simu', base_names)]
# base_name <- base_names[3]

for(base_name in base_names){

  counts<- read.table(file = paste0(base_name,"/counts.txt"), sep = "\t", fill=T)
  cellinfo<-read.table(file = paste0(base_name,"/cellinfo.txt"), sep = "\t", fill=T)
  geneinfo<-read.table(file = paste0(base_name,"/geneinfo.txt"), sep = "\t", fill=T)
  count_df<-counts
  # count_df %<>% extract(!is.na(.))
  count_df[is.na(count_df)] = 0.

  nGroup <- length(unique(cellinfo$Group))
  nBatch <- length(unique(cellinfo$Batch))
  collapsed_counts <- matrix(1:nBatch*nGroup*nrow(count_df), nrow=nrow(count_df), ncol = nBatch*nGroup) * 0

  tagBatch <- unique(cellinfo$Batch)
  tagGroup <- unique(cellinfo$Group)

  shuffle = F
  # # Shuffling replicates
  if (shuffle){
    tagBatch[tagGroup=='Group1'] <- sample(tagBatch[tagGroup=='Group1'], 
                                           length(tagBatch[tagGroup=='Group1']), 
                                           replace = F)
    tagBatch[tagGroup=='Group2'] <- sample(tagBatch[tagGroup=='Group2'], 
                                           length(tagBatch[tagGroup=='Group2']), 
                                           replace = F)
  }

  for (i in 1:nBatch){
    for (j in 1:nGroup){
      iBatch <- which(cellinfo$Batch == tagBatch[i])
      iGroup <- which(cellinfo$Group == tagGroup[j])
      idx <- intersect(iBatch, iGroup)
      tmp <- count_df[, idx]
      collapsed_counts[, (i-1)*nGroup + j] <- rowSums(tmp)
    }
  }
  collapsed_counts <- data.frame(collapsed_counts)
  rownames(collapsed_counts) <- rownames(count_df)

  Batch <- c()
  for (i in 1:nBatch) {Batch<-c(Batch, rep(tagBatch[i], nGroup))} 
  Group <- rep(c(tagGroup[1], tagGroup[2]), nBatch*nGroup%/%2) 
  

  collapsed_cellinfo <- data.frame(Batch=Batch, Group=Group)
  collapsed_cellinfo$Group <- factor(collapsed_cellinfo$Group)
  collapsed_cellinfo$Batch <- factor(collapsed_cellinfo$Batch)

  y <- DGEList(counts=collapsed_counts, group=collapsed_cellinfo$Group)
  y <- calcNormFactors(y)

  design <- model.matrix(~Group, data=collapsed_cellinfo)
  y <- estimateDisp(y, design)
  
  # #LRT
  # fit <- glmFit(y, design)
  # test <- glmLRT(fit)

  # # QLF: default
  fit = glmQLFit(y, design)
  test = glmQLFTest(fit, coef = -1)

  test <- topTags(test, n = Inf)
  result.table <- data.frame('pvalue' = test$table$PValue, 
                             'adjpvalue' = test$table$FDR, 
                             'logFC' = test$table$logFC)
  rownames(result.table) <- rownames(test)
  
  write.table(result.table, file=paste0(base_name,"/all_pseudobulk_edger_result.txt"), sep = "\t")
}

