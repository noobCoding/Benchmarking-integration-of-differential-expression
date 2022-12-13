rm(list=ls())
library(Seurat)

# METHODS
vect_method <- c('combat', 'limma_bec', 'mnn_opt', 'Seurat', 'scMerge',  'zinbwave', 'scvi', 'scgen', 'scanorama', 'risc')
name_method <- c('Combat', 'limma_BEC', 'MNN', 'Seurat_BEC', 'scMerge',  'ZINBWaVE_BEC', 'scVI', 'scGen', 'Scanorama', 'RISC')

vect_simu <- gsub('data/','',list.dirs('data', recursive=FALSE))
# vect_simu <- vect_simu[c(1, 2, 3, 4, 5, 6)]
vect_simu <- vect_simu[3]

all_logfc <- list()
mcolor <- c(
  "#009E73", "gold1", "#CC79A7", "dodgerblue2", "#FF7F00", "#FB9A99", "#0072B2", "green4", "#E31A1C", "#6A3D9A",
  "skyblue2", "#FDBF6F", 
  "palegreen2",
  "maroon", "orchid1", "deeppink1",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "#E69F00",  "#56B4E9"
)

# Set space for 2 rows and 5 columns
par(mfrow=c(2,5))

for (simu in vect_simu){
  # import data, sample
  counts <- read.table(file = paste0('data/',simu, '/', 'counts.txt'),header=T,row.names=1,check.names = F)
  cellinfo <- read.table(paste0("data/",simu,"/cellinfo.txt"),sep="\t",header=T,row.names=1)
  colnames(counts) <- rownames(cellinfo)
  de_genes <- read.table(paste0('data/',simu,'/de__genes.txt'), head=T, fill = T)
  de_genes <- de_genes$x
  
  #### 
  batch1 <- rownames(cellinfo[cellinfo$Batch=="Batch1",])
  batch2 <- rownames(cellinfo[cellinfo$Batch=="Batch2",])

  group1 <- rownames(cellinfo[cellinfo$Group=="Group1",])
  group2 <- rownames(cellinfo[cellinfo$Group=="Group2",])

  # batch1 <- rownames(cellinfo[cellinfo$Batch=="1",])
  # batch2 <- rownames(cellinfo[cellinfo$Batch=="2",])
  # group1 <- rownames(cellinfo[cellinfo$Group=="A",])
  # group2 <- rownames(cellinfo[cellinfo$Group=="B",])

  #batch1 2
  cog_1 <- counts[, group1]
  cog_2 <- counts[, group2]

  lib_1 <- colSums(cog_1)
  lib_2 <- colSums(cog_2)

  genes1 <- c()
  for (g in de_genes){
    tg_1 <- cog_1[rownames(cog_1)==g, ]
    tg_1 <- tg_1/lib_1
    tg_1 <- tg_1[tg_1 != 0]
    tg_1 <- log2(1 +  tg_1 * 10^4)
    g1 <- sum(tg_1)/length(tg_1)

    tg_2 <- cog_2[rownames(cog_2)==g, ]
    tg_2 <- tg_2/lib_2
    tg_2 <- tg_2[tg_2 != 0]
    tg_2 <- log2(1 +  tg_2 * 10^4)
    g2 <- sum(tg_2)/length(tg_2)

    genes1 <- c(genes1, g1 - g2)
  }
  names(genes1) <- de_genes
  logFC_raw <- genes1
  all_logfc$raw <- logFC_raw
  
  df <- data.frame()
  for (method in vect_method){
      # import data, sample
      if (method %in% c('scvi',  'scgen','scanorama')){
        corrected_data <- read.csv(file = paste0('data/',simu, '/', method , '_corrected_data.csv'),header=T,row.names=1,check.names = F)
        if (method=='scanorama'){
          corrected_data <- corrected_data * 10
        }
      } else if (method == 'risc') {
        load(paste0('data/',simu,'/first_processed/', 'splatter.0.05mincell_filtered_RISC_first_processed.RData'))
        corrected_data <- first_processed$count
        
      } else if (method == 'zinbwave') {
        corrected_data <- readRDS(paste0('data/',simu, '/zinbwave/splatter.0.05mincell_filtered_LC_output.rds'))
        corrected_data <- corrected_data@assays@data$normalizedValues
        
      } else if (method == 'mnn_opt') {
        load(paste0('data/',simu,'/first_processed/', 'splatter.0.05mincell_filtered_', method,'_first_processed.RData'))
        corrected_data <- first_processed$count* 20
        
      } else {
        load(paste0('data/',simu,'/first_processed/', 'splatter.0.05mincell_filtered_', method,'_first_processed.RData'))
        corrected_data <- first_processed$count
      }
    
     
      ####
     
      group1 <- which(first_processed$group=='Group1')
      group2 <- which(first_processed$group=='Group2')
   
      cog_1 <- corrected_data[, group1]
      cog_2 <- corrected_data[, group2]
      
      lib_1 <- colSums(cog_1)
      lib_2 <- colSums(cog_2)

      logFC_corrected <- c()
      if (method %in% c('risc', 'scvi')){
        for (g in de_genes){
          tg_1 <- cog_1[rownames(cog_1)==g, ]
          tg_1 <- tg_1/lib_1
          tg_1 <- tg_1[tg_1 != 0]
          tg_1 <- log2(1 +  tg_1 * 10^4)
          g1 <- sum(tg_1)/length(tg_1)

          tg_2 <- cog_2[rownames(cog_2)==g, ]
          tg_2 <- tg_2/lib_2
          tg_2 <- tg_2[tg_2 != 0]
          tg_2 <- log2(1 +  tg_2 * 10^4)
          g2 <- sum(tg_2)/length(tg_2)

          logFC_corrected <- c(logFC_corrected, g1 - g2)
        }
      }
      else {
        for (g in de_genes){
          tg_1 <- as.numeric(cog_1[rownames(cog_1)==g, ])
          g1 <- log(x=mean(x=expm1(x=tg_1)) + 1)

          tg_2 <- as.numeric(cog_2[rownames(cog_2)==g, ])
          g2 <- log(mean(expm1(tg_2)) + 1)

          logFC_corrected <- c(logFC_corrected, g1 - g2)
        }
      }
      names(logFC_corrected) <- de_genes

      logFC_corrected[is.na(logFC_corrected)] <- 0
      all_logfc[[method]] <- logFC_corrected

      res<-cor.test(logFC_raw, logFC_corrected, method='pearson')
      print(method)
      print(res$estimate)
      print(res$p.value)

      cos_dist = 0
      for (i in 1:length(logFC_raw)){
        cos_dist = cos_dist + 1 - sign(logFC_raw[i])* c(logFC_raw[i], logFC_corrected[i]) %*% c(1,1) / 
                                                    ( sqrt(2)*sqrt(logFC_raw[i]^2 + logFC_corrected[i]^2) )
      }
      cos_dist = cos_dist/length(logFC_raw)
      
      ### plotting
      plot(logFC_raw, logFC_corrected, col=mcolor[which(vect_method==method)],  pch = 16, 
           xlab = "logFC_raw", ylim=c(-1,1), xlim=c(-1,1),cex.lab=1.3, cex.axis=1.3, cex.main=2., cex.sub=1.3,
           cex=1.5, main=name_method[which(vect_method==method)], axes=FALSE)
      axis(1, pos=0)
      axis(2, pos=0)
      
      text(-.5, .8, paste0("Pearson: ", round(res$estimate, 3)), cex = 1.3, col='blue')
      text(-.5, .6, paste0("p-value: ", formatC(res$p.value, format='e', digits=2)), cex = 1.3, col='blue')
      text(-.5, .4, paste0("Angular Dist.: ", round(cos_dist,3)), cex = 1.3, col='blue')
      
  }
}



