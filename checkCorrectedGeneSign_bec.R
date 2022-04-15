library(Seurat)
library(lme4)
library(ggplot2)

dir.create('Corrected_genes')

base_name <- paste0('Corrected_genes/')
vect_HVG <- c('all')
vect_simu <- gsub('data/','',list.dirs('data', recursive=FALSE))
base_name <- paste0('Corrected_genes/')

# METHODS
vect_method <- c('raw', 'Combat', 'limma', 'MNN', 'seurat3', 'scmerge', 'zinbwave')
name_method <- c('Raw', 'Combat', 'limma', 'MNNCorrect', 'Seurat3', 'scMerge', 'ZINB_WaVE')

all_simu <- list()
for (simu in vect_simu){
  
  cellinfo <- read.table(paste0('data/',simu,'/cellinfo.txt'), head=T, fill=T)
  sample_batch12 <- read.table(paste0("data/",simu,"/cellinfo.txt"),sep="\t",header=T,row.names=1, fill=T)
  
  rawdata <- read.table(file = paste0('data/',simu,'/counts.txt'),sep="\t",header=T,row.names=1,check.names = F, fill = T)
  seu <- CreateSeuratObject(rawdata)
  seu <- SCTransform(seu)
  seu <- AddMetaData(seu, cellinfo$Batch, col.name = "Batch")
  seu <- AddMetaData(seu, cellinfo$Group, col.name = "Group")
  Idents(seu) <- seu@meta.data$Group
  raw_markers <- FindMarkers(object = seu, ident.1 = "Group1", ident.2 = "Group2",test.use='wilcox',
                             logfc.threshold = 0, 
                             min.cells.feature = 0,
                             min.cells.group = 0,
                             min.pct = 0,
                             only.pos = F)
  
  all_markers <- list()
  all_markers$raw <- raw_markers
  
  for (mt in vect_method[vect_method != 'raw']){
    tdata <- read.table(file = paste0('demo_',mt,'/',simu,'/','all/output.txt'),sep="\t",header=T,row.names=1,check.names = F, fill=T)  
    tdata[is.na(tdata)] = 0.
    method = name_method[which(vect_method==mt)]
    
    tseu <- CreateSeuratObject(tdata, min.cells = 0, project = base_name, min.features = 0)
    tseu <- AddMetaData(tseu, cellinfo$Batch, col.name = "Batch")
    tseu <- AddMetaData(tseu, cellinfo$Group, col.name = "Group")
    Idents(tseu) <- seu@meta.data$Group
    t_markers <- FindMarkers(object = tseu, ident.1 = "Group1", ident.2 = "Group2",test.use='wilcox',
                             logfc.threshold = 0, 
                             min.cells.feature = 0,
                             min.cells.group = 0,
                             min.pct = 0,
                             only.pos = F)
    all_markers[[mt]] <- t_markers
  }

  
  DEG <- read.table(paste0('data/',simu,'/de__genes.txt'), head=T, fill = T)
  up <- read.table(paste0('data/',simu,'/true_up_genes.txt'), head=T, fill = T)
  down <- read.table(paste0('data/',simu,'/true_down_genes.txt'), head=T, fill = T)
  
  hist_001 <- c()
  perc_001 <- c()
  hist_005 <- c()
  perc_005 <- c()
  hist <- c()
  perc <- c()
  
  for (i in 1:length(vect_method)){
    count_001 = 0
    count_005 = 0
    count = 0
    
    for (g in DEG$x){
      try(
        if (all_markers[[ vect_method[i] ]][g, 'avg_log2FC'] < 0 &&
            g %in% up$x){
          
          if (all_markers[[ vect_method[i] ]][g, 'p_val_adj'] < 0.01){
            count_001 = count_001 + 1
          }
          
          if (all_markers[[ vect_method[i] ]][g, 'p_val_adj'] < 0.05){
            count_005 = count_005 + 1
          }
          
          count = count + 1
        } 
        else if (all_markers[[ vect_method[i] ]][g, 'avg_log2FC'] > 0 &&
                 g %in% down$x){
          
          if (all_markers[[ vect_method[i] ]][g, 'p_val_adj'] < 0.01){
            count_001 = count_001 + 1
          }
          
          if (all_markers[[ vect_method[i] ]][g, 'p_val_adj'] < 0.05){
            count_005 = count_005 + 1
          }
          
          count = count + 1
        }
        , TRUE
      )
    }
    
    hist_001 <- c(hist_001, count_001)
    perc_001 <- c(perc_001, round(100*count_001/length(DEG$x), 2) )
    
    hist_005 <- c(hist_005, count_005)
    perc_005 <- c(perc_005, round(100*count_005/length(DEG$x), 2) )
    
    hist <- c(hist, count)
    perc <- c(perc, round(100*count/length(DEG$x), 2) )
  }
  
  df <- data.frame(percent=perc, percent_005=perc_005, percent_001=perc_001,  
                   count=hist, count_005=hist_005, count_001=hist_001)
  rownames(df) <- vect_method
  all_simu[[simu]] <- df
}

saveRDS(all_simu, file = paste0('4b_bec_deg_swap.rds'))