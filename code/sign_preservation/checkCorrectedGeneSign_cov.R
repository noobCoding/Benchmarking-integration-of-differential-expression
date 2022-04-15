########## summary confusion matrix
# rm(list=ls())
vect_simu <- gsub('data/','',list.dirs('data', recursive=FALSE))

# load libraries
cov_method <- c( "mast",  "mast_cov",  "deseq2",  "deseq2_cov",  "zinbwave_deseq2",  "zinbwave_deseq2_cov",
                 "edger_detrate",  "edger_detrate_cov",  "edger",  "edger_cov",  "zinbwave_edger",  "zinbwave_edger_cov",
                 "limma_voom",  "limma_voom_cov",  "limma_trend_false",  "limma_trend_false_cov",
                 "combat_limma_trend_false",  "mnn_limma_trend_false",  "scmerge_limma_trend_false")
all_simu <- list()


for (simu in vect_simu){
  DEG <- read.table(paste0('data/',simu,'/de__genes.txt'), head=T, fill = T)
  up <- read.table(paste0('data/',simu,'/true_up_genes.txt'), head=T, fill = T)
  down <- read.table(paste0('data/',simu,'/true_down_genes.txt'), head=T, fill = T)
  
  hist_001 <- c()
  perc_001 <- c()
  hist_005 <- c()
  perc_005 <- c()
  hist <- c()
  perc <- c()
  
  for (method in cov_method){
    print(method)
    if (method=='deseq2'){
      markers <-read.table(paste0('data/',simu,'/','all_deseq2_result_table.txt'), head=T,  sep='\t', fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
    } else if (method=='deseq2_cov'){
      markers <-read.table(paste0('data/',simu,'/','all_deseq2_cov_table.txt'), head=T,  sep='\t', fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
      
    } else if (method=='edger'){
      markers <-read.table(paste0('data/',simu,'/','all_edger_result_table.txt'), head=T,  sep='\t', fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
      
    } else if (method=='edger_cov'){
      markers <-read.table(paste0('data/',simu,'/','all_edger_cov_table.txt'), head=T,  sep='\t', fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
      
    } else if (method=='edger_detrate'){
      markers <-read.table(paste0('data/',simu,'/','all_edger_detrate_table.txt'), head=T,  sep='\t', fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
      
    }else if (method=='edger_detrate_cov'){
      markers <-read.table(paste0('data/',simu,'/','all_edger_detrate_cov_table.txt'), head=T,  sep='\t', fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
      
    }  else if (method=='limma_voom'){
      markers <-read.table(paste0('data/',simu,'/','voom_result_table.txt'), head=T,  sep='\t', fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
      
    } else if (method=='limma_voom_cov'){
      markers <-read.table(paste0('data/',simu,'/','voom_cov_table.txt'), head=T,  sep='\t', fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
      
    } else if (method=='limma_trend'){
      markers <-read.table(paste0('data/',simu,'/','limma_trend_result.txt'), head=T,  sep='\t', fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
      
    } else if (method=='limma_trend_cov'){
      markers <-read.table(paste0('data/',simu,'/','limma_trend_cov_result.txt'), head=T,  sep='\t', fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
      
    } else if (method=='limma_trend_false'){
      markers <-read.table(paste0('data/',simu,'/','limma_trend_result_false.txt'), head=T,  sep='\t', fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
      
    } else if (method=='limma_trend_false_cov'){
      markers <-read.table(paste0('data/',simu,'/','limma_trend_cov_false.txt'), head=T,  sep='\t', fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
      
    } else if (method=='mnn_limma_trend'){
      markers <-read.table(paste0('data/',simu,'/','limma_trend_MNN_result.txt'), head=T,  sep='\t', fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
      
    } else if (method=='mnn_limma_trend_false'){
      markers <-read.table(paste0('data/',simu,'/','limma_trend_MNN_result_false.txt'), head=T,  sep='\t', fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
      
    } else if (method=='scmerge_limma_trend'){
      markers <-read.table(paste0('data/',simu,'/','limma_trend_scMerge_result.txt'), head=T,  sep='\t', fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
      
    } else if (method=='scmerge_limma_trend_false'){
      markers <-read.table(paste0('data/',simu,'/','limma_trend_scMerge_result_false.txt'), head=T,  sep='\t', fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
      
    } else if (method=='combat_limma_trend'){
      markers <-read.table(paste0('data/',simu,'/','limma_trend_combat_result.txt'), head=T,  sep='\t', fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
      
    } else if (method=='combat_limma_trend_false'){
      markers <-read.table(paste0('data/',simu,'/','limma_trend_combat_result_false.txt'), head=T,  sep='\t', fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
      
    } else if (method=='zinbwave_deseq2'){
      markers <-read.table(paste0('data/',simu,'/','all_zinbwave_deseq2_table.txt'), head=T,  sep='\t', fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
      
    } else if (method=='zinbwave_deseq2_cov'){
      markers <-read.table(paste0('data/',simu,'/','all_zinbwave_deseq2_table_cov.txt'), head=T,  sep='\t', fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
      
    } else if (method=='zinbwave_edger'){
      markers <-read.table(paste0('data/',simu,'/','all_zinbwave_edger_table.txt'), head=T,  sep='\t', fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
      
    } else if (method=='zinbwave_edger_cov'){
      markers <-read.table(paste0('data/',simu,'/','all_zinbwave_edger_table_cov.txt'), head=T,  sep='\t', fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
      
    } 
    else if (method=='mast'){
      markers <-read.table(paste0('data/',simu,'/','all_MAST_result_table.txt'), head=T, fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
      
      markers <- markers[, c(1, 2, 6, 3)]
      colnames(markers) <- c("genes", "pvalue",	"adjpvalue",	"logFC")
      markers$logFC <- -1. * markers$logFC
      
    } else if (method=='mast_cov'){
      markers <-read.table(paste0('data/',simu,'/','all_MAST_Cov_table.txt'), head=T, fill=T)
      markers <- data.frame(genes = rownames(markers), markers)
      
      markers <- markers[, c(1, 2, 6, 3)]
      colnames(markers) <- c("genes", "pvalue",	"adjpvalue",	"logFC")
      markers$logFC <- -1. * markers$logFC
    }
    
  
    count_001 = 0
    count_005 = 0
    count = 0
    
    for (g in DEG$x){
      try(
        if (markers[g, 'logFC'] > 0 &&
            g %in% up$x){
          
          if (markers[g, 'adjpvalue'] < 0.01){
            count_001 = count_001 + 1
          }
          
          if (markers[g, 'adjpvalue'] < 0.05){
            count_005 = count_005 + 1
          }
          
          count = count + 1
        } 
        else if (markers[g, 'logFC'] < 0 &&
                 g %in% down$x){
          
          if (markers[g, 'adjpvalue'] < 0.01){
            count_001 = count_001 + 1
          }
          
          if (markers[g, 'adjpvalue'] < 0.05){
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
  rownames(df) <- cov_method
  all_simu[[simu]] <- df
}

  

saveRDS(all_simu, file = paste0('4b_cov_deg_swap.rds'))