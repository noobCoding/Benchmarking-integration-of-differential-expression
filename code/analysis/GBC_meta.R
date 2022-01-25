########## summary confusion matrix

rm(list=ls())

# load libraries
library(gridExtra)
library(openxlsx)
dir.create('Venn_diagram')
meta_method <- c(  "voom+REM"                                          ,"voom+FEM"                                         ,
                   "voom+modt_wFisher_sample+gene"                     ,  "voom+modt_Fisher"                                  ,
                   "edgeR_wFisher_sample+gene"                        ,
                   "deseq2_wFisher_sample+gene"                        ,
                   "LogNormalize+REM"                                 ,
                   "LogNormalize+FEM"                                  ,
                   "LogNormalize+limma.trend_wFisher_sample+gene"     ,
                   "voom+modt_Fisher"              ,  
                   "edgeR_Fisher"                   , 
                   "deseq2_Fisher"                   ,
                   "LogNormalize+limma.trend_Fisher" ,
                   "deseq2+FEM"                                        , "deseq2+REM"                
)

df_med <- c(      "voom+REM", "voom+FEM",
                  "LogNormalize+REM", "LogNormalize+FEM",
                  "deseq2+REM", "deseq2+FEM"
)

meta_med <- c(    "voom+modt_Fisher"              ,  
                  "edgeR_Fisher"                   , 
                  "deseq2_Fisher"                   ,
                  "LogNormalize+limma.trend_Fisher"
)

newmeta_med <- c( "voom+modt_wFisher_sample+gene"                     ,
                  "edgeR_wFisher_sample+gene"                        , 
                  "deseq2_wFisher_sample+gene"                        ,
                  "LogNormalize+limma.trend_wFisher_sample+gene"     
)


vect_method <- c('Combat', 'limma_bec',
                 'MNN', 'mnn_limma_trend_false',
                 'scmerge', 'scmerge_limma_trend_false',
                 'seurat3', 'zinbwave', 'zinbwave_deseq2', 'zinbwave_edger', 
                 'zinbwave_deseq2_cov', 'zinbwave_edger_cov',
                 'raw_data', 'mast', 'mast_cov', 
                 'deseq2', 'deseq2_cov', 
                 'edger', 'edger_cov', 'edger_detrate', 'edger_detrate_cov',
                 "limma_voom", "limma_voom_cov",
                 "limma_trend_false", "limma_trend_false_cov",
                 "combat_limma_trend_false",
                 
                 "voom+REM"                                          ,"voom+FEM"                                         ,
                 "LogNormalize+REM"                                  ,"LogNormalize+FEM"                                  ,
                 "deseq2+REM"                                        ,"deseq2+FEM"                ,
                 "voom+modt_Fisher"                                  ,"voom+modt_wFisher_sample+gene"                     ,
                 "edgeR_Fisher"                                      ,"edgeR_wFisher_sample+gene"                        ,
                 "deseq2_Fisher"                                     ,"deseq2_wFisher_sample+gene"                        ,
                 "LogNormalize+limma.trend_Fisher"                   ,"LogNormalize+limma.trend_wFisher_sample+gene"      
)

main_Fscore <- function(select){
  vect_HVG <- c('all')
  
  # SIMULATIONS
  vect_simu <- gsub('data/','',list.dirs('data', recursive=FALSE))
  
  Fscore_list <- lapply(vect_HVG,function(HVG){
    base_name <- paste0('Venn_diagram/')
    
    df_all <- lapply(vect_simu,function(simu){
      
      df <- sapply(vect_method,function(method){
        
        if(method=='raw_data'){
          S3 <- read.table(paste0('Seurat_DEGs_auc/',simu,'/','S3_batch12/',method,'_',HVG,'/degs_batch12__seurat_auc_DEG.txt'), head=T, sep='\t', fill=T)
          S3 <- data.frame(genes = rownames(S3), S3)
        } else {
          if (method=='deseq2'){
            S3 <-read.table(paste0('data/',simu,'/','all_deseq2_result_table.txt'), head=T,  sep='\t', fill=T)
            S3 <- data.frame(genes = row.names(S3), S3)
            
          } else if (method=='deseq2_cov'){
            S3 <-read.table(paste0('data/',simu,'/','all_deseq2_cov_table.txt'), head=T,  sep='\t', fill=T)
            S3 <- data.frame(genes = row.names(S3), S3)
            
          } else if (method=='edger'){
            S3 <-read.table(paste0('data/',simu,'/','all_edger_result_table.txt'), head=T,  sep='\t', fill=T)
            S3 <- data.frame(genes = row.names(S3), S3)
            
          } else if (method=='edger_cov'){
            S3 <-read.table(paste0('data/',simu,'/','all_edger_cov_table.txt'), head=T,  sep='\t', fill=T)
            S3 <- data.frame(genes = row.names(S3), S3)
            
          } else if (method=='edger_detrate'){
            S3 <-read.table(paste0('data/',simu,'/','all_edger_detrate_table.txt'), head=T,  sep='\t', fill=T)
            S3 <- data.frame(genes = row.names(S3), S3)
            
          }else if (method=='edger_detrate_cov'){
            S3 <-read.table(paste0('data/',simu,'/','all_edger_detrate_cov_table.txt'), head=T,  sep='\t', fill=T)
            S3 <- data.frame(genes = row.names(S3), S3)
            
          }  else if (method=='limma_voom'){
            S3 <-read.table(paste0('data/',simu,'/','voom_result_table.txt'), head=T,  sep='\t', fill=T)
            S3 <- data.frame(genes = row.names(S3), S3)
            
          } else if (method=='limma_voom_cov'){
            S3 <-read.table(paste0('data/',simu,'/','voom_cov_table.txt'), head=T,  sep='\t', fill=T)
            S3 <- data.frame(genes = row.names(S3), S3)
            
          } else if (method=='limma_trend'){
            S3 <-read.table(paste0('data/',simu,'/','limma_trend_result.txt'), head=T,  sep='\t', fill=T)
            S3 <- data.frame(genes = row.names(S3), S3)
            
          } else if (method=='limma_trend_cov'){
            S3 <-read.table(paste0('data/',simu,'/','limma_trend_cov_result.txt'), head=T,  sep='\t', fill=T)
            S3 <- data.frame(genes = row.names(S3), S3)
            
          } else if (method=='limma_trend_false'){
            S3 <-read.table(paste0('data/',simu,'/','limma_trend_result_false.txt'), head=T,  sep='\t', fill=T)
            S3 <- data.frame(genes = row.names(S3), S3)
            
          } else if (method=='limma_trend_false_cov'){
            S3 <-read.table(paste0('data/',simu,'/','limma_trend_cov_false.txt'), head=T,  sep='\t', fill=T)
            S3 <- data.frame(genes = row.names(S3), S3)
            
          } else if (method=='mnn_limma_trend'){
            S3 <-read.table(paste0('data/',simu,'/','limma_trend_MNN_result.txt'), head=T,  sep='\t', fill=T)
            S3 <- data.frame(genes = row.names(S3), S3)
            
          } else if (method=='mnn_limma_trend_false'){
            S3 <-read.table(paste0('data/',simu,'/','limma_trend_MNN_result_false.txt'), head=T,  sep='\t', fill=T)
            S3 <- data.frame(genes = row.names(S3), S3)
            
          } else if (method=='scmerge_limma_trend'){
            S3 <-read.table(paste0('data/',simu,'/','limma_trend_scMerge_result.txt'), head=T,  sep='\t', fill=T)
            S3 <- data.frame(genes = row.names(S3), S3)
            
          } else if (method=='scmerge_limma_trend_false'){
            S3 <-read.table(paste0('data/',simu,'/','limma_trend_scMerge_result_false.txt'), head=T,  sep='\t', fill=T)
            S3 <- data.frame(genes = row.names(S3), S3)
            
          } else if (method=='combat_limma_trend'){
            S3 <-read.table(paste0('data/',simu,'/','limma_trend_combat_result.txt'), head=T,  sep='\t', fill=T)
            S3 <- data.frame(genes = row.names(S3), S3)
            
          } else if (method=='combat_limma_trend_false'){
            S3 <-read.table(paste0('data/',simu,'/','limma_trend_combat_result_false.txt'), head=T,  sep='\t', fill=T)
            S3 <- data.frame(genes = row.names(S3), S3)
            
          } else if (method=='zinbwave_deseq2'){
            S3 <-read.table(paste0('data/',simu,'/','all_zinbwave_deseq2_table.txt'), head=T,  sep='\t', fill=T)
            S3 <- data.frame(genes = row.names(S3), S3)
            
          } else if (method=='zinbwave_deseq2_cov'){
            S3 <-read.table(paste0('data/',simu,'/','all_zinbwave_deseq2_table_cov.txt'), head=T,  sep='\t', fill=T)
            S3 <- data.frame(genes = row.names(S3), S3)
            
          } else if (method=='zinbwave_edger'){
            S3 <-read.table(paste0('data/',simu,'/','all_zinbwave_edger_table.txt'), head=T,  sep='\t', fill=T)
            S3 <- data.frame(genes = row.names(S3), S3)
            
          } else if (method=='zinbwave_edger_cov'){
            S3 <-read.table(paste0('data/',simu,'/','all_zinbwave_edger_table_cov.txt'), head=T,  sep='\t', fill=T)
            S3 <- data.frame(genes = row.names(S3), S3)
            
          } else if (method=='mast'){
            S3 <-read.table(paste0('data/',simu,'/','all_MAST_result_table.txt'), head=T, fill=T)
            S3 <- data.frame(X = row.names(S3), S3)
            
          } else if (method=='mast_cov'){
            S3 <-read.table(paste0('data/',simu,'/','all_MAST_Cov_table.txt'), head=T, fill=T)
            S3 <- data.frame(X = row.names(S3), S3)
            
          } else if (method %in% meta_method){
            S3 <- data.frame()
          } else if (method=='limma_bec'){
            S3 <- read.table(paste0('Seurat_DEGs_auc/',simu,'/','S3_batch12/','after_', 'limma','_',HVG,'/degs_batch12_seurat_auc_DEG.txt'), head=T, sep='\t', fill=T)
          }
            else {
            S3 <- read.table(paste0('Seurat_DEGs_auc/',simu,'/','S3_batch12/','after_',method,'_',HVG,'/degs_batch12_seurat_auc_DEG.txt'), head=T, sep='\t', fill=T)
          }
        }
        geneinfo <- read.table(paste0('data/',simu,'/geneinfo.txt'), head=T, fill=T)
        colnames(S3)[colnames(S3) == 'avg_log2FC'] <- 'avg_logFC'
        
        cutoff = 0.05
        if(select=='UP'){
          if (method %in% c('deseq2', 'edger', 'limma_voom', 'deseq2_cov', 'edger_cov',
                            'zinbwave_deseq2', 'zinbwave_edger', 'zinbwave_deseq2_cov', 'zinbwave_edger_cov', 
                            'edger_detrate', 'edger_detrate_cov',
                            "limma_voom_cov", "limma_trend_cov", "limma_trend_false_cov",
                            'limma_trend', 'scmerge_limma_trend', 'mnn_limma_trend',
                            "combat_limma_trend", "combat_limma_trend_false",
                            'limma_trend_false', 'scmerge_limma_trend_false', 'mnn_limma_trend_false')){
            S3 <- S3[S3$logFC<0,]
          }  else if (method %in% c(df_med, newmeta_med)){
            
            load(paste0('data/',simu,'/','2batch_righttailed_Meta_Res.RData'))
            tmp <- MetaDE.Res[[method]]
            
            if (is.null(rownames(tmp$FDR))){
              S3 <- data.frame(Gene=names(tmp$FDR), pval=tmp$FDR)
            } else {
              S3 <- data.frame(Gene=rownames(tmp$FDR), pval=tmp$FDR)
            }
            colnames(S3) <- c('Gene', 'pval') 
            
          }  else if (method %in% c(meta_med) ){
            load(paste0('data/',simu,'/','2batch_righttailed_Meta_Res.RData'))
            tmp <- MetaDE.Res[[method]]
            
            if (is.null(rownames(tmp$meta.analysis$FDR))){
              S3 <- data.frame(Gene=names(tmp$meta.analysis$FDR), pval=tmp$meta.analysis$FDR)
            } else {
              S3 <- data.frame(Gene=rownames(tmp$meta.analysis$FDR), pval=tmp$meta.analysis$FDR)
            }
            colnames(S3) <- c('Gene', 'pval')
            
          }
          else {
            S3 <- S3[S3$avg_logFC>0,]
          }
          
          if (method %in% df_med){
            up_genes <- read.table(paste0('data/',simu,'/true_down_genes.txt'), head=T, fill = T)
            S7 <- as.data.frame(geneinfo[which(geneinfo$x %in% up_genes$x),])
            # S7 <- as.data.frame(geneinfo[which(geneinfo$Gene %in% up_genes$x),])
            # rownames(S7) <- S7[[1]]
            rownames(S7) <- S7[[1]]
            colnames(S7) <- c('Gene')
          } else {
            up_genes <- read.table(paste0('data/',simu,'/true_up_genes.txt'), head=T, fill = T)
            S7 <- as.data.frame(geneinfo[which(geneinfo$x %in% up_genes$x),])
            # S7 <- as.data.frame(geneinfo[which(geneinfo$Gene %in% up_genes$x),])
            # rownames(S7) <- S7[[1]]
            rownames(S7) <- S7[[1]]
            colnames(S7) <- c('Gene')
          }
          
          
        } else if(select=='DOWN'){
          if (method %in% c('deseq2', 'edger', 'limma_voom', 'deseq2_cov', 'edger_cov',
                            'zinbwave_deseq2', 'zinbwave_edger',  'zinbwave_deseq2_cov', 'zinbwave_edger_cov', 
                            'edger_detrate', 'edger_detrate_cov',
                            "limma_voom_cov", "limma_trend_cov", "limma_trend_false_cov",
                            'limma_trend', 'scmerge_limma_trend', 'mnn_limma_trend',
                            "combat_limma_trend", "combat_limma_trend_false",
                            'limma_trend_false', 'scmerge_limma_trend_false', 'mnn_limma_trend_false')){
            S3 <- S3[S3$logFC>0,]
            
          } else if (method %in% c(df_med, newmeta_med)){
            
            load(paste0('data/',simu,'/','2batch_lefttailed_Meta_Res.RData'))  
            tmp <- MetaDE.Res[[method]]
            
            if (is.null(rownames(tmp$FDR))){
              S3 <- data.frame(Gene=names(tmp$FDR), pval=tmp$FDR)
            } else {
              S3 <- data.frame(Gene=rownames(tmp$FDR), pval=tmp$FDR)
            }
            colnames(S3) <- c('Gene', 'pval')
            
          } else if (method %in% meta_med){
            
            load(paste0('data/',simu,'/','2batch_lefttailed_Meta_Res.RData'))  
            tmp <- MetaDE.Res[[method]]
            
            if (is.null(rownames(tmp$meta.analysis$FDR))){
              S3 <- data.frame(Gene=names(tmp$meta.analysis$FDR), pval=tmp$meta.analysis$FDR)
            } else {
              S3 <- data.frame(Gene=rownames(tmp$meta.analysis$FDR), pval=tmp$meta.analysis$FDR)
            }
            colnames(S3) <- c('Gene', 'pval')
            
          } 
          else {
            S3 <- S3[S3$avg_logFC<0,]
          }
          
          if (method %in% df_med) {
            down_genes <- read.table(paste0('data/',simu,'/true_up_genes.txt'), head=T, fill = T)
            S7 <- as.data.frame(geneinfo[which(geneinfo$x %in% down_genes$x),])
            # S7 <- as.data.frame(geneinfo[which(geneinfo$Gene %in% down_genes$x),])
            rownames(S7) <- S7[[1]]
            colnames(S7) <- c('Gene')
          } else{
            down_genes <- read.table(paste0('data/',simu,'/true_down_genes.txt'), head=T, fill = T)
            S7 <- as.data.frame(geneinfo[which(geneinfo$x %in% down_genes$x),])
            # S7 <- as.data.frame(geneinfo[which(geneinfo$Gene %in% down_genes$x),])
            rownames(S7) <- S7[[1]]
            colnames(S7) <- c('Gene')
          }
          
        } else {
          stop('select UP or DOWN DEGs')
        }
        
        N <- geneinfo$x
        # N <- geneinfo$Gene
        
        if (method %in% c('deseq2', 'edger', 'limma_voom', 'deseq2_cov', 'edger_cov',
                          'zinbwave_deseq2', 'zinbwave_edger',  'zinbwave_deseq2_cov', 'zinbwave_edger_cov', 
                          'edger_detrate', 'edger_detrate_cov',
                          "limma_voom_cov", "limma_trend_cov", "limma_trend_false_cov",
                          'limma_trend', 'scmerge_limma_trend', 'mnn_limma_trend',
                          "combat_limma_trend", "combat_limma_trend_false",
                          'limma_trend_false', 'scmerge_limma_trend_false', 'mnn_limma_trend_false')){
          S3 <- S3[S3$adjpvalue<=cutoff,]
          norm <- S3$genes
          
        } else if (method %in% meta_method){
          S3 <- S3[S3$pval <= cutoff,]
          norm <- S3$Gene
        }
        else {
          S3 <- S3[S3$p_val_adj<=cutoff,]
          norm <- S3$X
        }
        
        GT <- S7$Gene
        # table
        TP <- sum(GT %in% norm)
        FP <- length(norm)-TP
        FN <- length(GT) - TP
        # TN <- length(N) - length(norm) - FN
        TN <- length(N) - TP - FP - FN
        
        # recall (sensitivity)
        TPR <- round(TP/(TP + FN),3)
        
        # precision (positive predictive value)
        PPV <- round(TP/(TP + FP),3)
        
        # F-beta 
        beta = 1/2
        beta_2= beta^2
        
        Fscore <- (1 + beta_2)*(PPV*TPR)/(beta_2*PPV + TPR)
        Fscore <- round(Fscore,3)
        
        # table matrix confusion 
        data <- data.frame(TP=TP,FN=FN,FP=FP,recall=TPR,precision=PPV,Fscore=Fscore, row.names = method)
        # data <- data.frame(TP=TP,FN=FN,FP=FP,recall=TPR,precision=PPV,Fscore=FP, row.names = method)
        data[is.na(data)]=0
        return(data)
      })
      
      mef <- t(df)
      rownames(mef) <-  vect_method 
      
      mef <- rbind(rep('',dim(mef)[2]),mef)
      rownames(mef)[rownames(mef)==""] <- simu
      # print (mef)
      return(list(mef,df['Fscore',]))
      
    })
    
    df_all_all <- do.call(rbind,lapply(df_all,function(l){l[[1]]}))
    
    # add Average block
    df_fscore <- do.call(rbind,lapply(df_all,function(l){unlist(l[[2]])}))
    
    write.table(df_all_all,paste0(base_name,'summary_confusion_matrix_',HVG,'_DEGs_',select,'.txt'),sep='\t', quote=F)
    return(df_fscore)
    
  })
  names(Fscore_list) <- vect_HVG
  return(Fscore_list)
}

######################### BOXPLOT OF FSCORE
library(ggplot2)
library(cowplot)

# Prepare the data
Fscore_up <- main_Fscore(select='UP') # comment write.table(df_all_all)
Fscore_down <- main_Fscore(select='DOWN') # comment write.table(df_all_all)

Fall <- rbind(Fscore_up$all, Fscore_down$all)
# F_all <- reshape2::melt(list(Fscore_up$all,Fscore_down$all), value.name = "F.score")
F_all <- reshape2::melt(list(Fall), value.name = "F.score")
F_all <- F_all[, -1]

F_all$Var2 <- factor(F_all$Var2, levels = unique(F_all$Var2))
F_all$Var2 <- relevel(F_all$Var2, ref = "Combat")
F_all$Var2 <- relevel(F_all$Var2, ref = "limma_bec")
F_all$Var2 <- relevel(F_all$Var2, ref = "MNN")
F_all$Var2 <- relevel(F_all$Var2, ref = "scmerge")
F_all$Var2 <- relevel(F_all$Var2, ref = "seurat3")
F_all$Var2 <- relevel(F_all$Var2, ref = "zinbwave")
F_all$Var2 <- relevel(F_all$Var2, ref = "raw_data")

F_all$Var2 <- relevel(F_all$Var2, ref = "mast")
F_all$Var2 <- relevel(F_all$Var2, ref = "mast_cov")

F_all$Var2 <- relevel(F_all$Var2, ref = "deseq2")
F_all$Var2 <- relevel(F_all$Var2, ref = "deseq2_cov")
F_all$Var2 <- relevel(F_all$Var2, ref = "zinbwave_deseq2")
F_all$Var2 <- relevel(F_all$Var2, ref = "zinbwave_deseq2_cov")
F_all$Var2 <- relevel(F_all$Var2, ref = "deseq2_Fisher")
F_all$Var2 <- relevel(F_all$Var2, ref = "deseq2_wFisher_sample+gene")

F_all$Var2 <- relevel(F_all$Var2, ref = "edger_detrate")
F_all$Var2 <- relevel(F_all$Var2, ref = "edger_detrate_cov")
F_all$Var2 <- relevel(F_all$Var2, ref = "edger")
F_all$Var2 <- relevel(F_all$Var2, ref = "edger_cov")
F_all$Var2 <- relevel(F_all$Var2, ref = "zinbwave_edger")
F_all$Var2 <- relevel(F_all$Var2, ref = "zinbwave_edger_cov")
F_all$Var2 <- relevel(F_all$Var2, ref = "edgeR_Fisher")
F_all$Var2 <- relevel(F_all$Var2, ref = "edgeR_wFisher_sample+gene")


F_all$Var2 <- relevel(F_all$Var2, ref = "limma_voom")
F_all$Var2 <- relevel(F_all$Var2, ref = "limma_voom_cov")
F_all$Var2 <- relevel(F_all$Var2, ref = "limma_trend_false")
F_all$Var2 <- relevel(F_all$Var2, ref = "limma_trend_false_cov")
F_all$Var2 <- relevel(F_all$Var2, ref = "combat_limma_trend_false")
F_all$Var2 <- relevel(F_all$Var2, ref = "mnn_limma_trend_false")
F_all$Var2 <- relevel(F_all$Var2, ref = "scmerge_limma_trend_false")
F_all$Var2 <- relevel(F_all$Var2, ref = "LogNormalize+limma.trend_Fisher")
F_all$Var2 <- relevel(F_all$Var2, ref = "LogNormalize+limma.trend_wFisher_sample+gene")


F_all$Var2 <- relevel(F_all$Var2, ref = "deseq2+FEM")
F_all$Var2 <- relevel(F_all$Var2, ref = "voom+FEM")
F_all$Var2 <- relevel(F_all$Var2, ref = "LogNormalize+FEM")

F_all$Var2 <- relevel(F_all$Var2, ref = "deseq2+REM")
F_all$Var2 <- relevel(F_all$Var2, ref = "voom+REM")
F_all$Var2 <- relevel(F_all$Var2, ref = "LogNormalize+REM")

F_all$Var2 <- relevel(F_all$Var2, ref = "voom+modt_Fisher")
F_all$Var2 <- relevel(F_all$Var2, ref = "voom+modt_wFisher_sample+gene")

t = length(unique(F_all$Var2))
meanplot = aggregate(F_all[2],list(rep(1:(nrow(F_all[2])%/%12+1),each=12,len=nrow(F_all))),median)[-1]

rawmedian = rep(c(meanplot$F.score[13], meanplot$F.score[13 + t],
                  meanplot$F.score[13 + 2*t], meanplot$F.score[13 + 3*t]),
                each=nrow(F_all[2])%/%4)

F_all$RawMedian = as.numeric(rawmedian)

library(plyr)
df <- ddply(F_all, .(L1), summarise, median=median(RawMedian, na.rm = TRUE))

RAW <- c('raw_data')
BEC <- c('Combat', 'limma_bec',
         'MNN', 'mnn_limma_trend_false', 
         'scmerge', 'scmerge_limma_trend_false',
         "combat_limma_trend_false",
         'seurat3', 'zinbwave', 'zinbwave_deseq2', 'zinbwave_edger',
         'zinbwave_deseq2_cov', 'zinbwave_edger_cov'
         )

COV <- c('mast', 'mast_cov', 
         'deseq2', 'deseq2_cov', 
         'edger', 'edger_cov', 'edger_detrate', 'edger_detrate_cov',
         "limma_voom", "limma_voom_cov",
         "limma_trend_false", "limma_trend_false_cov"
         )

META<- c( 
  "voom+REM"                                          ,"voom+FEM"                                         ,
  "LogNormalize+REM"                                  ,"LogNormalize+FEM"                                  ,
  "deseq2+REM"                                        ,"deseq2+FEM"                ,
  "voom+modt_Fisher"                                  ,"voom+modt_wFisher_sample+gene"                     ,
  "edgeR_Fisher"                                      ,"edgeR_wFisher_sample+gene"                        ,
  "deseq2_Fisher"                                     ,"deseq2_wFisher_sample+gene"                        ,
  "LogNormalize+limma.trend_Fisher"                   ,"LogNormalize+limma.trend_wFisher_sample+gene"      
)

mybreaks <- c("Combat", "limma_bec",  "MNN",  "scmerge",  "seurat3",  "zinbwave",  "raw_data",
              "mast",  "mast_cov",  "deseq2",  "deseq2_cov",  "zinbwave_deseq2",  "zinbwave_deseq2_cov",
              "deseq2_Fisher",  "deseq2_wFisher_sample+gene",
              "edger_detrate",  "edger_detrate_cov",  "edger",  "edger_cov",  "zinbwave_edger",  "zinbwave_edger_cov",
              "edgeR_Fisher",  "edgeR_wFisher_sample+gene",
              "limma_voom",  "limma_voom_cov",  "limma_trend_false",  "limma_trend_false_cov",
              "combat_limma_trend_false",  "mnn_limma_trend_false",  "scmerge_limma_trend_false",
              "LogNormalize+limma.trend_Fisher",  "LogNormalize+limma.trend_wFisher_sample+gene",
  
              "deseq2+FEM",  "voom+FEM",  "LogNormalize+FEM",  "deseq2+REM",  "voom+REM",  "LogNormalize+REM",
              "voom+modt_Fisher",  "voom+modt_wFisher_sample+gene")

a <- ifelse(mybreaks %in% BEC, 'red4', 
            ifelse(mybreaks %in% COV, 'blue3',
                   ifelse(mybreaks %in% META, 'darkgreen', 
                          ifelse(mybreaks %in% RAW, 'black', 'white'))))

p4 <- ggplot(F_all, aes(x=Var2, y=F.score, color=Var2)) + 
  geom_boxplot(outlier.size = 1) + coord_flip() +
  geom_hline(data=df, aes(yintercept=median),linetype="dashed", color='black') + 
  # labs(y = 'F-score') + theme(axis.text.y = element_text(size=14, color = a))
  labs(y = 'F-score (beta=0.5)') + theme(axis.text.y = element_text(size=14, color = a))
  # labs(y = 'False Positive Count') + theme(axis.text.y = element_text(size=14, color = a))

p4 <- p4 + theme(axis.line = element_line(color = "black", size = 0.5, linetype = "solid"),
                 axis.title.x = element_text(size=16),
                 axis.title.y = element_blank(),
                 axis.text.x = element_text(size=14, colour = 'black', hjust = 1),
                 panel.grid.major = element_line(colour = "grey86"),
                 panel.grid.major.x = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 panel.grid.minor.y = element_blank(),
                 panel.background = element_blank(),
                 panel.border = element_blank(),
                 panel.spacing.x = unit(0.5, "lines"),
                 panel.spacing.y = unit(1, "lines"),
                 strip.text = element_text(size=17, color="black"),
                 strip.background.x = element_rect(fill="#CDE8DF"),
                 # strip.background.x = element_blank(),
                 # strip.background.y = element_blank(),
                 legend.position="none") 


p4 <- p4 + scale_x_discrete(breaks=mybreaks,

                            labels=c("Combat_Wilcox", "limma_BEC_Wilcox","MNNCorrect_Wilcox",
                                     "scMerge_Wilcox", "Seurat_Wilcox", "ZINB-WaVE_Wilcox",
                                     "Raw_Wilcox", 'MAST', 'MAST_Cov',
                                     'DESeq2','DESeq2_Cov', "ZINB-WaVE_DESeq2", "ZINB-WaVE_DESeq2_Cov",
                                     "DESeq2_Fisher", "DESeq2_wFisher",
                                     
                                     'edgeR_DetRate', 'edgeR_DetRate_Cov', 
                                     'edgeR', 'edgeR_Cov', "ZINB-WaVE_edgeR",  "ZINB-WaVE_edgeR_Cov",
                                     "edgeR_Fisher", "edgeR_wFisher",
                                     
                                     "limma", 'limma_Cov', 'limmatrend','limmatrend_Cov',
                                     'Combat_limmatrend', 'MNNCorrect_limmatrend', 'scMerge_limmatrend',
                                     "LogNorm+limmatrend_Fisher", "LogNorm+limmatrend_wFisher", 
                                     
                                     "DESeq2+FEM", "LogNorm+FEM", "voom+FEM",
                                     "DESeq2+REM", "LogNorm+REM", "voom+REM", 
                                     
                                     "voom+modt_Fisher", "voom+modt_wFisher" ))

p5 <- p4 + geom_point() #geom_jitter(shape=16, position=position_jitter(0.2),size=1)
     

p5
# ggsave(filename=paste0('Venn_diagram/Fscore_boxplot_flip2.png'), plot=p4, width = 18, height = 16)
tiff(filename = "Venn_diagram/Fscore_allgenes.tiff", units="px", width=900, height=800)
dev.off()