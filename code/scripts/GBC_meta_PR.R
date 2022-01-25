########## summary confusion matrix

rm(list=ls())

# load libraries
library(gridExtra)
library(openxlsx)
dir.create('Venn_diagram_PR_allgenes')

p_threshold <- c(0, 10^-4, 0.001, 0.005, 0.01, 0.05, 0.0500001, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.1)

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


# vect_method <- c('Combat', 'limma',
#                  'MNN', 'mnn_limma_trend', 'mnn_limma_trend_false',
#                  'scmerge', 'scmerge_limma_trend', 'scmerge_limma_trend_false',
#                  'seurat3', 'zinbwave', 'zinbwave_deseq2', 'zinbwave_edger', 
#                  'zinbwave_deseq2_cov', 'zinbwave_edger_cov',    
#                  'raw_data', 'mast', 'mast_cov', 
#                  'deseq2', 'deseq2_cov', 
#                  'edger', 'edger_cov', 'edger_detrate', 'edger_detrate_cov',
#                  "limma_voom", "limma_voom_cov", "limma_trend", "limma_trend_cov", 
#                  "limma_trend_false", "limma_trend_false_cov",
#                  "combat_limma_trend_false",
#                  
#                  "voom+REM"                                          ,"voom+FEM"                                         ,
#                  "LogNormalize+REM"                                  ,"LogNormalize+FEM"                                  ,
#                  "deseq2+REM"                                        ,"deseq2+FEM"                ,
#                  "voom+modt_Fisher"                                  ,"voom+modt_wFisher_sample+gene"                     ,
#                  "edgeR_Fisher"                                      ,"edgeR_wFisher_sample+gene"                        ,
#                  "deseq2_Fisher"                                     ,"deseq2_wFisher_sample+gene"                        ,
#                  "LogNormalize+limma.trend_Fisher"                   ,"LogNormalize+limma.trend_wFisher_sample+gene"      
# )
vect_method <- c("Combat", "limma_bec",  "MNN",  "scmerge",  "seurat3",  "zinbwave",  "raw_data",
                 "mast",  "mast_cov",  "deseq2",  "deseq2_cov",  "zinbwave_deseq2",  "zinbwave_deseq2_cov",
                 "deseq2_Fisher",  "deseq2_wFisher_sample+gene",
                 "edger_detrate",  "edger_detrate_cov",  "edger",  "edger_cov",  "zinbwave_edger",  "zinbwave_edger_cov",
                 "edgeR_Fisher",  "edgeR_wFisher_sample+gene",
                 "limma_voom",  "limma_voom_cov",  "limma_trend_false",  "limma_trend_false_cov",
                 "combat_limma_trend_false",  "mnn_limma_trend_false",  "scmerge_limma_trend_false",
                 "LogNormalize+limma.trend_Fisher",  "LogNormalize+limma.trend_wFisher_sample+gene",
                 
                 "deseq2+FEM",   "LogNormalize+FEM", "voom+FEM",
                 "deseq2+REM",   "LogNormalize+REM", "voom+REM",
                 "voom+modt_Fisher",  "voom+modt_wFisher_sample+gene")
vect_method <- rev(vect_method)

main_Fscore <- function(select){
  # METHODS
  
  # HVG (all/as Seurat)
  vect_HVG <- c('all')
  
  # SIMULATIONS
  vect_simu <- gsub('data/','',list.dirs('data', recursive=FALSE))
  
  Fscore_list <- lapply(vect_HVG,function(HVG){
    base_name <- paste0('Venn_diagram_PR_allgenes/')
    
    df_all <- lapply(p_threshold, function(cutoff){
      
      df_simu<- lapply(vect_simu, function(simu){
        
        df <- sapply(vect_method, function(method){
          
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
              S3 <- read.table(paste0('Seurat_DEGs_auc/',simu,'/','S3_batch12/','after_','limma','_',HVG,'/degs_batch12_seurat_auc_DEG.txt'), head=T, sep='\t', fill=T)
            } else {
              S3 <- read.table(paste0('Seurat_DEGs_auc/',simu,'/','S3_batch12/','after_',method,'_',HVG,'/degs_batch12_seurat_auc_DEG.txt'), head=T, sep='\t', fill=T)
            }
          }
          geneinfo <- read.table(paste0('data/',simu,'/geneinfo.txt'), head=T, fill=T)
          colnames(S3)[colnames(S3) == 'avg_log2FC'] <- 'avg_logFC'
          
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
                S3 <- data.frame(Gene=names(tmp$FDR), pval=tmp$pval, qval=tmp$FDR)
              } else {
                S3 <- data.frame(Gene=rownames(tmp$FDR),pval=tmp$pval, qval=tmp$FDR)
              }
              colnames(S3) <- c('Gene', 'pval', 'qval')
              
            }  else if (method %in% c(meta_med) ){
              load(paste0('data/',simu,'/','2batch_righttailed_Meta_Res.RData'))
              tmp <- MetaDE.Res[[method]]
              
              if (is.null(rownames(tmp$meta.analysis$FDR))){
                S3 <- data.frame(Gene=names(tmp$meta.analysis$FDR), pval=tmp$meta.analysis$pval, qval=tmp$meta.analysis$FDR)
              } else {
                S3 <- data.frame(Gene=rownames(tmp$meta.analysis$FDR), pval=tmp$meta.analysis$pval, qval=tmp$meta.analysis$FDR)
              }
              colnames(S3) <- c('Gene', 'pval', 'qval')
              
            }
            else {
              S3 <- S3[S3$avg_logFC>0,]
            }
            
            if (method %in% df_med){
            # if (method %in% c( #"voom+REM", "voom+FEM",
            #                     # "voom+modt_Fisher", "edgeR_Fisher",
            #                    "LogNormalize+REM", "LogNormalize+FEM",
            #                    "deseq2+REM", "deseq2+FEM")){
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
                S3 <- data.frame(Gene=names(tmp$FDR), pval=tmp$pval, qval=tmp$FDR)
              } else {
                S3 <- data.frame(Gene=rownames(tmp$FDR),pval=tmp$pval, qval=tmp$FDR)
              }
              colnames(S3) <- c('Gene', 'pval', 'qval')
              
            } else if (method %in% meta_med){
              
              load(paste0('data/',simu,'/','2batch_lefttailed_Meta_Res.RData'))  
              tmp <- MetaDE.Res[[method]]
              
              if (is.null(rownames(tmp$meta.analysis$FDR))){
                S3 <- data.frame(Gene=names(tmp$meta.analysis$FDR), pval=tmp$meta.analysis$pval, qval=tmp$meta.analysis$FDR)
              } else {
                S3 <- data.frame(Gene=rownames(tmp$meta.analysis$FDR), pval=tmp$meta.analysis$pval, qval=tmp$meta.analysis$FDR)
              }
              colnames(S3) <- c('Gene', 'pval', 'qval')
              
            } 
            else {
              S3 <- S3[S3$avg_logFC<0,]
            }
            
            if (method %in% df_med) {
            # if (method %in% c( #"voom+REM", "voom+FEM",
            #   "LogNormalize+REM", "LogNormalize+FEM",
            #   "deseq2+REM", "deseq2+FEM")){
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
         
          if (cutoff == 0.05){
            if (method %in% c('deseq2', 'edger', 'limma_voom', 'deseq2_cov', 'edger_cov',
                              'zinbwave_deseq2', 'zinbwave_edger',  'zinbwave_deseq2_cov', 'zinbwave_edger_cov',
                              'edger_detrate', 'edger_detrate_cov',
                              "limma_voom_cov", "limma_trend_cov", "limma_trend_false_cov",
                              'limma_trend', 'scmerge_limma_trend', 'mnn_limma_trend',
                              "combat_limma_trend", "combat_limma_trend_false",
                              'limma_trend_false', 'scmerge_limma_trend_false', 'mnn_limma_trend_false')){
              currGenes <- S3$genes
              missingGenes <- N[!(N %in% currGenes)]
              ml <- length(missingGenes)
              mp <- seq(1, 1, len=ml)

              new_df <- data.frame(genes=currGenes, p_val=S3$adjpvalue)
              compliment_df <- data.frame(genes=missingGenes, p_val=mp)

              S3 <- rbind(new_df, compliment_df)
              colnames(S3) <- c('genes', 'p_val')

            } else if (method %in% meta_method){
              currGenes <- S3$Gene
              missingGenes <- N[!(N %in% currGenes)]
              ml <- length(missingGenes)
              mp <- seq(1, 1, len=ml)

              new_df <- data.frame(genes=currGenes, p_val=S3$qval)
              compliment_df <- data.frame(genes=missingGenes, p_val=mp)

              S3 <- rbind(new_df, compliment_df)
              colnames(S3) <- c('genes', 'p_val')

            } else {
              currGenes <- S3$X
              missingGenes <- N[!(N %in% currGenes)]
              ml <- length(missingGenes)
              mp <- seq(1, 1, len=ml)

              new_df <- data.frame(genes=currGenes, p_val=S3$p_val_adj)
              compliment_df <- data.frame(genes=missingGenes, p_val=mp)

              S3 <- rbind(new_df, compliment_df)
              colnames(S3) <- c('genes', 'p_val')
            }
          } else {
            if (method %in% c('deseq2', 'edger', 'limma_voom', 'deseq2_cov', 'edger_cov',
                              'zinbwave_deseq2', 'zinbwave_edger',  'zinbwave_deseq2_cov', 'zinbwave_edger_cov', 
                              'edger_detrate', 'edger_detrate_cov',
                              "limma_voom_cov", "limma_trend_cov", "limma_trend_false_cov",
                              'limma_trend', 'scmerge_limma_trend', 'mnn_limma_trend',
                              "combat_limma_trend", "combat_limma_trend_false",
                              'limma_trend_false', 'scmerge_limma_trend_false', 'mnn_limma_trend_false')){
              currGenes <- S3$genes
              missingGenes <- N[!(N %in% currGenes)]
              ml <- length(missingGenes)
              mp <- seq(1, 1, len=ml)
              
              new_df <- data.frame(genes=currGenes, p_val=S3$pvalue)
              compliment_df <- data.frame(genes=missingGenes, p_val=mp)
              
              S3 <- rbind(new_df, compliment_df)
              colnames(S3) <- c('genes', 'p_val')
              
            } else if (method %in% meta_method){
              currGenes <- S3$Gene
              missingGenes <- N[!(N %in% currGenes)]
              ml <- length(missingGenes)
              mp <- seq(1, 1, len=ml)
              
              new_df <- data.frame(genes=currGenes, p_val=S3$pval)
              compliment_df <- data.frame(genes=missingGenes, p_val=mp)
              
              S3 <- rbind(new_df, compliment_df)
              colnames(S3) <- c('genes', 'p_val')
              
            } else {
              currGenes <- S3$X
              missingGenes <- N[!(N %in% currGenes)]
              ml <- length(missingGenes)
              mp <- seq(1, 1, len=ml)
              
              new_df <- data.frame(genes=currGenes, p_val=S3$p_val)
              compliment_df <- data.frame(genes=missingGenes, p_val=mp)
              
              S3 <- rbind(new_df, compliment_df)
              colnames(S3) <- c('genes', 'p_val')
            }
          }
          S3 <- S3[S3$p_val<cutoff,]
          
          GT <- S7$Gene
          norm <- S3$genes
          
         # table
          TP <- sum(GT %in% norm)
          FP <- length(norm)-TP
          FN <- length(GT) - TP
          TN <- length(N) - TP - FP - FN
          
          # sensitivity/recall/true positive rate
          TPR <- round(TP/(TP + FN), 3)
          # specificity/ true negative rate
          TNR <- round(TN/(TN + FP), 3)
          # precision (positive predictive value)
          PPV <- round(TP/(TP + FP), 3)
          
          # F-score
          Fscore <- 2*((PPV*TPR)/(PPV+TPR))
          Fscore <- round(Fscore, 3)
          
          data <- data.frame(TP=TP, FN=FN, FP=FP, TN=TN, 
                             PPV=PPV, TPR=TPR, TNR=TNR, Fscore=Fscore,
                             cutoff=cutoff, simu=simu)
          data[is.na(data)] <- 0.
          return(data)
        })
        return(t(df))
      })
      
      df_all_simu <- do.call(rbind,df_simu)
      return(df_all_simu)
    })
    
    df_all_all <- do.call(rbind,df_all)
    write.table(df_all_all,paste0(base_name,'summary_confusion_matrix_',HVG,'_DEGs_',select,'.txt'),sep='\t', quote=F)
    return(df_all_all)
  })
  names(Fscore_list) <- vect_HVG
  return(Fscore_list)
}


######################### BOXPLOT OF FSCORE
library(ggplot2)

Score_up <- main_Fscore(select='UP') # comment write.table(df_all_all)
Score_down <- main_Fscore(select='DOWN') 

legend_titles <- c("Combat_Wilcox", "limma_BEC_Wilcox","MNNCorrect_Wilcox",
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
                   
                   "DESeq2+FEM",  "LogNormalize+FEM", "voom+FEM",
                   "DESeq2+REM",  "LogNormalize+REM", "voom+REM",
                   "voom+modt_Fisher", "voom+modt_wFisher" )
legend_titles <- rev(legend_titles)

for (mtd in vect_method){
  rownames(Score_up$all)[rownames(Score_up$all) == mtd] <- legend_titles[vect_method==mtd]
  rownames(Score_down$all)[rownames(Score_down$all) == mtd] <- legend_titles[vect_method==mtd]
}

Score_up$all <- data.frame(method = rownames(Score_up$all), Score_up$all)
Score_down$all <- data.frame(method = rownames(Score_down$all), Score_down$all)

Fall <- rbind(Score_up$all, Score_down$all)

library(ggplot2)

#############################################33
simple_auc <- function(PPV, TPR){
  # inputs already sorted, best scores first 
  dPPV <- c(diff(PPV), 0)
  dTPR <- c(diff(TPR), 0)
  return(sum(dTPR * PPV) + sum(dTPR * dPPV)/2)
}

# ######################## all
up_all <- data.frame(matrix(ncol = 4, nrow = 0))
up_all_005 <- data.frame(matrix(ncol = 4, nrow = 0))
tmp_up_all <- data.frame(matrix(ncol = 4, nrow = 0))

for (mtd in legend_titles){
  # mtd <- "edgeR_wFisher"
  # mtd <- ""
  tmp <- Fall[Fall$method==mtd,]
  ppv <- tpr <- c()
  for (cutoff in p_threshold){
    # cutoff <- p_threshold[1]
    tppv <- mean(unlist(tmp$PPV[tmp$cutoff==cutoff]))
    ttpr <- mean(unlist(tmp$TPR[tmp$cutoff==cutoff]))
    tppv[tppv==0] = 1   
    
    ppv <- c(ppv, tppv)
    tpr <- c(tpr, ttpr)
  }
  tmp_up_all <- rbind(tmp_up_all, 
                      data.frame(method=rep(mtd, length(ppv)), precision=ppv, recall=tpr, cutoff=p_threshold ))
}

colnames(tmp_up_all) <- c("method", "Precision", "Recall", "cutoff")
tmp_up_all$Precision<-as.numeric(tmp_up_all$Precision)
tmp_up_all$Recall<-as.numeric(tmp_up_all$Recall)

# ################################################ 
for (mtd in legend_titles){
  # mtd = 'ZINB-WaVE_DESeq2'
  dtList <- tmp_up_all[tmp_up_all$method==mtd,]
  prec <- dtList$Precision
  reca <- dtList$Recall
  cuto <- dtList$cutoff

  up_all <- rbind(up_all, dtList)
}
colnames(up_all) <- c("method", "Precision", "Recall", "cutoff")

################################################ 
Recall_cutoff = 0.5  ## TPR column, cutoff is changed to 1

half_curve <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(half_curve) <- c("method", "Precision", "Recall", "cutoff")

for (mtd in legend_titles){
  dtList <- up_all[up_all$method==mtd,]
  upperBound <- dtList[dtList$Recall > Recall_cutoff,]
  lowerBound <- dtList[dtList$Recall <= Recall_cutoff,]
  
  if (length(upperBound$Recall) == 0){
    recall <- c(rev(lowerBound$Recall)[1:2], 1)
    precision <- c(rev(lowerBound$Precision)[1:2], min(dtList$Precision),)
    ispl <- spline(recall, precision, xout = c(0.5), method = 'natural')
    
  } else {
    recall <- c(rev(lowerBound$Recall)[1:2], upperBound$Recall[1:2])
    precision <- c(rev(lowerBound$Precision)[1:2], upperBound$Precision[1:2])
    ispl <- spline(recall, precision, xout = c(0.5), method = 'natural')
  }
  
  lowerBound <- rbind(lowerBound, data.frame(method=mtd, Precision=ispl$y, Recall=ispl$x, cutoff=1))
  
  half_curve <- rbind(half_curve, lowerBound)
}

pAUPR_coef = 0.5 #0.5 # 0.25 for .75, 0.5 for 0.5

taupr <- c()
for (mtd in legend_titles){
  # mtd = 'ZINB-WaVE_DESeq2'
  tmp <- half_curve[half_curve$method==mtd,]
  ppv <- tmp[tmp$cutoff != 0.05,]$Precision
  tpr <- tmp[tmp$cutoff != 0.05,]$Recall
  taupr <- c(taupr, round(simple_auc(PPV = ppv, TPR = tpr) / pAUPR_coef, 3))
}
aupr_up_all <- data.frame(score= taupr, method=unique(up_all$method))

up_all$Precision<-as.numeric(up_all$Precision)
up_all$Recall<-as.numeric(up_all$Recall)

up_all_005 <- up_all[up_all$cutoff==0.05, ]

y <- seq(.975, 0.25, len=length(legend_titles))
# y <- seq(.625, 0., len=length(legend_titles))
sy <- y
sy[order(aupr_up_all$score, decreasing = T)] <- y
aupr_up_all <- data.frame(y=sy, aupr_up_all)

###  highlight cutoff = 0.05 - F-score
up_all$method <- factor(up_all$method, levels=unique(up_all$method))


# 
#
graphString <- paste0("ggplot(up_all, aes(x=Recall, y=Precision, color=method)) + geom_point() + geom_line() +",
                      " scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1)) + labs(y='Precision', x = 'Recall') +",
                      'theme(legend.position="none") +',
                      'geom_text( hjust = "right", size=3.5,
                                                    aes(
                                                          x = 1.,
                                                          y = .975 + (y[1]-y[2]),
                                                          label = "pAUPR",
                                                          color = NULL
                                                        )
                                                    )'
)

graphString <- paste0(graphString, '+geom_vline(aes(xintercept=0.5),linetype="dashed", color="black")')

for (mtd in legend_titles){
  graphString <- paste0 (graphString, '+ geom_text( hjust = "right", size=3.5,
                                                    aes(
                                                         x = 1.,
                                                          y = aupr_up_all$y[aupr_up_all$method=="',mtd,'"],
                                                          label = paste0("',mtd,'", " — " , format(round(aupr_up_all$score[aupr_up_all$method=="',mtd,'"], 3), nsmall=3) ),
                                                          color = "',mtd,'"
                                                        )
                                                    )'
  )
}
graphString <- paste0(graphString, '+geom_point(data=up_all_005, aes(x=Recall, y=Precision,color=method), shape=1, size=4)')
eval(parse(text=graphString)) 
