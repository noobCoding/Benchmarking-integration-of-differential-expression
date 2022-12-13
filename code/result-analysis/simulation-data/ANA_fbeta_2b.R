########## summary confusion matrix

rm(list=ls())

# load libraries
library(gridExtra)
library(openxlsx)

vect_method <-c(
  "Combat_Wilcox",
  "limma_BEC_Wilcox",
  "MNNCorrect_Wilcox",
  "scMerge_Wilcox",
  "Seurat_Wilcox",
  "ZW_Wilcox",
  "scVI_Wilcox",
  "scGen_Wilcox",
  "Scanorama_Wilcox",
  "Raw_Wilcox",
  
  "RISC_Wilcox",
  "RISC_QP",
  
  'Pseudobulk_DESeq2',
  'Pseudobulk_edgeR',
  'Pseudobulk_limma',
  'Pseudobulk_limmatrend',
  
  "MAST",
  "MAST_Cov",
  
  "DESeq2",
  "DESeq2_Cov",
  "ZW_DESeq2",
  "ZW_DESeq2_Cov",
  
  "edgeR_DetRate",
  "edgeR_DetRate_Cov",
  "edgeR",
  "edgeR_Cov",
  "ZW_edgeR",
  "ZW_edgeR_Cov",
  
  "limma",
  "limma_Cov",
  "limmatrend",
  "limmatrend_Cov",
  "Combat_limmatrend",
  "MNNCorrect_limmatrend",
  "scMerge_limmatrend",
  'scVI_limmatrend',
  'scGen_limmatrend',
  "Scanorama_limmatrend",
  "RISC_limmatrend",
  
  "DESeq2_FEM",
  "LogNorm_FEM",
  "DESeq2_REM",
  "LogNorm_REM",
  
  "DESeq2_wFisher",
  "edgeR_wFisher",
  "LogNorm+limmatrend_wFisher"
)

main_Fscore <- function(select){
# SIMULATIONS
  vect_simu <- gsub('data/','',list.dirs('data', recursive=FALSE))
  
  Fscore_list <- lapply(c('all'),function(HVG){
    
    df_all <- lapply(vect_simu,function(simu){
      # simu = 'simul1_dropout_37_b1_100_b2_150'
      result.list<-readRDS(paste0("data/sp80.", simu, '_full.RData'))
      
      df <- sapply(vect_method,function(method){
          
        cutoff = 0.05
        S3 <- result.list[[method]]
        
        if(select=='UP'){
          if (method %in% c( "Pseudobulk_DESeq2", "Pseudobulk_edgeR", "Pseudobulk_limma",
                             "Pseudobulk_limmatrend", "DESeq2", "DESeq2_Cov",
                             "ZW_DESeq2", "Pseudo_ZW_DESeq2_Cov", "edgeR",
                             "edgeR_Cov", "edgeR_DetRate", "edgeR_DetRate_Cov", "ZW_edgeR", "ZW_edgeR_Cov",            
                              "limma", "limma_Cov", "limmatrend",        
                             "limmatrend_Cov", "Combat_limmatrend", "MNNCorrect_limmatrend", "scMerge_limmatrend",  
                              "RISC_limmatrend", "RISC_NB", "RISC_QP", "scVI_limmatrend", "scGen_limmatrend",        
                             "Scanorama_limmatrend"))
            {
              S3 <- S3[S3$log2FC<0,]
            } else if (method %in% c("DESeq2_FEM", "voom_FEM", "LogNorm_FEM", "DESeq2_REM", "voom_REM", "LogNorm_REM")){
              S3 <- S3[S3$log2FC<0,]
              
            } else if (method %in% c('voom+modt_Fisher', "voom+modt_wFisher")){
              S3 <- S3[S3$log2FC>0,]
              
            } else if (method %in% c("DESeq2_wFisher","edgeR_wFisher", "LogNorm+limmatrend_wFisher")){
              
              S3 <- S3[S3$log2FC<0,]
            } else if (method %in% c('Scanorama_Wilcox', 'scVI_Wilcox', 'scGen_Wilcox')){
              S3 <- S3[S3$log2FC>0,]
            }
          else {
              
              S3 <- S3[S3$log2FC<0,]
            }
            S7 <-result.list$GeneInfo$Up  
          }
        
        else if(select=='DOWN'){
          
          if (method %in% c( "Pseudobulk_DESeq2", "Pseudobulk_edgeR", "Pseudobulk_limma",
                             "Pseudobulk_limmatrend", "DESeq2", "DESeq2_Cov",
                             "ZW_DESeq2", "Pseudo_ZW_DESeq2_Cov", "edgeR",
                             "edgeR_Cov", "edgeR_DetRate", "edgeR_DetRate_Cov", "ZW_edgeR", "ZW_edgeR_Cov",            
                             "limma", "limma_Cov", "limmatrend",        
                             "limmatrend_Cov", "Combat_limmatrend", "MNNCorrect_limmatrend", "scMerge_limmatrend",  
                             "RISC_limmatrend", "RISC_NB", "RISC_QP", "scVI_limmatrend", "scGen_limmatrend",        
                             "Scanorama_limmatrend"))
          {
            S3 <- S3[S3$log2FC>0,]
          } else if (method %in% c("DESeq2_FEM", "voom_FEM", "LogNorm_FEM", "DESeq2_REM", "voom_REM", "LogNorm_REM")){
            S3 <- S3[S3$log2FC>0,]
            
          } else if (method %in% c('voom+modt_Fisher', "voom+modt_wFisher")){
            S3 <- S3[S3$log2FC<0,]
            
          } else if (method %in% c("DESeq2_wFisher","edgeR_wFisher", "LogNorm+limmatrend_wFisher")){
            
            S3 <- S3[S3$log2FC>0,]
          } else if (method %in% c('Scanorama_Wilcox', 'scVI_Wilcox', 'scGen_Wilcox')){
            S3 <- S3[S3$log2FC<0,]
          }
          else {
            
            S3 <- S3[S3$log2FC>0,]
          }
          S7 <-result.list$GeneInfo$Down  
        }
        
        else {
          stop('select UP or DOWN DEGs')
        }
        
        N <- result.list$GeneInfo$All
        S3 <- S3[S3$adj.pval <= cutoff,]
        
        norm <- rownames(S3)
        GT <- S7
    
        # table
        TP <- sum(GT %in% norm)
        FP <- length(norm)-TP
        FN <- length(GT) - TP
        TN <- length(N) - TP - FP - FN
        
        # recall (sensitivity)
        TPR <- round(TP/(TP + FN),3)
        # specificity/ true negative rate
        TNR <- round(TN/(TN + FP), 3)
        # precision (positive predictive value)
        PPV <- round(TP/(TP + FP),3)
        
        # F-beta 
        beta = 1/2
        beta_2= beta^2
        
        Fscore <- (1 + beta_2)*(PPV*TPR)/(beta_2*PPV + TPR)
        Fscore <- round(Fscore,3)
        
        # table matrix confusion 
        data <- data.frame(TP=TP,FN=FN,FP=FP,precision=PPV,
                           recall=TPR, specificity=TNR,
                           Fscore=Fscore, row.names = method)
        data[is.na(data)]=0
        return(data)
      })
      
      mef <- t(df)
      rownames(mef) <-  vect_method 
      
      mef <- rbind(rep('',dim(mef)[2]),mef)
      rownames(mef)[rownames(mef)==""] <- simu
      return(list(mef,df['Fscore',]))
      
    })
    
    df_all_all <- do.call(rbind,lapply(df_all,function(l){l[[1]]}))
    
    # add Average block
    df_fscore <- do.call(rbind,lapply(df_all,function(l){unlist(l[[2]])}))
    
    write.table(df_all_all,paste0('summary_confusion_matrix_',select,'.txt'),sep='\t', quote=F)
    return(df_fscore)
    
  })
  names(Fscore_list) <- c('all')
  return(Fscore_list)
}

######################### BOXPLOT OF FSCORE
library(ggplot2)
library(cowplot)

# Prepare the data
Fscore_up <- main_Fscore(select='UP') # comment write.table(df_all_all)
Fscore_down <- main_Fscore(select='DOWN') # comment write.table(df_all_all)

Fall <- rbind(Fscore_up$all, Fscore_down$all)
F_all <- reshape2::melt(list(Fall), value.name = "F.score")
F_all <- F_all[, -1]

F_all$Var2 <- factor(F_all$Var2, levels = unique(F_all$Var2))
F_all$Var2 <- relevel(F_all$Var2, ref = "Combat_Wilcox")
F_all$Var2 <- relevel(F_all$Var2, ref = "limma_BEC_Wilcox")
F_all$Var2 <- relevel(F_all$Var2, ref = "MNNCorrect_Wilcox")
F_all$Var2 <- relevel(F_all$Var2, ref = "scMerge_Wilcox")
F_all$Var2 <- relevel(F_all$Var2, ref = "Seurat_Wilcox")
F_all$Var2 <- relevel(F_all$Var2, ref = "ZW_Wilcox")
F_all$Var2 <- relevel(F_all$Var2, ref = "scVI_Wilcox")
F_all$Var2 <- relevel(F_all$Var2, ref = "scGen_Wilcox")
F_all$Var2 <- relevel(F_all$Var2, ref = "Scanorama_Wilcox")
F_all$Var2 <- relevel(F_all$Var2, ref = "Raw_Wilcox")

F_all$Var2 <- relevel(F_all$Var2, ref = "RISC_Wilcox")
F_all$Var2 <- relevel(F_all$Var2, ref = "RISC_QP")

F_all$Var2 <- relevel(F_all$Var2, ref = 'Pseudobulk_DESeq2')
F_all$Var2 <- relevel(F_all$Var2, ref = 'Pseudobulk_edgeR')
F_all$Var2 <- relevel(F_all$Var2, ref = 'Pseudobulk_limma')
F_all$Var2 <- relevel(F_all$Var2, ref = 'Pseudobulk_limmatrend')

F_all$Var2 <- relevel(F_all$Var2, ref = "MAST")
F_all$Var2 <- relevel(F_all$Var2, ref = "MAST_Cov")

F_all$Var2 <- relevel(F_all$Var2, ref = "DESeq2")
F_all$Var2 <- relevel(F_all$Var2, ref = "DESeq2_Cov")
F_all$Var2 <- relevel(F_all$Var2, ref = "ZW_DESeq2")
F_all$Var2 <- relevel(F_all$Var2, ref = "ZW_DESeq2_Cov")

F_all$Var2 <- relevel(F_all$Var2, ref = "edgeR_DetRate")
F_all$Var2 <- relevel(F_all$Var2, ref = "edgeR_DetRate_Cov")
F_all$Var2 <- relevel(F_all$Var2, ref = "edgeR")
F_all$Var2 <- relevel(F_all$Var2, ref = "edgeR_Cov")
F_all$Var2 <- relevel(F_all$Var2, ref = "ZW_edgeR")
F_all$Var2 <- relevel(F_all$Var2, ref = "ZW_edgeR_Cov")

F_all$Var2 <- relevel(F_all$Var2, ref = "limma")
F_all$Var2 <- relevel(F_all$Var2, ref = "limma_Cov")
F_all$Var2 <- relevel(F_all$Var2, ref = "limmatrend")
F_all$Var2 <- relevel(F_all$Var2, ref = "limmatrend_Cov")
F_all$Var2 <- relevel(F_all$Var2, ref = "Combat_limmatrend")
F_all$Var2 <- relevel(F_all$Var2, ref = "MNNCorrect_limmatrend")
F_all$Var2 <- relevel(F_all$Var2, ref = "scMerge_limmatrend")
F_all$Var2 <- relevel(F_all$Var2, ref = 'scVI_limmatrend')
F_all$Var2 <- relevel(F_all$Var2, ref = 'scGen_limmatrend')
F_all$Var2 <- relevel(F_all$Var2, ref = "Scanorama_limmatrend")
F_all$Var2 <- relevel(F_all$Var2, ref = "RISC_limmatrend")

F_all$Var2 <- relevel(F_all$Var2, ref = "DESeq2_FEM")
F_all$Var2 <- relevel(F_all$Var2, ref = "LogNorm_FEM")

F_all$Var2 <- relevel(F_all$Var2, ref = "DESeq2_REM")
F_all$Var2 <- relevel(F_all$Var2, ref = "LogNorm_REM")

F_all$Var2 <- relevel(F_all$Var2, ref = "DESeq2_wFisher")
F_all$Var2 <- relevel(F_all$Var2, ref = "edgeR_wFisher")
F_all$Var2 <- relevel(F_all$Var2, ref = "LogNorm+limmatrend_wFisher")

t = length(unique(F_all$Var2))
meanplot = aggregate(F_all[2],list(rep(1:(nrow(F_all[2])%/%12+1),each=12,len=nrow(F_all))),median)[-1]

rawmedian = rep(c(meanplot$F.score[10], meanplot$F.score[10 + t],
                  meanplot$F.score[10 + 2*t], meanplot$F.score[10 + 3*t]),
                each=nrow(F_all[2])%/%4)
F_all$RawMedian = as.numeric(rawmedian)

library(plyr)
df <- ddply(F_all, .(L1), summarise, median=median(RawMedian, na.rm = TRUE))

RAW <- c('Raw_Wilcox')
BEC <- c("Combat_Wilcox", "limma_BEC_Wilcox", "MNNCorrect_Wilcox", "scMerge_Wilcox", "ZW_Wilcox", "Seurat_Wilcox",
         "Combat_limmatrend", "MNNCorrect_limmatrend", "scMerge_limmatrend", "RISC_limmatrend",
         "RISC_Wilcox", "RISC_NB", "RISC_QP", "scVI_limmatrend", "scGen_limmatrend",        
         "Scanorama_limmatrend", "scVI_Wilcox", "scGen_Wilcox", "Scanorama_Wilcox",
         "Pseudobulk_DESeq2", "Pseudobulk_edgeR", "Pseudobulk_limma", "Pseudobulk_limmatrend"
)

COV <- c("MAST", "MAST_Cov", "DESeq2", "DESeq2_Cov",
         "edgeR", "edgeR_Cov", "edgeR_DetRate", "edgeR_DetRate_Cov",
         "limma", "limma_Cov", "limmatrend","limmatrend_Cov",
         "ZW_edgeR", "ZW_edgeR_Cov", "ZW_DESeq2", "ZW_DESeq2_Cov"  
)

META<- c( "DESeq2_FEM", "voom_FEM", "LogNorm_FEM", "DESeq2_REM", "voom_REM", "LogNorm_REM", 
          "DESeq2_wFisher", "edgeR_wFisher", "LogNorm+limmatrend_wFisher","voom+modt_Fisher", "voom+modt_wFisher"
)

mybreaks <- vect_method

a <- ifelse(mybreaks %in% BEC, 'red4', 
            ifelse(mybreaks %in% COV, 'blue4',
                   ifelse(mybreaks %in% META, 'darkgreen', 
                          ifelse(mybreaks %in% RAW, 'black', 'white'))))

briterhex <- function(colors) {
  res <- c()
  for (i in 1:length(colors)) {
    v <- as.vector(col2rgb(colors[i])) * 1.05
    v <- sapply(v, function(i) {
      min(i, 255)
    })
    res[i] <- rgb(v[1], v[2], v[3], max = 255)
  }
  return(res)
}
Var2_uniq <- rev(unique(F_all$Var2))
Var2_col <- briterhex(scales::hue_pal(h = c(0, 340) + 15,
                                      c = 100, l = 65, h.start = 0,
                                      direction = 1)(length(Var2_uniq)))
names(Var2_col) <- Var2_uniq
Var2_col['Raw_Wilcox'] = 'black'

p4 <- ggplot(F_all, aes(x=Var2, y=F.score, color=Var2)) + scale_color_manual(values = Var2_col)+
  geom_boxplot(outlier.size = 1) + coord_flip() +
  geom_hline(data=df, aes(yintercept=median),linetype="dashed", color='black') + 
  labs(y = 'F-score (beta=0.5)') + theme(axis.text.y = element_text(size=16, color = a))

p4 <- p4 + theme(axis.line = element_line(color = "black", size = 0.5, linetype = "solid"),
                 axis.title.x = element_text(size=22),
                 axis.title.y = element_blank(),
                 axis.text.x = element_text(size=22, colour = 'black', hjust = 1),
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
                 legend.position="none") 


p4 <- p4 + scale_x_discrete(breaks=mybreaks,
                            
                            labels=c(
                              "Combat_Wilcox",
                              "limma_BEC_Wilcox",
                              "MNN_Wilcox",
                              "scMerge_Wilcox",
                              "Seurat_Wilcox",
                              "ZW_BEC_Wilcox",
                              "scVI_Wilcox",
                              "scGen_Wilcox",
                              "Scanorama_Wilcox",
                              "Raw_Wilcox",
                              
                              "RISC_Wilcox",
                              "RISC_QP",
                              
                              'Pseudobulk_DESeq2',
                              'Pseudobulk_edgeR',
                              'Pseudobulk_limma',
                              'Pseudobulk_limmatrend',
                              
                              "MAST",
                              "MAST_Cov",
                              
                              "DESeq2",
                              "DESeq2_Cov",
                              "ZW_DESeq2",
                              "ZW_DESeq2_Cov",
                              
                              "edgeR_DetRate",
                              "edgeR_DetRate_Cov",
                              "edgeR",
                              "edgeR_Cov",
                              "ZW_edgeR",
                              "ZW_edgeR_Cov",
                              
                              "limmavoom",
                              "limmavoom_Cov",
                              "limmatrend",
                              "limmatrend_Cov",
                              "Combat_limmatrend",
                              "MNN_limmatrend",
                              "scMerge_limmatrend",
                              'scVI_limmatrend',
                              'scGen_limmatrend',
                              "Scanorama_limmatrend",
                              "RISC_limmatrend",
                              
                              "DESeq2_FEM",
                              "LogN_FEM",
                              "DESeq2_REM",
                              "LogN_REM",
                              "DESeq2_wFisher",
                              "edgeR_wFisher",
                              "LogN+limmatrend_wFisher"))

p5 <- p4 + geom_point() #geom_jitter(shape=16, position=position_jitter(0.2),size=1)
p5
pdf("sp80_LBE_fbeta.pdf", width=10, height=10)
p5
dev.off()