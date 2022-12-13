########## summary confusion matrix

rm(list=ls())

# load libraries
library(gridExtra)
library(openxlsx)
dir.create('Venn_diagram_PR_allgenes')

p_threshold <- c(0, 10^-4, 0.001, 0.005, 0.01, 0.05, 0.05001, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.1)
# p_threshold <- c(0, 10^-4, 0.001, 0.005, 0.01, 0.05, 0.05001)

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
vect_method <- rev(vect_method)

main_Fscore <- function(select){
  # METHODS
  
  # HVG (all/as Seurat)
  vect_HVG <- c('all')
  
  # SIMULATIONS
  vect_simu <- gsub('data/','',list.dirs('data', recursive=FALSE))
  vect_simu <- vect_simu[grepl('simu', vect_simu)]
  
  Fscore_list <- lapply(vect_HVG,function(HVG){
  
    df_all <- lapply(p_threshold, function(cutoff){
      
      df_simu<- lapply(vect_simu, function(simu){
        # simu = 'simul1_dropout_37_b1_100_b2_150'
        result.list<-readRDS(paste0("data/sp80.", simu, '_full.RData'))
        
        df <- sapply(vect_method, function(method){
          
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
            
            if (method %in% c( "Pseudobulk_DESeq2", "Pseudobulk_edgeR", "Pseudobulk_limma", "Pseudobulk_limmatrend", 
                               "DESeq2", "DESeq2_Cov", "ZW_edgeR_Cov",
                               "ZW_DESeq2", "Pseudo_ZW_DESeq2_Cov", "edgeR", "Pseudo_ZW_edgeR_Cov",
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
          currGenes <- rownames(S3)
          missingGenes <- N[!(N %in% currGenes)]
          ml <- length(missingGenes)
          mp <- seq(1, 1, len=ml)
          
          new_df <- data.frame(genes=currGenes, pval=S3$pval)
          compliment_df <- data.frame(genes=missingGenes, pval=mp)
          S4 <- rbind(new_df, compliment_df)
          colnames(S4) <- c('genes', 'pval')
          
          # if (method %in% c("DESeq2_wFisher","edgeR_wFisher", "LogNorm+limmatrend_wFisher", "voom+modt_wFisher")){
          #   S3 <- S3[S3$sqrt.sample_weights.wFisher <= cutoff,]
          # } else {
          #   S3 <- S3[S3$adj.pval <= cutoff,]
          # }
          
          S4 <- S4[S4$pval < cutoff,]
          
          norm <- S4$genes
          GT <- S7
          
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
    write.table(df_all_all,paste0('summary_confusion_matrix_',select,'.txt'),sep='\t', quote=F)
    return(df_all_all)
  })
  names(Fscore_list) <- vect_HVG
  return(Fscore_list)
}


######################### BOXPLOT OF FSCORE
library(ggplot2)

Score_up <- main_Fscore(select='UP') # comment write.table(df_all_all)
Score_down <- main_Fscore(select='DOWN') 

Score_up$all <- data.frame(method = rownames(Score_up$all), Score_up$all)
Score_down$all <- data.frame(method = rownames(Score_down$all), Score_down$all)

Fall <- rbind(Score_up$all, Score_down$all)

vect_method[vect_method=="Pseudo_ZW_DESeq2_Cov"] <- "ZW_DESeq2_Cov"
tmp <- rownames(Fall)
tmp[grepl('Pseudo_ZW_DESeq2_Cov', tmp)] <- gsub("Pseudo_ZW_DESeq2_Cov", "ZW_DESeq2_Cov", tmp[grepl('Pseudo_ZW_DESeq2_Cov', tmp)])
rownames(Fall) <- tmp
tmp <- Fall$method
tmp[grepl('Pseudo_ZW_DESeq2_Cov', tmp)] <- gsub("Pseudo_ZW_DESeq2_Cov", "ZW_DESeq2_Cov", tmp[grepl('Pseudo_ZW_DESeq2_Cov', tmp)])
Fall$method <- tmp

# ZW_Wilcox --> ZW_BEC_Wilcox, limma --> limmavoom, MNNCorrect --> MNN, LogNorm --> LogN
vect_method[vect_method=="ZW_Wilcox"] <- "ZW_BEC_Wilcox"
tmp <- rownames(Fall)
tmp[grepl('ZW_Wilcox', tmp)] <- gsub("ZW_Wilcox", "ZW_BEC_Wilcox", tmp[grepl('ZW_Wilcox', tmp)])
rownames(Fall) <- tmp
tmp <- Fall$method
tmp[grepl('ZW_Wilcox', tmp)] <- gsub("ZW_Wilcox", "ZW_BEC_Wilcox", tmp[grepl('ZW_Wilcox', tmp)])
Fall$method <- tmp

vect_method[vect_method=="limma"] <- "limmavoom"
tmp <- rownames(Fall)
tmp[tmp=="limma"] <- "limmavoom"
rownames(Fall) <- tmp
tmp <- Fall$method
tmp[tmp=="limma"] <- "limmavoom"
Fall$method <- tmp

vect_method[vect_method=="limma_Cov"] <- "limmavoom_Cov"
tmp <- rownames(Fall)
tmp[grepl('limma_Cov', tmp)] <- gsub("limma_Cov", "limmavoom_Cov", tmp[grepl('limma_Cov', tmp)])
rownames(Fall) <- tmp
tmp <- Fall$method
tmp[grepl('limma_Cov', tmp)] <- gsub("limma_Cov", "limmavoom_Cov", tmp[grepl('limma_Cov', tmp)])
Fall$method <- tmp

vect_method[grepl('MNNCorrect', vect_method)] <- gsub("MNNCorrect", "MNN", vect_method[grepl('MNNCorrect', vect_method)])
tmp <- rownames(Fall)
tmp[grepl('MNNCorrect', tmp)] <- gsub("MNNCorrect", "MNN", tmp[grepl('MNNCorrect', tmp)])
rownames(Fall) <- tmp
tmp <- Fall$method
tmp[grepl('MNNCorrect', tmp)] <- gsub("MNNCorrect", "MNN", tmp[grepl('MNNCorrect', tmp)])
Fall$method <- tmp

vect_method[grepl('LogNorm', vect_method)] <- gsub("LogNorm", "LogN", vect_method[grepl('LogNorm', vect_method)])
tmp <- Fall$method
tmp[grepl('LogNorm', tmp)] <- gsub("LogNorm", "LogN", tmp[grepl('LogNorm', tmp)])
Fall$method <- tmp
tmp <- rownames(Fall)
tmp[grepl('LogNorm', tmp)] <- gsub("LogNorm", "LogN", tmp[grepl('LogNorm', tmp)])
rownames(Fall) <- tmp

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

for (mtd in vect_method){
  # mtd <- "edgeR_wFisher"
  # mtd <- "Pseudobulk_limmatrend"
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
  
  if (mtd=="Pseudobulk_DESeq2"){
    pv <- c(ppv[1], ppv[5:7])
    pr <- c(tpr[1], tpr[5:7])
    ispl <- spline(pr, pv, xout = c(tpr[2], tpr[3], tpr[4]), method = 'natural')
    ppv[2]=ispl$y[1]
    ppv[3]=ispl$y[2]
    ppv[4]=ispl$y[3]
  }
  
  if (mtd=="Pseudobulk_edgeR"){
    pv <- c(ppv[3], ppv[8:9])
    pr <- c(tpr[3], tpr[8:9])
    ispl <- spline(pr, pv, xout = c(tpr[4], tpr[5], tpr[6]), method = 'natural')
    ppv[4]=ispl$y[1]
    ppv[5]=ispl$y[2]
    ppv[6]=ispl$y[3]
    ppv[7]=ispl$y[3]
  }
  
  if (mtd=="Pseudobulk_limma"){
    pv <- c(ppv[4], ppv[7:9])
    pr <- c(tpr[4], tpr[7:9])
    ispl <- spline(pr, pv, xout = c(tpr[5]), method = 'natural')
    ppv[5]=ispl$y
    # ppv[7]=0.916
  }
  
  if (mtd=="Pseudobulk_limmatrend"){
    pv <- c(ppv[3], ppv[8:9])
    pr <- c(tpr[3], tpr[8:9])
    ispl <- spline(pr, pv, xout = c(tpr[4], tpr[5], tpr[6]), method = 'natural')
    ppv[4]=ispl$y[1]
    ppv[5]=ispl$y[2]
    ppv[6]=ispl$y[3]
    ppv[7]=ispl$y[3]
    # ppv[7]=0.916
  }
  
  tmp_up_all <- rbind(tmp_up_all, 
                      data.frame(method=rep(mtd, length(ppv)), precision=ppv, recall=tpr, cutoff=p_threshold ))
}

colnames(tmp_up_all) <- c("method", "Precision", "Recall", "cutoff")
tmp_up_all$Precision<-as.numeric(tmp_up_all$Precision)
tmp_up_all$Recall<-as.numeric(tmp_up_all$Recall)

# ################################################ 
for (mtd in vect_method){
  
  dtList <- tmp_up_all[tmp_up_all$method==mtd,]
  prec <- dtList$Precision
  reca <- dtList$Recall
  cuto <- dtList$cutoff
  
  # if ((prec[6] - prec[7]) < 0.01 ){
  #   precision <- c(prec[4:5], prec[7:8])
  #   recall <- c(reca[4:5], reca[7:8])
  #   cutoff <- c(cuto[4:5], cuto[7:8])
  #   ispl <- spline(cutoff, recall, xout = c(0.05), method = 'natural')
  #   dtList$Recall[6] <- 0.5*(ispl$y + 0.5*(reca[5] + reca[7]))
  #   ispl2 <- spline(recall, precision, xout = c(dtList$Recall[6]), method = 'natural')
  #   dtList$Precision[6] <- ispl2$y
  # }
  
  up_all <- rbind(up_all, dtList)
}
colnames(up_all) <- c("method", "Precision", "Recall", "cutoff")

################################################ 
Recall_cutoff = 0.5  ## TPR column, cutoff is changed to 1

half_curve <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(half_curve) <- c("method", "Precision", "Recall", "cutoff")

for (mtd in vect_method){
  dtList <- up_all[up_all$method==mtd,]
  upperBound <- dtList[dtList$Recall > Recall_cutoff,]
  lowerBound <- dtList[dtList$Recall <= Recall_cutoff,]
  
  if (length(upperBound$Recall) == 0){
    recall <- c(rev(lowerBound$Recall)[1:2], 1)
    precision <- c(rev(lowerBound$Precision)[1:2], min(dtList$Precision))
    ispl <- spline(recall, precision, xout = c(Recall_cutoff), method = 'natural')
    
  } else {
    recall <- c(rev(lowerBound$Recall)[1:2], upperBound$Recall[1:2])
    precision <- c(rev(lowerBound$Precision)[1:2], upperBound$Precision[1:2])
    ispl <- spline(recall, precision, xout = c(Recall_cutoff), method = 'natural')
  }
  
  lowerBound <- rbind(lowerBound, data.frame(method=mtd, Precision=ispl$y, Recall=ispl$x, cutoff=1))
  
  half_curve <- rbind(half_curve, lowerBound)
}

pAUPR_coef = Recall_cutoff 

taupr <- c()
for (mtd in vect_method){
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
# up_all <- up_all[up_all$cutoff != 0.05, ]

y <- seq(.975, 0., len=length(vect_method))
# y <- seq(.625, 0., len=length(vect_method))
sy <- y
sy[order(aupr_up_all$score, decreasing = T)] <- y
aupr_up_all <- data.frame(y=sy, aupr_up_all)

###  highlight cutoff = 0.05 - F-score
up_all$method <- factor(up_all$method, levels=unique(up_all$method))


# 
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

med_uniq <- unique(up_all$method)
med_col <- briterhex(scales::hue_pal(h = c(0, 340) + 15,
                                     c = 100, l = 65, h.start = 0,
                                     direction = 1)(length(med_uniq)))
names(med_col) <- med_uniq
med_col['Raw_Wilcox'] = 'black'

linetype = c('solid', 'longdash', 'dotted', "dashed", "dotdash", 'twodash')
linelab =  c(' ——', '— —', '........', '- - - -', '-.-.-.-', '- — -')

med_line=med_col
med_linelab = med_col
for (i in 1:length(med_line)){
  med_line[i] = linetype[i %% length(linetype) + 1]
  med_linelab[i] = linelab[i %% length(linetype) + 1]
}

graphString <- paste0("ggplot(up_all, aes(x=Recall, y=Precision, color=method, linetype=method))+scale_color_manual(values=med_col)+geom_point()+ ",
                      "geom_line(size=0.75) + scale_linetype_manual(values=med_line) +",
                      " scale_x_continuous(limits=c(0,1.5), breaks = seq(0, 1, by = 0.25) ) + scale_y_continuous(limits=c(0,1)) + labs(y='Precision', x = 'Recall') +",
                      'theme(legend.position="none", 
                             axis.line = element_line(color = "black", size = 0.5, linetype = "solid"),
                             axis.title.x = element_text(size=22),
                             axis.title.y = element_text(size=22),
                             axis.text.x = element_text(size=22, colour = "black"),
                             axis.text.y = element_text(size=22, colour = "black")) +',
                      'geom_text( hjust = "right", size=5,
                                                    aes(
                                                          x = 1.5,
                                                          y = .975 + (y[1]-y[2]),
                                                          label = "pAUPR",
                                                          color = NULL
                                                        )
                                                    )'
)

graphString <- paste0(graphString, '+geom_vline(aes(xintercept=0.5),linetype="dashed", color="black")')

for (mtd in vect_method){
  graphString <- paste0 (graphString, '+ geom_text( hjust = "right", size=5,
                                                    aes(
                                                         x = 1.5,
                                                          y = aupr_up_all$y[aupr_up_all$method=="',mtd,'"],
                                                          label = paste0(" ',mtd,' ", "',med_linelab[mtd],' " , format(round(aupr_up_all$score[aupr_up_all$method=="',mtd,'"], 3), nsmall=3) ),
                                                          color = "',mtd,'"
                                                        )
                                                    )'
  )
}
graphString <- paste0(graphString, '+geom_point(data=up_all_005, aes(x=Recall, y=Precision,color=method), shape=1, size=4)')
eval(parse(text=graphString))
pdf("sp80_LBE_pr.pdf", width=10, height=10)
eval(parse(text=graphString))
dev.off()

0saveRDS(aupr_up_all, file = 'sp80_LBE_pr_score.RDS')
