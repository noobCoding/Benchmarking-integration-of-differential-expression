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

all_simu <- list()
Fall <- data.frame()
vect_simu <- gsub('data/','',list.dirs('data', recursive=FALSE))

for (simu in vect_simu){
  DEG <- read.table(paste0('data/',simu,'/de__genes.txt'), head=T, fill = T)
  up <- read.table(paste0('data/',simu,'/true_up_genes.txt'), head=T, fill = T)
  down <- read.table(paste0('data/',simu,'/true_down_genes.txt'), head=T, fill = T)
  gene <- read.table(paste0('data/',simu,'/geneinfo.txt'), head=T, fill = T)
  
  counts <- read.table(file = paste0('data/',simu, '/', 'counts.txt'),header=T,row.names=1,check.names = F)
  cellinfo <- read.table(paste0("data/",simu,"/cellinfo.txt"),sep="\t",header=T,row.names=1)
  colnames(counts) <- rownames(cellinfo)
  #### 
  
  batch1 <- rownames(cellinfo[cellinfo$Batch=="Batch1",])
  batch2 <- rownames(cellinfo[cellinfo$Batch=="Batch2",])
  
  group1 <- rownames(cellinfo[cellinfo$Group=="Group1",])
  group2 <- rownames(cellinfo[cellinfo$Group=="Group2",])
  
  # batch1 <- rownames(cellinfo[cellinfo$Batch=="1",])
  # batch2 <- rownames(cellinfo[cellinfo$Batch=="2",])
  #
  # group1 <- rownames(cellinfo[cellinfo$Group=="A",])
  # group2 <- rownames(cellinfo[cellinfo$Group=="B",])
  
  #batch1 2
  # cog_1 <- counts[, group1]
  # cog_2 <- counts[, group2]
  # 
  # lib_1 <- colSums(cog_1)
  # lib_2 <- colSums(cog_2)
  # 
  # genes1 <- c()
  # for (g in DEG$x){
  #   tg_1 <- cog_1[rownames(cog_1)==g, ]
  #   tg_1 <- tg_1/lib_1
  #   tg_1 <- tg_1[tg_1 != 0]
  #   tg_1 <- log2(1 +  tg_1 * 10^4)
  #   g1 <- sum(tg_1)/length(tg_1)
  #   
  #   tg_2 <- cog_2[rownames(cog_2)==g, ]
  #   tg_2 <- tg_2/lib_2
  #   tg_2 <- tg_2[tg_2 != 0]
  #   tg_2 <- log2(1 +  tg_2 * 10^4)
  #   g2 <- sum(tg_2)/length(tg_2)
  #   
  #   genes1 <- c(genes1, g1 - g2)
  # }
  # names(genes1) <- DEG$x
  # logFC_raw <- genes1
  
  # simu = 'simul5_dropout_37'
  result.list<-readRDS(paste0("data/sp80.", simu, '_full.RData'))
  
  hist_005 <- c()
  perc_005 <- c()
  hist <- c()
  perc <- c()
  
  hist_qval <- c()
  perc_qval <- c()
  
  hist_both <- c()
  perc_both <- c()
  
  for (method in vect_method){
    print(method)
    
    count_qval = 0
    count_logfc = 0
    count_both = 0
    count = 0
    
    significant_logfc = 0
    significant_qval = 0
    significant_both = 0
    
    cutoff = 0.05
    
    # method = "scVI_Wilcox"
    raw <- result.list[['Raw_Wilcox']]
    S3 <- result.list[[method]]
    tmp <- names(S3)
    tmp[ which(tmp=="logFC") ] <- 'log2FC'
    names(S3) <- tmp
    if (method %in% c( "Pseudobulk_DESeq2", "Pseudobulk_edgeR", "Pseudobulk_limma",
                       "Pseudobulk_limmatrend", "DESeq2", "DESeq2_Cov",
                       "ZW_DESeq2", "Pseudo_ZW_DESeq2_Cov", "edgeR",
                       "edgeR_Cov", "edgeR_DetRate", "edgeR_DetRate_Cov", "ZW_edgeR", "ZW_edgeR_Cov",            
                       "limma", "limma_Cov", "limmatrend",        
                       "limmatrend_Cov", "Combat_limmatrend", "MNNCorrect_limmatrend", "scMerge_limmatrend",  
                       "RISC_limmatrend", "RISC_NB", "RISC_QP", "scVI_limmatrend", "scGen_limmatrend",        
                       "Scanorama_limmatrend"))
    {
    } else if (method %in% c("DESeq2_FEM", "voom_FEM", "LogNorm_FEM", "DESeq2_REM", "voom_REM", "LogNorm_REM")){
     
    } else if (method %in% c("DESeq2_wFisher","edgeR_wFisher", "LogNorm+limmatrend_wFisher")){
     
    } else if (method %in% c('Scanorama_Wilcox', 'scVI_Wilcox', 'scGen_Wilcox' )){
      S3$log2FC <- S3$log2FC * -1.
    }
    
    
    for (g in up$x){
      if ( S3$log2FC[rownames(S3)==g] > 0){ 
        
        if (S3$log2FC[rownames(S3)==g] > 0.5 ){ # logFC > 0.5
          significant_logfc = significant_logfc + 1
        }
        
        if (S3$adj.pval[rownames(S3)==g] <= 0.05 ){ # q-value cutoff 0.05
          significant_qval = significant_qval + 1
        }
        
        if (S3$adj.pval[rownames(S3)==g] <= 0.05 ){ # q-value cutoff 0.05
          if (S3$log2FC[rownames(S3)==g] > 0.5 ){ # logFC > 0.5
            significant_both = significant_both + 1
          }
        }
        
        count = count + 1
      }
    }
    
    for (g in down$x){
      if (S3$log2FC[rownames(S3)==g] < 0){
        if (S3$log2FC[rownames(S3)==g] < -0.5 ){ # logFC > 0.5
          significant_logfc = significant_logfc + 1
        }
        
        if (S3$adj.pval[rownames(S3)==g] <= 0.05 ){ # q-value cutoff 0.05
          significant_qval = significant_qval + 1
        }
        
        if (S3$adj.pval[rownames(S3)==g] <= 0.05 ){ # q-value cutoff 0.05
          if (S3$log2FC[rownames(S3)==g] < -0.5 ){ # logFC > 0.5
            significant_both = significant_both + 1
          }
        }
        
        count = count + 1
      }
    }
    
    perc_005 <- c(perc_005, round(100*significant_logfc/length(DEG$x), 2) )
    perc_qval <- c(perc_qval, round(100*significant_qval/length(DEG$x), 2) )
    perc_both <- c(perc_both, round(100*significant_both/length(DEG$x), 2) )
    
    hist <- c(hist, count)
    perc <- c(perc, round(100*count/length(DEG$x), 2) )
    
    }
  
  df <- data.frame(percent=perc, percent_005=perc_005,
                   percent_qval=perc_qval,
                   perc_both = perc_both, count=hist, 
                   simu=simu, method=method)
  
  rownames(df) <- vect_method
  df$method <- vect_method
  all_simu[[simu]] <- df
  Fall <- rbind(Fall, df)
}
saveRDS(Fall, file = paste0('sp80_deg_sign_check.rds'))

#     
######################### BOXPLOT OF FSCORE
library(ggplot2)
library(cowplot)

F_all <- Fall

F_all$method <- factor(F_all$method, levels = unique(F_all$method))

F_all$method <- relevel(F_all$method, ref = "Combat_Wilcox")
F_all$method <- relevel(F_all$method, ref = "limma_BEC_Wilcox")
F_all$method <- relevel(F_all$method, ref = "MNNCorrect_Wilcox")
F_all$method <- relevel(F_all$method, ref = "scMerge_Wilcox")
F_all$method <- relevel(F_all$method, ref = "Seurat_Wilcox")
F_all$method <- relevel(F_all$method, ref = "ZW_Wilcox")
F_all$method <- relevel(F_all$method, ref = "scVI_Wilcox")
F_all$method <- relevel(F_all$method, ref = "scGen_Wilcox")
F_all$method <- relevel(F_all$method, ref = "Scanorama_Wilcox")
F_all$method <- relevel(F_all$method, ref = "Raw_Wilcox")

F_all$method <- relevel(F_all$method, ref = "RISC_Wilcox")
F_all$method <- relevel(F_all$method, ref = "RISC_QP")

F_all$method <- relevel(F_all$method, ref = 'Pseudobulk_DESeq2')
F_all$method <- relevel(F_all$method, ref = 'Pseudobulk_edgeR')
F_all$method <- relevel(F_all$method, ref = 'Pseudobulk_limma')
F_all$method <- relevel(F_all$method, ref = 'Pseudobulk_limmatrend')

F_all$method <- relevel(F_all$method, ref = "MAST")
F_all$method <- relevel(F_all$method, ref = "MAST_Cov")

F_all$method <- relevel(F_all$method, ref = "DESeq2")
F_all$method <- relevel(F_all$method, ref = "DESeq2_Cov")
F_all$method <- relevel(F_all$method, ref = "ZW_DESeq2")
F_all$method <- relevel(F_all$method, ref = "ZW_DESeq2_Cov")

F_all$method <- relevel(F_all$method, ref = "edgeR_DetRate")
F_all$method <- relevel(F_all$method, ref = "edgeR_DetRate_Cov")
F_all$method <- relevel(F_all$method, ref = "edgeR")
F_all$method <- relevel(F_all$method, ref = "edgeR_Cov")
F_all$method <- relevel(F_all$method, ref = "ZW_edgeR")
F_all$method <- relevel(F_all$method, ref = "ZW_edgeR_Cov")

F_all$method <- relevel(F_all$method, ref = "limma")
F_all$method <- relevel(F_all$method, ref = "limma_Cov")
F_all$method <- relevel(F_all$method, ref = "limmatrend")
F_all$method <- relevel(F_all$method, ref = "limmatrend_Cov")
F_all$method <- relevel(F_all$method, ref = "Combat_limmatrend")
F_all$method <- relevel(F_all$method, ref = "MNNCorrect_limmatrend")
F_all$method <- relevel(F_all$method, ref = "scMerge_limmatrend")
F_all$method <- relevel(F_all$method, ref = 'scVI_limmatrend')
F_all$method <- relevel(F_all$method, ref = 'scGen_limmatrend')
F_all$method <- relevel(F_all$method, ref = "Scanorama_limmatrend")
F_all$method <- relevel(F_all$method, ref = "RISC_limmatrend")

F_all$method <- relevel(F_all$method, ref = "DESeq2_FEM")
F_all$method <- relevel(F_all$method, ref = "LogNorm_FEM")

F_all$method <- relevel(F_all$method, ref = "DESeq2_REM")
F_all$method <- relevel(F_all$method, ref = "LogNorm_REM")

F_all$method <- relevel(F_all$method, ref = "DESeq2_wFisher")
F_all$method <- relevel(F_all$method, ref = "edgeR_wFisher")
F_all$method <- relevel(F_all$method, ref = "LogNorm+limmatrend_wFisher")

library(plyr)

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
method_uniq <- rev(unique(F_all$method))
method_col <- briterhex(scales::hue_pal(h = c(0, 340) + 15,
                                      c = 100, l = 65, h.start = 0,
                                      direction = 1)(length(method_uniq)))
names(method_col) <- method_uniq
method_col['Raw_Wilcox'] = 'black'

t = length(unique(F_all$method))
rawmedian = median(F_all$perc_both[F_all$method=='Raw_Wilcox'])
rawmedian = rep(median(rawmedian), length(F_all$method))
F_all$RawMedian = as.numeric(rawmedian)
df <- ddply(F_all, .(perc_both), summarise, median=median(RawMedian, na.rm = TRUE))
p4 <- ggplot(F_all, aes(x=method, y=perc_both, color=method)) + scale_color_manual(values = method_col)+
  geom_boxplot(outlier.size = 1) + coord_flip() + ylim(c(0, 2))+
  geom_hline(data=df, aes(yintercept=median),linetype="dashed", color='black') +
  # labs(y = 'F-score') + theme(axis.text.y = element_text(size=14, color = a))
  labs(y = 'Error Ratio(%): abs(logFC)>0.5 and q-value≤0.05') + theme(axis.text.y = element_text(size=16, color = a))

# t = length(unique(F_all$method))
# rawmedian = median(F_all$percent_005[F_all$method=='Raw_Wilcox'])
# rawmedian = rep(median(rawmedian), length(F_all$method))
# F_all$RawMedian = as.numeric(rawmedian)
# df <- ddply(F_all, .(percent_005), summarise, median=median(RawMedian, na.rm = TRUE))

# t = length(unique(F_all$method))
# rawmedian = median(F_all$percent_qval[F_all$method=='Raw_Wilcox'])
# rawmedian = rep(median(rawmedian), length(F_all$method))
# F_all$RawMedian = as.numeric(rawmedian)
# df <- ddply(F_all, .(percent_qval), summarise, median=median(RawMedian, na.rm = TRUE))
# p4 <- ggplot(F_all, aes(x=method, y=percent_qval, color=method)) + scale_color_manual(values = method_col)+
# 
# # p4 <- ggplot(F_all, aes(x=method, y=percent_005, color=method)) + scale_color_manual(values = method_col)+
#   geom_boxplot(outlier.size = 1) + coord_flip() + #ylim(c(0, 2))+
#   geom_hline(data=df, aes(yintercept=median),linetype="dashed", color='black') +
#   # labs(y = 'F-score') + theme(axis.text.y = element_text(size=14, color = a))
#   labs(y = 'Error Ratio (%)') + theme(axis.text.y = element_text(size=16, color = a))

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
                 # strip.background.x = element_blank(),
                 # strip.background.y = element_blank(),
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
print(p5)
pdf("sp80_deg_signs.pdf", width=10, height=10)
p5
dev.off()