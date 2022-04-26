bec<-readRDS(file = '4b_bec_deg_swap.rds')
cov<-readRDS(file = '4b_cov_deg_swap.rds')
meta<-readRDS(file = '4b_meta_deg_sign_swap_new.rds')

# df <- rbind(data.frame(bec), data.frame(cov), data.frame(meta))
# # View(df)
# df <- df[, c(1, 7, 13, 19, 25, 31)]
# df <- df[, grepl('percent', colnames(df), fixed = TRUE)]
bec <- data.frame(bec)
cov <- data.frame(cov)
meta <- data.frame(meta)

# # ## all
# bec <- bec[, c(1, 7, 13, 19, 25, 31)]
# cov <- cov[, c(1, 7, 13, 19, 25, 31)]
# meta <- meta[, c(1, 7, 13, 19, 25, 31)]

# 5%
bec <- bec[, c(2, 8, 14, 20, 26, 32)]
cov <- cov[, c(2, 8, 14, 20, 26, 32)]
meta <- meta[, c(2, 8, 14, 20, 26, 32)]

colnames(bec) <- c(1, 2, 3, 4, 5, 6)
colnames(cov) <- c(1, 2, 3, 4, 5, 6)
colnames(meta) <- c(1, 2, 3, 4, 5, 6)
df <- rbind(bec, cov, meta)

tmp <- rownames(df)
tmp[tmp=='limma'] = 'limma_bec'
rownames(df) <- tmp

#######################################################################################
allpval <- reshape2::melt(unlist(df), value.name = "percent_pval")
# View(allpval)
allpval$method <- rep(tmp, 6)
allpval$method <- factor(allpval$method, levels = unique(allpval$method))

allpval$method <- relevel(allpval$method, ref = "Combat")
allpval$method <- relevel(allpval$method, ref = "limma_bec")
allpval$method <- relevel(allpval$method, ref = "MNN")
allpval$method <- relevel(allpval$method, ref = "scmerge")
allpval$method <- relevel(allpval$method, ref = "seurat3")
allpval$method <- relevel(allpval$method, ref = "zinbwave")
allpval$method <- relevel(allpval$method, ref = "raw")

allpval$method <- relevel(allpval$method, ref = "mast")
allpval$method <- relevel(allpval$method, ref = "mast_cov")

allpval$method <- relevel(allpval$method, ref = "deseq2")
allpval$method <- relevel(allpval$method, ref = "deseq2_cov")
allpval$method <- relevel(allpval$method, ref = "zinbwave_deseq2")
allpval$method <- relevel(allpval$method, ref = "zinbwave_deseq2_cov")
allpval$method <- relevel(allpval$method, ref = "deseq2_Fisher")
allpval$method <- relevel(allpval$method, ref = "deseq2_wFisher_sample+gene")

allpval$method <- relevel(allpval$method, ref = "edger_detrate")
allpval$method <- relevel(allpval$method, ref = "edger_detrate_cov")
allpval$method <- relevel(allpval$method, ref = "edger")
allpval$method <- relevel(allpval$method, ref = "edger_cov")
allpval$method <- relevel(allpval$method, ref = "zinbwave_edger")
allpval$method <- relevel(allpval$method, ref = "zinbwave_edger_cov")
allpval$method <- relevel(allpval$method, ref = "edgeR_Fisher")
allpval$method <- relevel(allpval$method, ref = "edgeR_wFisher_sample+gene")


allpval$method <- relevel(allpval$method, ref = "limma_voom")
allpval$method <- relevel(allpval$method, ref = "limma_voom_cov")
allpval$method <- relevel(allpval$method, ref = "limma_trend_false")
allpval$method <- relevel(allpval$method, ref = "limma_trend_false_cov")
allpval$method <- relevel(allpval$method, ref = "combat_limma_trend_false")
allpval$method <- relevel(allpval$method, ref = "mnn_limma_trend_false")
allpval$method <- relevel(allpval$method, ref = "scmerge_limma_trend_false")
allpval$method <- relevel(allpval$method, ref = "LogNormalize+limma.trend_Fisher")
allpval$method <- relevel(allpval$method, ref = "LogNormalize+limma.trend_wFisher_sample+gene")


allpval$method <- relevel(allpval$method, ref = "deseq2+FEM")
allpval$method <- relevel(allpval$method, ref = "voom+FEM")
allpval$method <- relevel(allpval$method, ref = "LogNormalize+FEM")

allpval$method <- relevel(allpval$method, ref = "deseq2+REM")
allpval$method <- relevel(allpval$method, ref = "voom+REM")
allpval$method <- relevel(allpval$method, ref = "LogNormalize+REM")

allpval$method <- relevel(allpval$method, ref = "voom+modt_Fisher")
allpval$method <- relevel(allpval$method, ref = "voom+modt_wFisher_sample+gene")


mybreaks <- c("Combat", "limma_bec",  "MNN",  "scmerge",  "seurat3",  "zinbwave",  "raw",
              "mast",  "mast_cov",  "deseq2",  "deseq2_cov",  "zinbwave_deseq2",  "zinbwave_deseq2_cov",
              "deseq2_Fisher",  "deseq2_wFisher_sample+gene",
              "edger_detrate",  "edger_detrate_cov",  "edger",  "edger_cov",  "zinbwave_edger",  "zinbwave_edger_cov",
              "edgeR_Fisher",  "edgeR_wFisher_sample+gene",
              "limma_voom",  "limma_voom_cov",  "limma_trend_false",  "limma_trend_false_cov",
              "combat_limma_trend_false",  "mnn_limma_trend_false",  "scmerge_limma_trend_false",
              "LogNormalize+limma.trend_Fisher",  "LogNormalize+limma.trend_wFisher_sample+gene",
              
              "deseq2+FEM",   "LogNormalize+FEM",  "voom+FEM",
              "deseq2+REM",    "LogNormalize+REM", "voom+REM",
              "voom+modt_Fisher",  "voom+modt_wFisher_sample+gene")

RAW <- c('raw')
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
a <- ifelse(mybreaks %in% BEC, 'red4',
            ifelse(mybreaks %in% COV, 'blue4',
                   ifelse(mybreaks %in% META, 'darkgreen',
                          ifelse(mybreaks %in% RAW, 'black', 'white'))))

# briterhex <- function(colors) {
#   res <- c()
#   for (i in 1:length(colors)) {
#     v <- as.vector(col2rgb(colors[i])) * 1.05
#     v <- sapply(v, function(i) {
#       min(i, 255)
#     })
#     res[i] <- rgb(v[1], v[2], v[3], max = 255)
#   }
#   return(res)
# }
# 
method_uniq <- rev(unique(allpval$method))
method_col <- briterhex(scales::hue_pal(h = c(0, 360) + 15,
                                        c = 100, l = 65, h.start = 0,
                                        direction = 1)(length(method_uniq)))
names(method_col) <- method_uniq
method_col['raw'] = 'black'

# View(allpval)
t = length(unique(allpval$method))
medianraw = median(allpval$percent_pval[allpval$method=='raw'])
rawmedian = rep(median(medianraw), length(allpval$method))
allpval$RawMedian = as.numeric(rawmedian)

library(plyr)
df <- ddply(allpval, .(percent_pval), summarise, median=median(RawMedian, na.rm = TRUE))

library(ggplot2)
p <- ggplot(allpval, aes(x=method, y=percent_pval, color=method)) + geom_boxplot(outlier.size = 1, width=0.6) + geom_point(size=2)
# p <- p + scale_color_manual(values = method_col)
p <- p + labs(x= NULL, y= "Percentage") 
p <- p + theme(legend.position=0) + coord_flip()
p <- p + theme(axis.title.y = element_text(size = 15, vjust= 0.5))
p <- p + theme(axis.text = element_text(size = 12))
p <- p + geom_hline(data=df, aes(yintercept=median),linetype="dashed", color='black')
p <- p + theme(axis.text.y = element_text(size=14, color = a))


p <- p + theme(axis.line = element_line(color = "black", size = 0.5, linetype = "solid"),
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


p <- p + scale_x_discrete(breaks=mybreaks,
                          
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

p