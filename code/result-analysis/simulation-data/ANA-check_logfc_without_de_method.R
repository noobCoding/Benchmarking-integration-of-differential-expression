rm(list=ls())

# load libraries
library(gridExtra)
library(openxlsx)

vect_method <- c('combat', 'limma_bec', 'mnn_opt', 'Seurat', 'scMerge',  'zinbwave', 'scvi', 'scgen', 'scanorama', 'risc')
name_method <- c('Combat', 'limma_BEC', 'MNN', 'Seurat', 'scMerge',  'ZINBWaVE_BEC', 'scVI', 'scGen', 'Scanorama', 'RISC')

mcolor <- c(
  "#009E73", "gold1", "#CC79A7", "dodgerblue2", "#FF7F00", "#FB9A99", "#0072B2", "green4", "#E31A1C", "#6A3D9A",
  "skyblue2", "#FDBF6F", 
  "palegreen2",
  "maroon", "orchid1", "deeppink1",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "#E69F00",  "#56B4E9"
)

all_simu <- list()
Fall <- data.frame()

vect_simu <- gsub('data/','',list.dirs('data', recursive=FALSE))
vect_simu <- vect_simu[c(1, 2, 3, 4, 5, 6)]


all_logfc <- list()

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

    #batch1 2
  cog_1 <- counts[, group1]
  cog_2 <- counts[, group2]
  
  lib_1 <- colSums(cog_1)
  lib_2 <- colSums(cog_2)
  
  genes1 <- c()
  for (g in DEG$x){
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
  names(genes1) <- DEG$x
  logFC_raw <- genes1
  # all_logfc$raw <- logFC_raw
  
  df <- data.frame()
  perc_005 <- c()
  perc <- c()
  cosign <- c()
  
  for (method in vect_method){
    # import data, sample
    if (method %in% c('scvi',  'scgen','scanorama')){
      corrected_data <- read.csv(file = paste0('data/',simu, '/', method , '_corrected_data.csv'),header=T,row.names=1,check.names = F)
      if (method == 'scvi'){
      } else
        if (method == 'scanorama'){
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
        corrected_data <- first_processed$count * 20
        
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
        for (g in DEG$x){
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
        for (g in DEG$x){
          tg_1 <- as.numeric(cog_1[rownames(cog_1)==g, ])
          g1 <- log(x=mean(x=expm1(x=tg_1)) + 1)
          
          tg_2 <- as.numeric(cog_2[rownames(cog_2)==g, ])
          g2 <- log(mean(expm1(tg_2)) + 1)
          
          logFC_corrected <- c(logFC_corrected, g1 - g2)
        }
      }
      names(logFC_corrected) <- DEG$x
      
      logFC_corrected[is.na(logFC_corrected)] <- 0
      all_logfc[[method]] <- logFC_corrected
      
      perc <- c(perc, 100*round(sum(logFC_corrected*logFC_raw<0)/length(logFC_raw), 2) )
      
      sig_deg <- names(logFC_raw[abs(logFC_raw) > 0.5 ])
      logfc_sig_raw <- logFC_raw[abs(logFC_raw) > 0.5 ]
      logfc_sig_cor <- logFC_corrected[names(logFC_corrected) %in% sig_deg]
      
      perc_005 <- c(perc_005, 100*round(sum(logfc_sig_cor*logfc_sig_raw<0)/length(logFC_raw), 2) )
      
      cos_dist = 0
      for (i in 1:length(logFC_raw)){
        cos_dist = cos_dist + 1 - sign(logFC_raw[i])* c(logFC_raw[i], logFC_corrected[i]) %*% c(1,1) / 
          ( sqrt(2)*sqrt(logFC_raw[i]^2 + logFC_corrected[i]^2) )
      }
      cos_dist = cos_dist/length(logFC_raw)
      
      cosign <- c(cosign, cos_dist)
  }
  
  df <- data.frame(percent=perc, 
                   percent_005=perc_005, 
                   cosign = cosign,
                   simu=simu, 
                   method=method)
  
  rownames(df) <- vect_method
  df$method <- vect_method
  all_simu[[simu]] <- df
  Fall <- rbind(Fall, df)
}
saveRDS(Fall, file = paste0('sp80_deg_sign_check_b4.rds'))

#     
######################### BOXPLOT OF FSCORE
library(ggplot2)
library(cowplot)

F_all <- Fall

F_all$method <- factor(F_all$method, levels = unique(F_all$method))

F_all$method <- relevel(F_all$method, ref = "combat")
F_all$method <- relevel(F_all$method, ref = "limma_bec")
F_all$method <- relevel(F_all$method, ref = "mnn_opt")
F_all$method <- relevel(F_all$method, ref = "Seurat")
F_all$method <- relevel(F_all$method, ref = "scMerge")
F_all$method <- relevel(F_all$method, ref = "zinbwave")
F_all$method <- relevel(F_all$method, ref = "scvi")
F_all$method <- relevel(F_all$method, ref = "scgen")
F_all$method <- relevel(F_all$method, ref = "scanorama")
F_all$method <- relevel(F_all$method, ref = "risc")

mybreaks <- vect_method
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

p4 <- ggplot(F_all, aes(x=method, y=cosign, color=method)) + scale_color_manual(values = method_col)+
  geom_boxplot(outlier.size = 1) + coord_flip() + ylim(0, 0.6)+
  labs(y = 'Average Cosign Distance') + theme(axis.text.y = element_text(size=16, color = 'black'))

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
                            labels=name_method)

p5 <- p4 + geom_point() #geom_jitter(shape=16, position=position_jitter(0.2),size=1)
p5
pdf("sp80_deg_signs_b4.pdf", width=10, height=10)
p5
dev.off()