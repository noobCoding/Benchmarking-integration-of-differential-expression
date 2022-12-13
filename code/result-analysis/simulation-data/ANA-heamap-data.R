rm(list=ls())
heatmap <- c()

### false positive
df1<- read.csv("~/Downloads/Forheatmap_FP_moddepth.txt", sep=' ')
df2<- read.csv("~/Downloads/Forheatmap_FP_lowdepth.txt", sep=' ')

perf_1 <- df1$performance
perf_2 <- df2$performance

perf_1[perf_1=='top_good'] <- 1
perf_2[perf_2=='top_good'] <- 1

perf_1[perf_1=='good'] <- 2
perf_2[perf_2=='good'] <- 2

perf_1[perf_1=='bad'] <- 3
perf_2[perf_2=='bad'] <- 3

perf_1 <- as.numeric(perf_1)
perf_2 <- as.numeric(perf_2)

perf_total <- ceiling((perf_1 + perf_2)/2)
names(perf_total) <- df1$Methods

heatmap <- rbind(heatmap, perf_total)

########## sign preservation  
dtname <- c('pan', 'bcell', 'tcell', 'sp80', '7b', 'avg4_LBE')

library(dplyr)

all_olp <- list()
alop <- c()
for (dt in dtname){
  df <- readRDS(file = paste0(dt, '_deg_sign_check.rds'))
  df <- df[,c('percent','method')]
  
  bpdt <- list()
  for (mt in unique(df$method)){
    # mt = 'Combat_Wilcox'
    df %>% filter(method==mt) %>% select('percent') -> tmp
    bpdt[[mt]] <- tmp$percent
  }
  
  bp<-boxplot(bpdt) 
  bps<-bp$stats
  
  colnames(bps)<-names(bpdt)
  rownames(bps)<-c("Min","First Quartile","Median","Third Quartile","Maximum")
  iqr <- bps['Third Quartile',] - bps['First Quartile',]
  m1 <- abs(bps['Third Quartile',] - bps['First Quartile', 'Raw_Wilcox'])
  m2 <- abs(bps['First Quartile',] - bps['First Quartile', 'Raw_Wilcox'])
  m3 <- abs(bps['First Quartile',] - bps['Third Quartile', 'Raw_Wilcox'])
  m4 <- abs(bps['Third Quartile',] - bps['Third Quartile', 'Raw_Wilcox'])
  
  ovs <- ifelse (m1 > m2, m1, m2)
  ovs <- ifelse (ovs > m3, ovs, m3)
  ovs <- ifelse (ovs > m4, ovs, m4)
  
  dbm <- abs(bps['Median',] - bps['Median', 'Raw_Wilcox'])
  olp <- dbm/ovs*100
  bps <- rbind(bps, iqr, dbm, ovs, olp)
  bps <- round(bps, 2)
  
  all_olp[[dt]] <- bps
  alop <- rbind(alop, bps[c('olp'),])
  
}

top_error <- c()
for (i in 1:46){
  top_error <- c(top_error, mean(alop[,i]))
}
names(top_error) <- colnames(alop)
ranking <- rank(top_error)
ranking[top_error < 30] <- 1
ranking[top_error >= 30 & top_error <=60] <- 2
ranking[top_error > 60] <- 3

#
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
  'Pseudobulk_limmavoom',
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
  "LogN+limmatrend_wFisher")

names(ranking) <- labels
match(names(ranking), colnames(heatmap))

# re-order the columns
ranking<- ranking[colnames(heatmap)]
heatmap <- rbind(heatmap, ranking)

##  speed 
df<- read.csv("~/Downloads/Forheatmap_cpu_time_Covid19+LUAD.txt", sep=' ')
avglogtime <- as.numeric(df$avglogtime)

m1<-min(c(mean(avglogtime), median(avglogtime)))
m2<-max(c(mean(avglogtime), median(avglogtime)))

ranking <- rank(avglogtime)
ranking[avglogtime < m1] <- 1
ranking[avglogtime >= m1 & avglogtime <= m2] <- 2
ranking[avglogtime > m2] <- 3
heatmap <- rbind(heatmap, ranking)

##    scalability
#           cells x genes
# LUAD	    7764	  10278 
# Covid-19	100361	7242
luad <- 7764*10278
covid<- 100361*7242

c<-covid^(1/2)
l<-luad^(1/2)

s <- c(c,l)

###### MSD
pN <- c()
for (m in df$methods){
  time <- as.numeric(df[df$methods==m,][3:4])
  regm <- lm(time~s)
  
  pN <- c(pN, regm$coefficients[2])
  summary(regm)
}
names(pN) <- df$methods
print(pN)

ranking <- rank(pN)
ranking[pN < 1] <- 1
ranking[pN >=1 & pN <=2] <- 2
ranking[pN > 2] <- 3
heatmap <- rbind(heatmap, ranking)

rownames(heatmap) <- c('False_positive/Discovery', 'Sign_preservation', 'Speed', 'Scalability')
saveRDS(heatmap, "heatmapdata.RDS")
