# rm(list=ls())

pan<-readRDS(file = 'pan_deg_swap_bps.rds')
bcell<-readRDS(file = 'bcell_deg_swap_bps.rds')
tcell<-readRDS(file = 'tcell_deg_swap_bps.rds')
sp2b<-readRDS(file = '2b_deg_swap_bps.rds')
sp4b<-readRDS(file = '4b_deg_swap_bps.rds')

olpall <- data.frame(sp2b=sp2b['olp',],
                     sp4b=sp4b['olp',],
                     pan=pan['olp',],
                     # bcell=bcell['olp',],
                     tcell=tcell['olp',])

######################################################
df <- t(olpall)
# boxplot(df, horizontal=TRUE, )
bp<-boxplot(df, xaxt = "n", yaxt = "n") 
abline(h=33, col = "blue", lty = 5, lwd = 3)
# abline(h=49.5, col = "green", lty = 5, lwd = 3)
abline(h=66, col = "red", lty = 5, lwd = 3)

## Draw x-axis without labels.
axis(side = 2, labels = FALSE)

## Draw y-axis.
axis(side = 2,
     ## Rotate labels perpendicular to y-axis.
     las = 2,
     ## Adjust y-axis label positions.
     mgp = c(3, 0.75, 0))

## Draw the x-axis labels.
text(x = 1:length(colnames(df)),
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.45,
     ## Use names from the data list.
     labels = colnames(df),
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 45,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.965,
     ## Increase label size.
     cex = 1.2)

bps<-bp$stats

colnames(bps)<-colnames(sp2b)
rownames(bps)<-c("Min","First Quartile","Median","Third Quartile","Maximum")
iqr <- bps['Third Quartile',] - bps['First Quartile',]
ovs1 <- abs(bps['Third Quartile',] - bps['First Quartile', 'raw'])
ovs2 <- abs(bps['First Quartile',] - bps['Third Quartile', 'raw'])
ovs <- ifelse(ovs1 > ovs2, ovs1, ovs2)

dbm <- bps['Median',] - bps['Median', 'raw']
olp <- dbm/ovs*100

bps <- rbind(bps, iqr, dbm, ovs)
bps <- rbind(bps, olp)
bps <- round(bps, 2)

bps <- t(bps)
bps[is.na(bps)]=0
sum(bps[,'dbm']<=33)
sum(bps[,'dbm']<=66) - sum(bps['dbm']<=33)
sum(bps[,'dbm']>66)
# bps <- t(bps)

# rownames(bps)[bps[,'dbm']<=33]
# rownames(bps)[bps[,'dbm']>33 & bps[,'dbm']<=66]

#########################
# install.packages("factoextra")
library(factoextra)

mylabels=c("Raw_Wilcox", "Combat_Wilcox", "limma_BEC_Wilcox","MNNCorrect_Wilcox",
         "Seurat_Wilcox","scMerge_Wilcox", "ZINB-WaVE_Wilcox",
         
         'MAST', 'MAST_Cov',
         'DESeq2','DESeq2_Cov', "ZINB-WaVE_DESeq2", "ZINB-WaVE_DESeq2_Cov",
         
         'edgeR_DetRate', 'edgeR_DetRate_Cov', 
         'edgeR', 'edgeR_Cov', "ZINB-WaVE_edgeR",  "ZINB-WaVE_edgeR_Cov",
         
         "limma", 'limma_Cov', 'limmatrend','limmatrend_Cov',
         'Combat_limmatrend', 'MNNCorrect_limmatrend', 'scMerge_limmatrend',
         
         "voom+REM", "voom+FEM",
         "voom+modt_wFisher", "edgeR_wFisher", "DESeq2_wFisher",
         
         "LogNorm+FEM", "LogNorm+REM", "LogNorm+limmatrend_wFisher",  
         
         "voom+modt_Fisher", "edgeR_Fisher", "DESeq2_Fisher", "LogNorm+limmatrend_Fisher", 
         
         "DESeq2+FEM", "DESeq2+REM"
          )

# Loading dataset
# df <- olpall
# df <- bps[,c('dbm','ovs')]
# df <- bps[,c('dbm')]
df <- bps[,!colnames(bps) %in% c('olp')]
df <- bps[,!colnames(bps) %in% c('olp', 'dbm', 'ovs')]
rownames(df) <- mylabels

# Omitting any NA values
# df <- na.omit(df)

# Scaling dataset
df <- scale(df)

# output to be present as PNG file
# png(file = "KMeansExample.png")

km <- kmeans(df, centers = 3, nstart = 19)
# km$cluster[km$cluster==2]

# Visualize the clusters
print(fviz_cluster(km, data = df))
# # saving the file
# dev.off()



# # Installing Packages
# install.packages("fpc")
# # Loading package
# library(fpc)
# 
# #Loading data
# data(iris)
# # Structure 
# str(iris)
# # Fitting DBScan clustering Model 
# # to training dataset
# iris_1 <- iris[-5]
# 
# # Fitting DBScan clustering Model 
# # to training dataset
# set.seed(220)  # Setting seed
# Dbscan_cl <- dbscan(iris_1, eps = 0.45, MinPts = 5)
# Dbscan_cl
# 
# # Checking cluster
# Dbscan_cl$cluster
# 
# # Table
# table(Dbscan_cl$cluster, iris$Species)
# 
# # Plotting Cluster
# plot(Dbscan_cl, iris_1, main = "DBScan")
# plot(Dbscan_cl, iris_1, main = "Petal Width vs Sepal Length")