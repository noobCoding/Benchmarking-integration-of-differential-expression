library(magrittr)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(cluster)
library(factoextra)
library(dendextend)
methods_to_plot<-c('Pseudobulk_DESeq2','Pseudobulk_edgeR','Pseudobulk_limmavoom','Pseudobulk_limmatrend','Combat_Wilcox','limma_BEC_Wilcox','MNN_Wilcox','scMerge_Wilcox','Seurat_Wilcox','ZW_BEC_Wilcox','Raw_Wilcox','RISC_Wilcox' , "Scanorama_Wilcox", "scGen_Wilcox", "scVI_Wilcox",'MAST','MAST_Cov','DESeq2','DESeq2_Cov','ZW_DESeq2','ZW_DESeq2_Cov','DESeq2_wFisher','edgeR_DetRate','edgeR_DetRate_Cov','edgeR','edgeR_Cov','ZW_edgeR','ZW_edgeR_Cov','edgeR_wFisher','limmavoom','limmavoom_Cov','limmatrend','limmatrend_Cov','Combat_limmatrend','MNN_limmatrend','scMerge_limmatrend','LogN+limmatrend_wFisher','RISC_limmatrend',"Scanorama_limmatrend", "scGen_limmatrend", "scVI_limmatrend",'DESeq2_FEM','LogN_FEM','DESeq2_REM','LogN_REM','RISC_QP')
{
  # analyzed_data.dir=c('scLUAD_Epithelial','scLUAD_Myeloid','scLUAD_T-NK','COVID-19_Monocyte')
  # select two data for comparison. Run analysis for first data and then run same for the second data.
  analyzed_data.dir='scLUAD_Epithelial' ##example folder name
  analyze_with='none'
  sort.meth='pvalue'
  # fdr_cut=0.05
  fdr_cut='none'
  load(paste0(analyzed_data.dir,'/DE_result_analyzed_with_',analyze_with,'_',sort.meth,'_sorted_fdr_',fdr_cut,'_cut.rda'))
}
geneNames<-result.list$gene


for_correlation.list=list()
for(sort.meth in c('pvalue','logFC')){
  for_correlation.list[[sort.meth]]=list()
  for(meth in methods_to_plot){
    if(sort.meth=='pvalue'){
      #just avoid log2(pval) to be -inf
      res.temp=(result.list[[meth]]$pval*0.999)+0.001
      res.temp=(-log2(res.temp))*sign(result.list[[meth]]$logFC)
    }else if(sort.meth=='logFC'){
      res.temp=result.list[[meth]]$logFC
    }
    names(res.temp)=rownames(result.list[[meth]])
    res.temp<-res.temp[geneNames]
    for_correlation.list[[sort.meth]][[meth]]=res.temp
  }
}


save(for_correlation.list,file=paste0(analyzed_data.dir,'/corr_heatmap_',analyze_with,'.RData'))
load(paste0(analyzed_data.dir,'/corr_heatmap_',analyze_with,'.RData'))

corr_meth='spearman'
for(sort.meth in c('logFC','pvalue')){
  dir.create(paste0(analyzed_data.dir,'/',sort.meth),recursive = T,showWarnings = F)
  meth.anals<-setdiff(names(result.list),c('genes'))
  if(sort.meth=='pvalue'){
  }else if(sort.meth=='logFC'){
    meth.anals<-meth.anals[str_detect(meth.anals,pattern='Fisher',negate = T)]
  }
  meth.anals<-intersect(methods_to_plot,meth.anals)
  
  cor_matrix<-matrix(NA,ncol=length(meth.anals),nrow=length(meth.anals))
  colnames(cor_matrix)=rownames(cor_matrix)=meth.anals
  for(ma1 in meth.anals){
    for(ma2 in meth.anals){
      cor_matrix[ma1,ma2]=(cor.test( x=for_correlation.list[[sort.meth]][[ma1]],y=for_correlation.list[[sort.meth]][[ma2]],
                                     method = corr_meth,
                                     continuity = FALSE,
                                     conf.level = 0.95)$estimate)
      
    }
  }
  {
    write.table(cor_matrix,file=paste0(analyzed_data.dir,'/',sort.meth,'/corr_heatmap ',corr_meth,' table',sort.meth,' sorted_',analyze_with,'.txt'))
    colors = c(seq(0,1,length=1000))
    my_palette <- colorRampPalette(c("white",brewer.pal(9,"Blues")))(n = 999)
    library(pheatmap)
    callback=function(mat, ...){
      as.hclust(agnes(mat, method = "ward"))
    }
    d<-as.dist(1-cor_matrix)
    for(nk in c(5:15)){
      pdf(paste0(analyzed_data.dir,'/',sort.meth,'/corr_heatmap ',corr_meth,' ',sort.meth,' sorted_',analyze_with,'_agnes_',nk,'_clusters.pdf'),width = 20,height=20)
      hc_cut<-cutree(as.hclust(agnes(d, method = "ward")), k = nk)
      hc_cut_df<-data.frame(Var1=factor(hc_cut))
      pheatmap(mat=cor_matrix,
               color=my_palette,cluster_cols = as.hclust(agnes(d, method = "ward")),
               cluster_rows = as.hclust(agnes(d, method = "ward")),
               cutree_rows = nk,
               cutree_cols = nk,
               annotation_row = hc_cut_df,
               annotation_col = hc_cut_df,fontsize = 20)
      dev.off()
    }
  }
}

sort.meth='pvalue'
cor_matrix<-read.table(file=paste0(analyzed_data.dir,'/',sort.meth,'/corr_heatmap ',corr_meth,' table',sort.meth,' sorted_',analyze_with,'.txt'))
d<-as.dist(1-cor_matrix)
hc3 <- agnes(d, method = "ward")
pdf(paste0(analyzed_data.dir,'/',sort.meth,'/corr_cluster ',corr_meth,' ',sort.meth,' sorted_',analyze_with,'_agnes_',nk,'_clusters.pdf'),width = 16,height=10)
pltree(hc3, cex = 1.8, hang = -2, main = "agnes hcluster(ward)")
rect.hclust(as.hclust(hc3), k = 12, border = 1:13)
dev.off()


##After running correlation analysis for two data, we can compare their clusters.
##Compare two correlation clusters, analysis of scluad epithelial cells data and covid19 monocytes.
comp_dend<-list()
compare_two_dat=c('scLUAD_Epithelial','COVID-19_Monocyte')
for(analyzed_data.dir in compare_two_dat){
  cor_matrix=read.table(file=paste0(analyzed_data.dir,'/',sort.meth,'/corr_heatmap ',corr_meth,' table',sort.meth,' sorted_',analyze_with,'.txt'))
  d<-as.dist(1-cor_matrix)
  hc1 <- agnes(d, method = "ward")
  comp_dend[[analyzed_data.dir]] <- as.dendrogram (hc1)
}
the_cor1=1
the_cor2=cor_bakers_gamma(comp_dend[[1]], comp_dend[[2]])
set.seed(23235)
R <- 1000000
cor_bakers_gamma_results <- numeric(R)
dend_mixed <- comp_dend[[1]]

for(i in 1:R){
  dend_mixed <- sample.dendrogram(dend_mixed, replace = FALSE)
  cor_bakers_gamma_results[i]=cor_bakers_gamma(comp_dend[[1]], dend_mixed)
}
df <- data.frame(
  corr=cor_bakers_gamma_results)
library(ggplot2)
# Basic density
p <- ggplot(df, aes(x=corr)) + 
  geom_density()+
  xlim(-1,1)

p<-p+geom_vline(aes(xintercept=the_cor2,
                    color="compared_correlation"), linetype="dashed", size=1)+
  scale_color_manual(name = "Baker's gamma correlation", values = c(`compared_correlation` = "red"))+
  theme_bw(base_size=15)
p<-p+labs(x="correlation", y="Density")+ggtitle("Baker's gamma correlation")
ggsave(filename = paste0('bakers_gamma_results',paste0(compare_two_dat,collapse = '-'),'.pdf'),width = 13,height = 10)
dev.off()