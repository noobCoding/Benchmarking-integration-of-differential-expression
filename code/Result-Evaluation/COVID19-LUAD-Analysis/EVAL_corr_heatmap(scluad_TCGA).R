library(magrittr)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(cluster)
library(factoextra)
library(dendextend)
base_dir='set base directory/'##Directory contains scLUAD&TCGA DE analysis results.
setwd(base_dir)

methods_to_plot<-c('TCGA_DESeq2','TCGA_edgeR','TCGA_limmavoom','TCGA_limmatrend','Pseudobulk_DESeq2','Pseudobulk_edgeR','Pseudobulk_limmavoom','Pseudobulk_limmatrend','Combat_Wilcox','limma_BEC_Wilcox','MNN_Wilcox','scMerge_Wilcox','Seurat_Wilcox','ZW_BEC_Wilcox','Raw_Wilcox','RISC_Wilcox' , "Scanorama_Wilcox", "scGen_Wilcox", "scVI_Wilcox",'MAST','MAST_Cov','DESeq2','DESeq2_Cov','ZW_DESeq2','ZW_DESeq2_Cov','DESeq2_wFisher','edgeR_DetRate','edgeR_DetRate_Cov','edgeR','edgeR_Cov','ZW_edgeR','ZW_edgeR_Cov','edgeR_wFisher','limmavoom','limmavoom_Cov','limmatrend','limmatrend_Cov','Combat_limmatrend','MNN_limmatrend','scMerge_limmatrend','LogN+limmatrend_wFisher','RISC_limmatrend',"Scanorama_limmatrend", "scGen_limmatrend", "scVI_limmatrend",'DESeq2_FEM','LogN_FEM','DESeq2_REM','LogN_REM','RISC_QP')
{
  analyze_with='TCGA'
  sort.meth='pvalue'
  # fdr_cut=0.05
  fdr_cut='none'
  load(paste0('DE_result_analyzed_with_',analyze_with,'_',sort.meth,'_sorted_fdr_',fdr_cut,'_cut.rda'))
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


save(for_correlation.list,file='corr_heatmap_sclung_TCGA.RData')
load('corr_heatmap_sclung_TCGA.RData')

corr_meth='spearman'
for(sort.meth in c('logFC','pvalue')){
  dir.create(paste0(sort.meth),recursive = T,showWarnings = F)
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
    write.table(cor_matrix,file=paste0(sort.meth,'/corr_heatmap ',corr_meth,' table',sort.meth,' sorted_sclung_TCGA.txt'))
    colors = c(seq(0,1,length=1000))
    my_palette <- colorRampPalette(c("white",brewer.pal(9,"Blues")))(n = 999)
    library(pheatmap)
    callback=function(mat, ...){
      as.hclust(agnes(mat, method = "ward"))
    }
    d<-as.dist(1-cor_matrix)
    for(nk in c(5:15)){
      pdf(paste0(sort.meth,'/corr_heatmap ',corr_meth,' ',sort.meth,' sorted_sclung_TCGA_agnes_',nk,'_clusters.pdf'),width = 20,height=20)
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
cor_matrix<-read.table(file=paste0(sort.meth,'/corr_heatmap ',corr_meth,' table',sort.meth,' sorted_sclung_TCGA.txt'))
d<-as.dist(1-cor_matrix)
hc3 <- agnes(d, method = "ward")
pdf(paste0(sort.meth,'/corr_cluster ',corr_meth,' ',sort.meth,' sorted_sclung_TCGA_agnes_',nk,'_clusters.pdf'),width = 16,height=10)
pltree(hc3, cex = 1.8, hang = -2, main = "sclung TCGA agnes hcluster(ward)")
rect.hclust(as.hclust(hc3), k = 12, border = 1:13)
dev.off()
