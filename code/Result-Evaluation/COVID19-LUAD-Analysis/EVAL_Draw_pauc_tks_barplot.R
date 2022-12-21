library(reshape2)
library(magrittr)
library(dplyr)
library(stringr)
library(ggplot2)
base_dir='set base directory/'
setwd(base_dir)
analysis_script_dir=base_dir##Can set different directory. Directory contains scripts for real data analysis(visualization, clustering, correlation, gsea, ...).
source(paste0(analysis_script_dir,'EVAL_Functions for visualization.R'))
width=25
height=16
size_element.text=35

methods_to_plot=c('TCGA_DESeq2','TCGA_edgeR','TCGA_limmavoom','TCGA_limmatrend','GSE43458.limma_Cov','GSE10072.limma_Cov','GSE31210.limma_Cov','Pseudobulk_DESeq2','Pseudobulk_edgeR','Pseudobulk_limmavoom','Pseudobulk_limmatrend','Combat_Wilcox','limma_BEC_Wilcox','MNN_Wilcox','scMerge_Wilcox','Seurat_Wilcox','ZW_BEC_Wilcox','Raw_Wilcox','RISC_Wilcox' , "Scanorama_Wilcox", "scGen_Wilcox", "scVI_Wilcox",'MAST','MAST_Cov','DESeq2','DESeq2_Cov','ZW_DESeq2','ZW_DESeq2_Cov','DESeq2_wFisher','edgeR_DetRate','edgeR_DetRate_Cov','edgeR','edgeR_Cov','ZW_edgeR','ZW_edgeR_Cov','edgeR_wFisher','limmavoom','limmavoom_Cov','limmatrend','limmatrend_Cov','Combat_limmatrend','MNN_limmatrend','scMerge_limmatrend','LogN+limmatrend_wFisher','RISC_limmatrend',"Scanorama_limmatrend", "scGen_limmatrend", "scVI_limmatrend",'DESeq2_FEM','LogN_FEM','DESeq2_REM','LogN_REM','RISC_QP')
fdr_cut=0.05
sort.meth='pvalue'
for(ref in c('Known disease genes','Prognostic genes')){
  pauc<-read.table(paste0(ref,' detected in top ',topn,' genes_',sort.meth, 'sorted_fdr_',fdr_cut,'_cut_pauc.txt'))
  tks<-read.table(paste0('DE_result_',ref,'detection_weighted_ks_',sort.meth,'_sorted_top_',topn_ks*100,'percent_genes_fdr_',fdr_cut,'_cut.txt'))
  
  draw_plot_bar_pauc(pauc,methods_to_plot,ref,sort.meth,top.prop=0.2, width, height, size_element.text)
  draw_plot_bar_pauc(tks,methods_to_plot, ref,sort.meth,top.prop=0.2, width, height, size_element.text)
}
