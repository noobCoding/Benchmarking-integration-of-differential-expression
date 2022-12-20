# setwd('~/')
result.list=list()
geneNames<-c()
for(i in c("combat+limmatrend.rda", "combat+wilcox.rda", "deseq2_Cov.rda", "deseq2_sep_sep+FEM.rda", "deseq2_sep_sep+REM.rda", "deseq2_sep_sep+wfisher.rda",
           "deseq2.rda", "edgeR_Cov.rda", "edgeR_Detrate_Cov.rda", "edgeR_Detrate.rda", "edgeR_sep_sep+wfisher.rda", "edgeR.rda",
           "limma_bec+wilcox.rda", "limmatrend_Cov.rda", "limmatrend.rda", "limmavoom_Cov.rda", "limmavoom.rda", "LogNormalize_sep_sep+es_sep_sep+FEM.rda",
           "LogNormalize_sep_sep+es_sep_sep+REM.rda", "LogNormalize_sep_sep+limmatrend_sep_sep+wfisher.rda","MAST_Cov.rda", "MAST.rda", "mnn+wilcox.rda",'mnn+limmatrend.rda',
           "pseudobulk+deseq2.rda", "pseudobulk+edgeR.rda", "pseudobulk+limmatrend.rda", "pseudobulk+limmavoom.rda", "risc+limmatrend.rda", "risc+QP.rda", "risc+wilcox.rda",
           "scmerge+limmatrend.rda", "scmerge+wilcox.rda", "Seurat+wilcox.rda", "wilcox.rda", "zinbwave+deseq2_Cov.rda", "zinbwave+deseq2.rda", "zinbwave+edgeR_Cov.rda",
           "zinbwave+edgeR.rda", "zinbwave+wilcox.rda" )){
  # , "scanorama+wilcox", "scanorama+limmatrend","scgen+wilcox","scgen+limmatrend","scvi+wilcox","scvi+limmatrend",'edgeR_sep.rda'
  load(i)
  i.use<-i%>%gsub(pattern='[.]rda$',replacement='')
  methods<-str_split(i.use,pattern='[+]')[[1]]
  
  if(methods[length(methods)]=='wilcox'){
    result.table=data.frame(pval=res$p_val,qval=res$p_val_adj,logFC=res$avg_log2FC,row.names = rownames(res))
  }else if(methods[length(methods)]%in%c('MAST','MAST_Cov')){
    result.table=data.frame(pval=res$p_val,qval=res$p_val_adj,logFC=res$avg_log2FC,row.names = rownames(res))
  }else if(methods[length(methods)]%in%c('limmatrend','limmatrend_Cov','limmavoom','limmavoom_Cov')){
    result.table=data.frame(pval=res$pvalue,qval=res$adjpvalue,logFC=res$logFC,row.names = rownames(res))
  }else if(methods[length(methods)]%in%c('deseq2','deseq2_Cov')){
    result.table=data.frame(pval=res$pvalue,qval=res$adjpvalue,logFC=res$logFC,row.names = rownames(res))
  }else if(methods[length(methods)]%in%c('edgeR','edgeR_Cov','edgeR_Detrate','edgeR_Detrate_Cov')){
    result.table=data.frame(pval=res$pvalue,qval=res$adjpvalue,logFC=res$logFC,row.names = rownames(res))
  }else if(methods[length(methods)]%in%c('REM','FEM')){
    result.table=data.frame(pval=res$pval,qval=res$FDR,logFC=res$mu.hat,row.names = rownames(res))
  }else if(methods[length(methods)]=='wfisher'){
    res$direction[res$direction=='+']=1
    res$direction[res$direction=='-']=-1
    res$direction%<>%as.numeric()
    result.table=data.frame(pval=res$pval,qval=res$FDR,logFC=res$direction,row.names = rownames(res))
  }else if(methods[length(methods)]=='QP'){
    result.table=data.frame(pval=res$Pvalue,qval=res$Padj,logFC=res$logFC,row.names = rownames(res))
  }else if(str_detect(methods[length(methods)],pattern='[_]sep$')){
    for(ind.batch in colnames(res$abs)){
      result.table=data.frame(pval=res$abs[,ind.batch],qval=res$abs[,ind.batch],logFC=res$logFC[,ind.batch],row.names=rownames(res))
      result.table$qval=p.adjust(p=result.table$qval,method = 'BH')
      geneNames%<>%c(rownames(result.table))
      geneNames%<>%unique()
      result.list[[meth.name.convert(paste0(i.use,'_ind.analysis_',ind.batch))]]=result.table
    }
  }
  
  if(!str_detect(methods[length(methods)],pattern='[_]sep$')){
    geneNames%<>%c(rownames(result.table))
    geneNames%<>%unique()
    result.list[[meth.name.convert(i.use)]]=result.table
  }
}
geneNames%<>%sort()
for(i in names(result.list)){
  result.table<-result.list[[i]]
  result.table<-result.table[geneNames,]
  rownames(result.table)=geneNames
  result.table[is.na(result.table$pval),]=c(1,1,0)
  result.list[[i]]<-result.table
}
# result.list <- lapply(result.list, function(x) c('x', x))
# result.list[['gene']]=geneNames
result.list<-append(list(gene=geneNames),result.list)
save(result.list,file='DE_result.rda')




meth.name.convert<-function(x){
  inputnames<-c("combat+limmatrend", "combat+wilcox", "deseq2_Cov", "deseq2_sep_sep+FEM", "deseq2_sep_sep+REM", "deseq2_sep_sep+wfisher",
                "deseq2", "edgeR_Cov", "edgeR_Detrate_Cov", "edgeR_Detrate", "edgeR_sep_sep+wfisher", "edgeR",
                "limma_bec+wilcox", "limmatrend_Cov", "limmatrend", "limmavoom_Cov", "limmavoom", "LogNormalize_sep_sep+es_sep_sep+FEM",
                "LogNormalize_sep_sep+es_sep_sep+REM", "LogNormalize_sep_sep+limmatrend_sep_sep+wfisher","MAST_Cov", "MAST", "mnn+wilcox",'mnn+limmatrend',
                "pseudobulk+deseq2", "pseudobulk+edgeR", "pseudobulk+limmatrend", "pseudobulk+limmavoom", "risc+limmatrend", "risc+QP", "risc+wilcox",
                "scmerge+limmatrend", "scmerge+wilcox", "Seurat+wilcox", "wilcox", "zinbwave+deseq2_Cov", "zinbwave+deseq2", "zinbwave+edgeR_Cov",
                "zinbwave+edgeR", "zinbwave+wilcox", "scanorama+wilcox", "scanorama+limmatrend","scgen+wilcox","scgen+limmatrend","scvi+wilcox","scvi+limmatrend")
  outputnames<-c('Combat_limmatrend','Combat_Wilcox','DESeq2_Cov','DESeq2_FEM','DESeq2_REM','DESeq2_wFisher',
                 'DESeq2','edgeR_Cov','edgeR_DetRate_Cov','edgeR_DetRate','edgeR_wFisher','edgeR',
                 'limma_BEC_Wilcox','limmatrend_Cov','limmatrend','limmavoom_Cov','limmavoom','LogN_FEM',
                 'LogN_REM','LogN+limmatrend_wFisher','MAST_Cov','MAST','MNN_Wilcox','MNN_limmatrend',
                 'Pseudobulk_DESeq2','Pseudobulk_edgeR','Pseudobulk_limmatrend','Pseudobulk_limmavoom','RISC_limmatrend','RISC_QP','RISC_Wilcox',
                 'scMerge_limmatrend','scMerge_Wilcox','Seurat_Wilcox','Raw_Wilcox','ZW_DESeq2_Cov','ZW_DESeq2','ZW_edgeR_Cov',
                 'ZW_edgeR','ZW_BEC_Wilcox', "Scanorama_Wilcox","Scanorama_limmatrend", "scGen_Wilcox", "scGen_limmatrend", "scVI_Wilcox", "scVI_limmatrend")
  if(any(inputnames==x)){
    return(outputnames[which(inputnames==x)])
  }else{
    return(x)
  }
  
}
