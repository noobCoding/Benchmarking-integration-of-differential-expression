base_dir='set base directory/'#Directory contains DE_result.rda file.
setwd(base_dir)
data_analyze_with_dir=base_dir##Can set different directory. Directory contains TCGA/GSE43458/GSE31210 analysis results.
for(analyze_with in c('none','TCGA','GSE43458',"GSE31210")){
  genes_analyze_with=data_analyze_with=NULL
  if(analyze_with=='TCGA'){
    data_analyze_with=list()
    for(meth in c('DESeq2','edgeR','limmatrend','limmavoom')){
      result.table.with<-read.table(paste0(data_analyze_with_dir,'TCGA-LUAD.rnaseq.processed.htseq_count.',meth,'_Cov.txt'))
      result.table=data.frame(pval=result.table.with$pvalue,qval=result.table.with$adjpvalue,logFC=result.table.with$logFC,row.names = rownames(result.table.with))
      data_analyze_with[[paste0(analyze_with,'_',meth)]]=result.table
      if(!is.null(genes_analyze_with)){
        genes_analyze_with<-intersect(genes_analyze_with,rownames(result.table))
      }else{
        genes_analyze_with<-rownames(result.table)
      }
      rm(result.table.with,result.table)
    }
  }else if(analyze_with%in%c('GSE43458',"GSE31210")){
    data_analyze_with=list()
    result.table.with<-read.table(paste0(data_analyze_with_dir,'LUAD-',analyze_with,'.microarray.limma.cov.txt'))
    result.table=data.frame(pval=result.table.with$pvalue,qval=result.table.with$adjpvalue,logFC=result.table.with$logFC,row.names = rownames(result.table.with))
    data_analyze_with[[paste0(analyze_with,'.limmatrend_Cov')]]=result.table
    genes_analyze_with<-rownames(result.table)
    rm(result.table.with,result.table)
  }
  for(sort.meth in c('abs(logFC)','pvalue')){
    for(fdr_cut in c(0.01,0.05,'none')){
      load('DE_result.rda')
      geneNames<-result.list$gene
      result.list[['gene']]=NULL
      if(!is.null(data_analyze_with)){
        geneNames<-intersect(geneNames,genes_analyze_with)
        result.list<-append(data_analyze_with,result.list)
        result.list<-lapply(result.list,FUN=function(x){x[geneNames,]})
      }
      if(fdr_cut!='none'){
        result.list<-lapply(result.list,FUN=function(x){x%>%dplyr::filter(qval<0.05)})
      }
      if(sort.meth=='abs(logFC)'){
        result.list<-lapply(result.list,FUN=function(x){
          sort=abs(x$logFC)
          sort2=rank(-(sort),ties.method='max')
          x$ranks=sort2
          x%<>%dplyr::arrange(ranks)
        })
      }else if(sort.meth=='pvalue'){
        result.list<-lapply(result.list,FUN=function(x){
          sort=x$pval
          abs.logfc<-abs(x$logFC)
          sort<-rank(sort,ties.method = 'max')
          sort2<-sort+(1-(abs.logfc/max(abs.logfc)))*0.5
          sort2<-rank(sort2,ties.method = 'max')
          x$ranks=sort2
          x%<>%dplyr::arrange(ranks)
        })
      }
      result.list<-append(list(gene=geneNames),result.list)
      save(result.list,file=paste0('DE_result_analyzed_with_',analyze_with,'_',sort.meth,'_sorted_fdr_',fdr_cut,'_cut.rda'))
    }
  }
}



