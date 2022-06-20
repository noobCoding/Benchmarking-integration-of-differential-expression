library(readr)
library(magrittr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(stringr)

base_dir='~/'
source(paste0(base_dir,'script_bk/do.analysis.R'))
source(paste0(base_dir,'script_bk/visualize.sources.R'))
setwd(base_dir)
meth.comb.list<-get.methcomblist()
for(ct in c('Epithelial cells','Myeloid cells', 'T-NK cells')){
  for(filter.rate in c(0.05)){
    {
      sub.ct=NULL
      filt_dir=ifelse(is.null(filter.rate),'',paste0('/filtered/',filter.rate,ifelse(is.null(sub.ct),'',paste0('_',sub.ct))))
      `%ni%`<-Negate(`%in%`)
      up.genes='both'
      adj.p_thr='none'
    }
    for(sort.meth in c('logFC','pvalue')){
      saved_dir.origin<-paste0(base_dir,"data/analyzed/scLUAD/",ct,ifelse(is.null(sub.ct),'',paste0('_',sub.ct)),'/',sort.meth,'/')
      saved_dir<-saved_dir.origin
      write_ks_df(meth.comb.tot=c(meth.comb.list[['Total']]),
                  ct=ct,base_dir=base_dir,filter.rate=filter.rate,adj.p_thr=adj.p_thr,ref=ref.use,ref_df_refined=ref_df_refined,
                  sub.ct=sub.ct,filt_dir=filt_dir,fig_dir=saved_dir, sort.meth=sort.meth,up.genes=up.genes)
      write.ks.adjp_thr(saved_dir=saved_dir,ct=ct,base_dir=base_dir,filter.rate=filter.rate,filt_dir=filt_dir,adj.p_thr=0.05,sub.ct=sub.ct,sort.meth=sort.meth,up.genes=up.genes)

      saved_dir<-gsub(pattern='[/]scLUAD[/]',replacement='/scLUAD-TCGA/',saved_dir.origin)
      write_ks_df(meth.comb.tot=c(meth.comb.list[['Total+TCGA']]),
                  ct=ct,base_dir=base_dir,filter.rate=filter.rate,adj.p_thr=adj.p_thr,ref=ref.use,ref_df_refined=ref_df_refined,
                  sub.ct=sub.ct,filt_dir=filt_dir,fig_dir=saved_dir, sort.meth=sort.meth,up.genes=up.genes)
      write.ks.adjp_thr(saved_dir=saved_dir,ct=ct,base_dir=base_dir,filter.rate=filter.rate,filt_dir=filt_dir,adj.p_thr=0.05,sub.ct=sub.ct,sort.meth=sort.meth,up.genes=up.genes)

      saved_dir<-gsub(pattern='[/]scLUAD[/]',replacement='/scLUAD-GSE43458/',saved_dir.origin)
      dir.create(saved_dir,showWarnings = T,recursive = T)
      write_ks_df(meth.comb.tot=c(meth.comb.list[['Total+GSE43458']]),
                  ct=ct,base_dir=base_dir,filter.rate=filter.rate,adj.p_thr=adj.p_thr,ref=ref.use,ref_df_refined=ref_df_refined,
                  sub.ct=sub.ct,filt_dir=filt_dir,fig_dir=saved_dir, sort.meth=sort.meth,up.genes=up.genes)
      write.ks.adjp_thr(saved_dir=saved_dir,ct=ct,base_dir=base_dir,filter.rate=filter.rate,filt_dir=filt_dir,adj.p_thr=0.05,sub.ct=sub.ct,sort.meth=sort.meth,up.genes=up.genes)

      
      saved_dir<-gsub(pattern='[/]scLUAD[/]',replacement='/scLUAD-GSE31210/',saved_dir.origin)
      dir.create(saved_dir,showWarnings = T,recursive = T)
      write_ks_df(meth.comb.tot=c(meth.comb.list[['Total+GSE31210']]),
                  ct=ct,base_dir=base_dir,filter.rate=filter.rate,adj.p_thr=adj.p_thr,ref=ref.use,ref_df_refined=ref_df_refined,
                  sub.ct=sub.ct,filt_dir=filt_dir,fig_dir=saved_dir, sort.meth=sort.meth,up.genes=up.genes)
      write.ks.adjp_thr(saved_dir=saved_dir,ct=ct,base_dir=base_dir,filter.rate=filter.rate,filt_dir=filt_dir,adj.p_thr=0.05,sub.ct=sub.ct,sort.meth=sort.meth,up.genes=up.genes)
    }
  }
}