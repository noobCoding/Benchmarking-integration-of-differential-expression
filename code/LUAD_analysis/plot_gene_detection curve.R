library(readr)
library(magrittr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(stringr)

base_dir='~/'
source('/hdd2/SC lung/script_bk/do.analysis.R')
source('/hdd2/SC lung/script_bk/visualize.sources.R')
setwd(base_dir)
FP='none'
meth.comb.list<-get.methcomblist()

for(ct in c('Epithelial cells','Myeloid cells','T-NK cells')){
  for(filter.rate in c(0.05)){
    {
      FDR_thr=0.05 #prognostic gene
      sub.ct=NULL
      if(FP=='tumor'){
        filt_dir_origin='/fp_filtered/'
      }else if(FP=='normal'){
        filt_dir_origin='/fp_filtered2/'
      }else if(FP=='none'){
        filt_dir_origin='/filtered/'
      }
      filt_dir=ifelse(is.null(filter.rate),'',paste0(filt_dir_origin,filter.rate,ifelse(is.null(sub.ct),'',paste0('_',sub.ct))))
      
      fig_dir_origin=paste0(base_dir,'data/analyzed/gene_detection_curve')
      
      dir.create(fig_dir_origin,showWarnings = F)
      fig_dir_origin<-file.path(fig_dir_origin,'changetomethods')
      
      `%ni%`<-Negate(`%in%`)
      up.genes='both'
      
    }
    
    for(ref in c('Known disease genes','Prognostic genes')){
      if(ref=='Known disease genes'){
        ref_df<-read.delim(paste0(base_dir,'data/original/C0152013_disease_gda_summary.tsv'))
        ref_df%<>%dplyr::filter(Score_gda>=0.3)
        ref_df%<>%dplyr::arrange(desc(Score_gda))
        ref_genes.1<-ref_df$Gene
        ref_weight.1<-ref_df$Score_gda
        ref_df<-read.delim(paste0(base_dir,'data/original/CTD_D000077192_genes_20210713220641.tsv'))
        ref_df%<>%dplyr::filter(Direct.Evidence=="marker/mechanism")
        ref_genes.2<-ref_df$Gene.Symbol
        both.in<-intersect(ref_genes.1,ref_genes.2)
        # median(ref_weight.1[ref_genes.1%in%both.in])
        # mean(ref_weight.1[ref_genes.1%in%both.in])
        ref_genes<-c(ref_genes.1,setdiff(ref_genes.2,ref_genes.1))
        ref_weight<-c(ref_weight.1, rep(median(ref_weight.1[ref_genes.1%in%both.in]), setdiff(ref_genes.2,ref_genes.1)%>%length()))
      }else if(ref=='Prognostic genes'){
        FDR_thr=0.05
        meta.res<-read.table(paste0(base_dir,'data/analyzed/for surv/meta.clinic.integrated.filtered.txt'))
        ref_genes<-meta.res%>%dplyr::filter(FDR<FDR_thr)%>%rownames()
        ref_weight=rep(1,ref_genes%>%length())
      }
      rm(ref_df)
      
      for(adj.p_thr in c('none',0.05)){
        
        for(sort.meth in c('pvalue')){
          fig_dir<-paste0(fig_dir_origin,'/',sort.meth,'/')
          
          for(compmethods in c('Total',"Total+TCGA","Total+GSE43458","Total+GSE31210")){
            compmethods.all<-list()
            compmethods.all[[compmethods]]=meth.comb.list[[compmethods]]
            dir.create(gsub(pattern='changetomethods',replacement =compmethods,x=fig_dir_origin), showWarnings = F)
            dir.create(gsub(pattern='changetomethods',replacement =compmethods,x=fig_dir), showWarnings = F)
            
            saved_dir=paste0(base_dir,"data/analyzed/",c('scLUAD','scLUAD-TCGA','scLUAD-GSE43458','scLUAD-GSE31210')[which(compmethods==c('Total',"Total+TCGA","Total+GSE43458","Total+GSE31210"))],"/",ct,ifelse(is.null(sub.ct),'',paste0('_',sub.ct)),'/',sort.meth,'/')
            draw_plot(saved_dir=saved_dir,base_dir=base_dir,meth.comb.tot=meth.comb.tot,meth.comb.use = compmethods.all,
                      ct=ct,filter.rate=filter.rate,adj.p_thr=adj.p_thr,ref=ref,ref_genes=ref_genes,ref_weight=ref_weight,
                      sub.ct=sub.ct,filt_dir=filt_dir,fig_dir=fig_dir, sort.meth=sort.meth,up.genes=up.genes,plot_topn_genes=c(0.2))
          }
          
          
          
        }
        
      }
      
    }
  }
}


