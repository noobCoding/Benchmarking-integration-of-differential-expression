library(readr)
library(magrittr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(stringr)

base_dir='/hdd2/SC lung/data/'
source('/hdd2/SC lung/script/do.analysis.R')
source('/hdd2/SC lung/script/run_ROC_plot.sources.R')
setwd(base_dir)

FP='none'
meth.comb.list<-get.methcomblist()
args = commandArgs(trailingOnly = T)
if (length(args) < 2)
  stop("Must set ct and subct", call. = F)
ct=args[1]
sub.ct=args[2]
ct%<>%gsub(pattern='_',replacement=' ')
Cell.types<-ct
if(sub.ct=='none'){
  sub.ct=NULL
}else{
  sub.ct%<>%gsub(pattern='_',replacement=' ')
}

for(filter.rate in c(0.05)){
  {
    FDR_thr=0.05 #prognostic gene
    if(FP=='tumor'){
      filt_dir_origin='/fp_filtered/'
    }else if(FP=='normal'){
      filt_dir_origin='/fp_filtered2/'
    }else if(FP=='none'){
      filt_dir_origin='/filtered/'
    }
    
    filt_dir=ifelse(is.null(filter.rate),'',paste0(filt_dir_origin,filter.rate,ifelse(is.null(sub.ct),'',paste0('_',sub.ct))))
    
    fig_dir_origin=paste0(base_dir,'/Tot_result/ps_curve_plots')
    
    dir.create(fig_dir_origin,showWarnings = F)
    fig_dir_origin<-file.path(fig_dir_origin,'changetomethods')
    
    `%ni%`<-Negate(`%in%`)
    up.genes='both'
    
  }
  # for(ref in c('CTD+disgenet','prognostic')){
  for(ref in c('CTD+disgenet')){
    # for(ref in c('CTD','disgenet', 'CTD+disgenet','clin')){
    
    if(ref=='CTD+disgenet'){
      ref_df<-read.delim('/hdd2/SC lung/data/C0152013_disease_gda_summary.tsv')
      ref_df%<>%dplyr::filter(Score_gda>=0.3)
      ref_df%<>%dplyr::arrange(desc(Score_gda))
      ref_genes.1<-ref_df$Gene
      ref_weight.1<-ref_df$Score_gda
      ref_df<-read.delim('/hdd2/SC lung/data/CTD_D000077192_genes_20210713220641.tsv')
      ref_df%<>%dplyr::filter(Direct.Evidence!="")
      ref_genes.2<-ref_df$Gene.Symbol
      both.in<-intersect(ref_genes.1,ref_genes.2)
      # union(ref_genes.1,ref_genes.2)%>%length()
      # median(ref_weight.1[ref_genes.1%in%both.in])
      # mean(ref_weight.1[ref_genes.1%in%both.in])
      ref_genes<-c(ref_genes.1,setdiff(ref_genes.2,ref_genes.1))
      
      ref_weight<-c(ref_weight.1, rep(median(ref_weight.1[ref_genes.1%in%both.in]), setdiff(ref_genes.2,ref_genes.1)%>%length()))
    }else if(ref=='prognostic_new'){
      FDR_thr=0.05
      meta.res<-read.table('meta.clinic.integrated.filtered_new_max.txt')
      ref_genes<-meta.res%>%dplyr::filter(FDR<FDR_thr)%>%rownames()
      ref_weight=rep(1,ref_genes%>%length())
      # ref_weight=meta.res[ref_genes,'logHR']
    }else if(ref=='prognostic_new.pos'){
      FDR_thr=0.05
      meta.res<-read.table('meta.clinic.integrated.filtered_new_max.txt')
      ref_genes<-meta.res%>%dplyr::filter(FDR<FDR_thr)%>%dplyr::filter(logHR>0)%>%rownames()
      ref_weight=rep(1,ref_genes%>%length())
    }else if(ref=='prognostic_new.mean'){
      FDR_thr=0.05
      meta.res<-read.table('meta.clinic.integrated.filtered_new_mean.txt')
      ref_genes<-meta.res%>%dplyr::filter(FDR<FDR_thr)%>%rownames()
      ref_weight=rep(1,ref_genes%>%length())
    }else if(ref=='prognostic_new.mean.pos'){
      FDR_thr=0.05
      meta.res<-read.table('meta.clinic.integrated.filtered_new_mean.txt')
      ref_genes<-meta.res%>%dplyr::filter(FDR<FDR_thr)%>%dplyr::filter(logHR>0)%>%rownames()
      ref_weight=rep(1,ref_genes%>%length())
    }else if(ref=='prognostic_old.mean'){
      FDR_thr=0.05
      meta.res<-read.table('meta.clinic.integrated.filtered_old.mean.txt')
      ref_genes<-meta.res%>%dplyr::filter(FDR<FDR_thr)%>%rownames()
      ref_weight=rep(1,ref_genes%>%length())
    }else if(ref=='prognostic_old.max'){
      FDR_thr=0.05
      meta.res<-read.table('meta.clinic.integrated.filtered_old.max.txt')
      ref_genes<-meta.res%>%dplyr::filter(FDR<FDR_thr)%>%rownames()
      ref_weight=rep(1,ref_genes%>%length())
    }else if(ref=='prognostic_old.mean.pos'){
      FDR_thr=0.05
      meta.res<-read.table('meta.clinic.integrated.filtered_old.mean.txt')
      ref_genes<-meta.res%>%dplyr::filter(FDR<FDR_thr)%>%dplyr::filter(logHR>0)%>%rownames()
      ref_weight=rep(1,ref_genes%>%length())
    }else if(ref=='prognostic_old.max.pos'){
      FDR_thr=0.05
      meta.res<-read.table('meta.clinic.integrated.filtered_old.max.txt')
      ref_genes<-meta.res%>%dplyr::filter(FDR<FDR_thr)%>%dplyr::filter(logHR>0)%>%rownames()
      ref_weight=rep(1,ref_genes%>%length())
    }else if(ref=='prognostic'){
      FDR_thr=0.05
      meta.res<-read.table('meta.clinic.integrated.filtered.txt')
      ref_genes<-meta.res%>%dplyr::filter(FDR<FDR_thr)%>%rownames()
      ref_weight=rep(1,ref_genes%>%length())
    }
    rm(ref_df)
    # for(adj.p_thr in c(0.01,0.05,'none')){
    for(adj.p_thr in c('none',0.01,0.05)){
      # for(adj.p_thr in c(0.01)){
      # for(sort.meth in c('logFC','adj.pvalue','pvalue')){
      # for(sort.meth in c('logFC','pvalue')){
      for(sort.meth in c('pvalue')){
        fig_dir<-paste0(fig_dir_origin,'/',sort.meth,'/')
        # dir.create(fig_dir, showWarnings = F)
        
        if(str_detect(string=ref,pattern='prognostic.*cutoff')){
          ref.use='prognostic'
        }else{
          ref.use=ref
        }
        
        
        # saved_dir=paste0(base_dir,"/Tot_result/",ifelse(any(str_detect(meth.comb.tot,pattern = 'TCGA')),'sc+TCGA','sc'),"/",ct,ifelse(is.null(sub.ct),'',paste0('_',sub.ct)),'/',sort.meth,'/')
        
        # compmethods.all<-list()
        # for(compmethods in c('DESeq2+TCGA','DESeq2+pseudobulk', 'edgeR+TCGA','edgeR+pseudobulk','Total+TCGA','Total+pseudobulk',"added_comp")){
        # for(compmethods in c("Total+TCGA+zinbcov")){
        
        # for(compmethods in c("Total+Hai_comp","Total+TCGA+Hai_comp","Total+Hai_comp_add","Total+Hai_comp_add2","Total+Hai_comp_add3")){
        for(compmethods in c('Main_figure+sc','Main_figure','Main_figure+add','Main_figure+add2','Main_figure+add3')){
          compmethods.all<-list()
          #지금 compmethods별로 input folder가 달라져서 여기서 list 실행해야함.
          # ,"Total+TCGA+Hai_comp+add","Total+TCGA+Hai_comp+add2","Total+TCGA+Hai_comp+add3"
          compmethods.all[[compmethods]]=meth.comb.list[[compmethods]]
          # fig_dir<-paste0(fig_dir_origin,'/',sort.meth,'/')
          dir.create(gsub(pattern='changetomethods',replacement =compmethods,x=fig_dir_origin), showWarnings = F)
          dir.create(gsub(pattern='changetomethods',replacement =compmethods,x=fig_dir), showWarnings = F)
          
          # saved_dir=paste0(base_dir,"/Tot_result/",c('sc','sc+TCGA','sc+add','sc+add2','sc+add3')[which(compmethods==c("Total+Hai_comp","Total+TCGA+Hai_comp","Total+Hai_comp_add","Total+Hai_comp_add2","Total+Hai_comp_add3"))],"/",ct,ifelse(is.null(sub.ct),'',paste0('_',sub.ct)),'/',sort.meth,'/')
          saved_dir=paste0(base_dir,"/Tot_result/",c('sc','sc+TCGA','sc+add','sc+add2','sc+add3')[which(compmethods==c('Main_figure+sc','Main_figure','Main_figure+add','Main_figure+add2','Main_figure+add3'))],"/",ct,ifelse(is.null(sub.ct),'',paste0('_',sub.ct)),'/',sort.meth,'/')
          draw_plot(saved_dir=saved_dir,base_dir=base_dir,meth.comb.tot=meth.comb.tot,meth.comb.use = compmethods.all,
                    ct=ct,filter.rate=filter.rate,adj.p_thr=adj.p_thr,ref=ref,ref_genes=ref_genes,ref_weight=ref_weight,
                    sub.ct=sub.ct,filt_dir=filt_dir,fig_dir=fig_dir, sort.meth=sort.meth,up.genes=up.genes,plot_topn_genes=c(0.1,0.2,'all'))
        }
        # saved_dir=paste0(base_dir,"/Tot_result/",ifelse(any(str_detect(unique(unlist(compmethods.all)),pattern = 'TCGA')),'sc+TCGA','sc'),"/",ct,ifelse(is.null(sub.ct),'',paste0('_',sub.ct)),'/',sort.meth,'/')
        # draw_plot(saved_dir=saved_dir,base_dir=base_dir,meth.comb.tot=meth.comb.tot,meth.comb.use = compmethods.all,
        #           ct=ct,filter.rate=filter.rate,adj.p_thr=adj.p_thr,ref=ref,ref_genes=ref_genes,ref_weight=ref_weight,
        #           sub.ct=sub.ct,filt_dir=filt_dir,fig_dir=fig_dir, sort.meth=sort.meth,up.genes=up.genes,plot_topn_genes=c(0.1,0.2,'all'))
        
        
      }
      
    }
    
  }
}