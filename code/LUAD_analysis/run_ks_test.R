library(dplyr)
library(tidyverse)
library(magrittr)
source('/hdd2/SC lung/script/run_ROC_plot.sources.R')
source('/hdd2/SC lung/script/ks_wbk.R')
base_dir='/hdd2/SC lung/data/'
setwd(base_dir)
meth.comb.list<-get.methcomblist()
for(ct in c('Myeloid cells','T-NK cells','Epithelial cells')){
  for(filter.rate in c(0.05)){
    for(adj.p_thr in c(0.01,0.05,'none')){
      {
        sub.ct=NULL
        `%ni%`<-Negate(`%in%`)
        up.genes='both'
        # adj.p_thr='none'
        ###ks params###
        topn_ks=0.2 #top 20% truncated ks test
      }
      for(ref in c('Known disease genes','Prognostic genes')){
        if(ref=='Known disease genes'){
          ref_df<-read.delim(paste0(base_dir,'/original/C0152013_disease_gda_summary.tsv'))
          ref_df%<>%dplyr::filter(Score_gda>=0.3)
          ref_df%<>%dplyr::arrange(desc(Score_gda))
          ref_genes.1<-ref_df$Gene
          ref_weight.1<-ref_df$Score_gda
          ref_df<-read.delim(paste0(base_dir,'/original/CTD_D000077192_genes_20210713220641.tsv'))
          ref_df%<>%dplyr::filter(Direct.Evidence=="marker/mechanism")
          ref_genes.2<-ref_df$Gene.Symbol
          both.in<-intersect(ref_genes.1,ref_genes.2)
          # median(ref_weight.1[ref_genes.1%in%both.in])
          # mean(ref_weight.1[ref_genes.1%in%both.in])
          ref_genes<-c(ref_genes.1,setdiff(ref_genes.2,ref_genes.1))
          ref_weight<-c(ref_weight.1, rep(median(ref_weight.1[ref_genes.1%in%both.in]), setdiff(ref_genes.2,ref_genes.1)%>%length()))
        }else if(ref=='Prognostic genes'){
          FDR_thr=0.05
          meta.res<-read.table(paste0(base_dir,'/for surv/meta.clinic.integrated.filtered.txt'))
          ref_genes<-meta.res%>%dplyr::filter(FDR<FDR_thr)%>%rownames()
          ref_weight=rep(1,ref_genes%>%length())
        }
        for(sort.meth in c('logFC','pvalue')){
          for(comp.with in c('Total','Total+TCGA','Total+GSE43458','Total+GSE31210')){
            if(comp.with=='Total'){
              saved_dir<-paste0(base_dir,"/analyzed/scLUAD/",ct,ifelse(is.null(sub.ct),'',paste0('_',sub.ct)),'/',sort.meth,'/')
              meth.comb.tot=meth.comb.list[['Total']]
            }else if(comp.with=='Total+TCGA'){
              saved_dir<-paste0(base_dir,"/analyzed/scLUAD-TCGA/",ct,ifelse(is.null(sub.ct),'',paste0('_',sub.ct)),'/',sort.meth,'/')
              meth.comb.tot=meth.comb.list[['Total+TCGA']]
            }else if(comp.with=='Total+GSE43458'){
              saved_dir<-paste0(base_dir,"/analyzed/scLUAD-GSE43458/",ct,ifelse(is.null(sub.ct),'',paste0('_',sub.ct)),'/',sort.meth,'/')
              meth.comb.tot=meth.comb.list[['Total+GSE43458']]
            }else if(comp.with=='Total+GSE31210'){
              saved_dir<-paste0(base_dir,"/analyzed/scLUAD-GSE31210/",ct,ifelse(is.null(sub.ct),'',paste0('_',sub.ct)),'/',sort.meth,'/')
              meth.comb.tot=meth.comb.list[['Total+GSE31210']]
            }

            load(paste0(saved_dir,up.genes,' | ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),'gene threshold : expressed in  ',filter.rate*100,'percent cells | adj.p_threshold : ',adj.p_thr, '|',sort.meth, ' sorted.ks.RData'))
            res_df<-sapply(setdiff(names(res_list),'genes'),FUN = function(x){match(res_list[['genes']],res_list[[x]][['genes']])})
            rownames(res_df)=res_list[['genes']]
            res_df%<>%as.data.frame()
            res_df_ks<-res_df[which(rownames(res_df)%in%ref_genes),]
            final_df_ks<-matrix(NA,nrow=ncol(res_df_ks),ncol=1)
            rownames(final_df_ks)=colnames(res_df_ks)
            for(i in 1:ncol(res_df_ks)){
              meth.comb=colnames(res_df_ks)[i]
              print(meth.comb)
              
              if(!is.null(topn_ks)){
                if(topn_ks>1){
                  topn_ks.use=as.numeric(topn_ks)
                  topn_ks.write=paste0('top',as.numeric(topn_ks))
                }else{
                  topn_ks.use=round(as.numeric(topn_ks)*nrow(res_df))
                  topn_ks.write=paste0('top',as.numeric(topn_ks)*100,'percent')
                }
                if(sum(!is.na(res_df[,meth.comb]))<topn_ks.use){
                  topn_ks.use=sum(!is.na(res_df[,meth.comb]))
                }
                
                topn_ord_ks<-res_df_ks[,meth.comb]
                
                
                topn_ks_last<-topn_ks.use+(1:(length(topn_ord_ks)-sum(topn_ord_ks<=topn_ks.use, na.rm = T)))*(nrow(res_df)-topn_ks.use)/(length(topn_ord_ks)-sum(topn_ord_ks<=topn_ks.use, na.rm = T))
                
                w.x=c(ref_weight[rownames(res_df_ks)][which(topn_ord_ks<=topn_ks.use)]%>%unname(),rep(mean(ref_weight[rownames(res_df_ks)][setdiff(c(1:nrow(res_df_ks)),which(topn_ord_ks<=topn_ks.use))]),topn_ks_last%>%length()))
                topn_ord_ks<-topn_ord_ks[which(topn_ord_ks<=topn_ks.use)]
                r1<-ks.test.truncated(c(topn_ord_ks,topn_ks_last),w.x=w.x,y='punif',min=0,max=nrow(res_df),alternative = 'gr',y.min=0, y.max=nrow(res_df))
              }else{
                if(sum(is.na(res_df_ks[,meth.comb]))>0){
                  print('in here')
                  topn_ks.use=sum(!is.na(res_df[,meth.comb]))
                  topn_ord_ks<-res_df_ks[,meth.comb]
                  topn_ks_last<-topn_ks.use+(1:(length(topn_ord_ks)-sum(topn_ord_ks<=topn_ks.use, na.rm = T)))*(nrow(res_df)-topn_ks.use)/(length(topn_ord_ks)-sum(topn_ord_ks<=topn_ks.use, na.rm = T))
                  w.x=c(ref_weight[rownames(res_df_ks)][which(topn_ord_ks<=topn_ks.use)]%>%unname(),rep(mean(ref_weight[rownames(res_df_ks)][setdiff(c(1:nrow(res_df_ks)),which(topn_ord_ks<=topn_ks.use))]),topn_ks_last%>%length()))
                  topn_ord_ks<-topn_ord_ks[which(topn_ord_ks<=topn_ks.use)]
                  r1<-ks.test.truncated(c(topn_ord_ks,topn_ks_last),w.x=w.x,y='punif',min=0,max=nrow(res_df),alternative = 'gr',y.min=0, y.max=nrow(res_df))
                }else{
                  #just run full ks test
                  w.x=ref_weight[rownames(res_df_ks)]
                  r1<-ks.test.truncated(res_df_ks[,meth.comb],w.x=w.x,y='punif',min=0,max=nrow(res_df),alternative = 'gr',y.min=0, y.max=nrow(res_df))
                  # r1<-ks.test(res_df_ks[,meth.comb],y='punif',min=0,max=nrow(res_df),alternative = 'gr')
                }
                
                
              }
              print(r1$p.value)
              final_df_ks[i,1]=r1$p.value
            }
            write.table(final_df_ks,paste0(saved_dir,ref,'.',up.genes,' | ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),'gene threshold : expressed in  ',filter.rate*100,'percent cells | adj.p_threshold : ',adj.p_thr, '|',sort.meth, ' sorted.',ifelse(is.null(topn_ks),'',topn_ks.write),'weighted_ks.txt'),sep = '\t')
            
            rm(res_list,res_df,res_df_ks,final_df_ks)
          }
        }
      }
      
    }
  }
}

