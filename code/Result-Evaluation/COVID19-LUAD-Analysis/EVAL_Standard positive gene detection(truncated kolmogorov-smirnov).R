analysis_script_dir='~/'
source(paste0(analysis_script_dir,'KS_truncated.R'))
library(dplyr)
library(tidyverse)
library(magrittr)
meth.comb.list<-get.methcomblist()

for(fdr_cut in c(0.01,0.05,'none')){
  {
    `%ni%`<-Negate(`%in%`)
    ###ks params###
    topn_ks=0.2 #top 20% truncated ks test
  }
  for(ref in c('Known disease genes','Prognostic genes')){
    if(ref=='Known disease genes'){
      source(paste0(analysis_script_dir,'/Known disease gene selection.R'))
      ref_genes<-Known_disease_genes
      ref_weight<-Known_disease_genes_weight
    }else if(ref=='Prognostic genes'){
      # Need to run 'Prognostic gene selection(survival_analysis).R' first
      meta.res<-read.table('survival_analysis_result.txt')
      Prognostic_genes<-rownames(meta.res)[meta.res$FDR<0.05]
      Prognostic_genes_weight=rep(1,length(Prognostic_genes))
      ref_genes<-Prognostic_genes
      ref_weight<-Prognostic_genes_weight
    }else if(ref=='Standard Positive genes(GO:0051607)'){
      source(paste0(analysis_script_dir,'/Standard positive gene selection(GO:0051607).R'))
      ref_genes<-Standard_Positive_gene_COVID19
      ref_weight<-Standard_Positive_gene_COVID19_weight
    }
    for(sort.meth in c('pvalue')){
      load(paste0('DE_result_',sort.meth,'_sorted_fdr_',fdr_cut,'_cut.rda'))
      res_df<-sapply(setdiff(names(result.list),'gene'),FUN=function(x){result.list[[x]][match(result.list[['gene']],rownames(result.list[[x]])),'ranks']})
      rownames(res_df)=res_list[['gene']]
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
          
          
          ##1108 fix
          ref_weight.use<-ref_weight
          topn_ord_ks_use<-topn_ord_ks[which(topn_ord_ks<=topn_ks.use)]
          ref_genes.use<-rownames(res_df_ks)[which(topn_ord_ks<=topn_ks.use)]
          if(any(table(topn_ord_ks_use)>1)){
            for(tie_rank in names(table(topn_ord_ks_use))[table(topn_ord_ks_use)>1]){
              ref_weight.use[ref_genes.use[which(topn_ord_ks_use==tie_rank)]]=mean(ref_weight.use[ref_genes.use[which(topn_ord_ks_use==tie_rank)]])
            }
            for(tie_rank in names(table(topn_ord_ks_use))[table(topn_ord_ks_use)>1]){
              sorted_toku<-sort(unique(topn_ord_ks_use))
              
              if(as.numeric(tie_rank)==min(topn_ord_ks_use)){
                topn_ord_ks_use[topn_ord_ks_use==tie_rank]=((1:sum(topn_ord_ks_use==tie_rank))/sum(topn_ord_ks_use==tie_rank))*as.numeric(tie_rank)+0
              }else{
                topn_ord_ks_use[topn_ord_ks_use==tie_rank]=((1:sum(topn_ord_ks_use==tie_rank))/sum(topn_ord_ks_use==tie_rank))*(as.numeric(tie_rank)-sorted_toku[which(sorted_toku==tie_rank)-1])+sorted_toku[which(sorted_toku==tie_rank)-1]
              }
            }
          }
          w.x=c(ref_weight.use[ref_genes.use]%>%unname(),rep(mean(ref_weight.use[setdiff(rownames(res_df_ks),ref_genes.use)]),topn_ks_last%>%length()))
          r1<-ks.test.truncated(c(topn_ord_ks_use,topn_ks_last),w.x=w.x,y='punif',min=0,max=nrow(res_df),alternative = 'gr',y.min=0, y.max=nrow(res_df))
        }else{
          if(sum(is.na(res_df_ks[,meth.comb]))>0){
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
          }
        }
        final_df_ks[i,1]=r1$p.value
      }
      
      write.table(final_df_ks,paste0('DE_result_',ref,'detection_weighted_ks_',sort.meth,'_sorted_top_',topn_ks*100,'percent_genes_fdr_',fdr_cut,'_cut.txt'),sep='\t')
    }
  }
}


