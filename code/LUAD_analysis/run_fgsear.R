library(fgsea)
library(dplyr)
library(magrittr)
library(readr)
library(stringr)
`%ni%`<-Negate(`%in%`)
base_dir='~/'
source(paste0(base_dir,'script_bk/visualize.sources.R'))


meth.comb.selection=list()
meth.comb.selection[['Epithelial cells']]=c('Raw_Wilcox','MAST','DESeq2','ZW_DESeq2_Cov','edgeR_Cov','ZW_edgeR_Cov','limmatrend','MNNCorrect_limmatrend')
surfix='wikipath'
ct='Epithelial cells'
dir.create(paste0(base_dir,'data/gsea/'),recursive = T,showWarnings = F)
for(filter.rate in c(0.05)){
  for(sort.meth in c('pvalue')){
    for(gene.used in c('scLUAD-TCGA')){
      if(gene.used=='scLUAD-TCGA'){
        meth.anals=c(meth.comb.selection[[ct]],'TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend')
      }else{
        meth.anals=meth.comb.selection[[ct]]
      }
      for(meth.anal in meth.anals){
        res.fgsear<-list()
        for(ctsubct in c('E')){
          {
            ctsubct.list<-get.ctsubct(ctsubct)
            ct=ctsubct.list$ct
            sub.ct=ctsubct.list$sub.ct
            filt_dir=ifelse(is.null(filter.rate),'',paste0('/filtered/',filter.rate,ifelse(is.null(sub.ct),'',paste0('_',sub.ct))))
            dat_dir.saved=paste0(base_dir,'data/analyzed/',gene.used,'/',ct,ifelse(is.null(sub.ct),'',paste0('_',sub.ct)),'/',sort.meth,'/')
            genes_dir=paste0(base_dir,'data',filt_dir,'/third_processed/')
          }
          {
            load(paste0(dat_dir.saved,'both | ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),'gene threshold : expressed in  ',filter.rate*100,'percent cells | adj.p_threshold : none|',sort.meth,' sorted.ks.RData'))
            
            if(sort.meth=='logFC'){
              res_df<-sapply(setdiff(names(res_list),'genes'),FUN = function(x){res_list[[x]][['log2FC']][match(res_list[['genes']],res_list[[x]][['genes']])]})
              
              colnames(res_df)<-sapply(colnames(res_df),FUN=function(x){select_bk_df(x)})%>%as.vector()
              if(any(setdiff(meth.comb.selection[[ct]],c('Pseudobulk','TCGA','GSE'))%ni%colnames(res_df))){
                stop('wrong columns in res_df')
              }
              res_df<-res_df%>%as.data.frame()
              # res_df$MAST=res_df$MAST*-1
              # res_df$MAST_Cov=res_df$MAST_Cov*-1
              
              res_df<-res_df%>%as.data.frame()%>%dplyr::select(meth.anal)
              rownames(res_df)=res_list[['genes']]
              rank_mean<-apply(res_df,1,mean)
            }else if(sort.meth=='pvalue'){
              res_df<-sapply(setdiff(names(res_list),'genes'),FUN = function(x){match(res_list[['genes']],res_list[[x]][['genes']])*sign(res_list[[x]][['log2FC']])[match(res_list[['genes']],res_list[[x]][['genes']])]})
              colnames(res_df)<-sapply(colnames(res_df),FUN=function(x){select_bk_df(x)})%>%as.vector()
              if(any(setdiff(meth.comb.selection[[ct]],c('Pseudobulk','TCGA','GSE'))%ni%colnames(res_df))){
                stop('wrong columns in res_df')
              }
              
              res_df<-res_df%>%as.data.frame()
              # res_df$MAST=res_df$MAST*-1
              # res_df$MAST_Cov=res_df$MAST_Cov*-1
              res_df<-res_df%>%dplyr::select(meth.anal)
              rownames(res_df)=res_list[['genes']]
              for(k in 1:ncol(res_df)){
                k1<-res_df[,k]
                k1.pos<-rank(k1[k1>=0]*(-1))
                k1.neg<-rank(k1[k1<0])*(-1)
                k1[k1>=0]=k1.pos
                k1[k1<0]=k1.neg
                res_df[,k]=k1
              }
              rank_mean<-apply(res_df,1,mean)
            }
            res_df%<>%as.data.frame()
            rank_mean%<>%sort()
            dbs.use<-c('WikiPathway_2021_Human')
            res.temp<-list()
            for(dbtb in dbs.use){
              dbf<-read_lines(file=paste0(base_dir,'data/original/',dbtb,'.txt'))
              dbl<-sapply(1:length(dbf),FUN=function(x){
                rt<-str_split(dbf[x],pattern='\t')[[1]]
              })
              names(dbl)=sapply(dbl,FUN=function(x){x[1]})
              dbl2<-sapply(names(dbl),FUN=function(x){
                setdiff(dbl[[x]],c(x,''))
              })
              res2<-fgsea(dbl2,(rank_mean), nPermSimple=100000,scoreType='std', minSize=15, maxSize=500, nproc=8)
              res2%<>%arrange(pval)
              res.temp[[dbtb]]=res2
            }
            
          }
          res.fgsear[[paste0(ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)))]]=res.temp
          
          meta.data=list()
          meta.data[[paste0(ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),'_genes')]]=res_list$genes
        }
        save(res.fgsear,meta.data,file=paste0(base_dir,'data/gsea/',filter.rate,"_filtered_gsea",'_',sort.meth,'_sorted_',gene.used,'_',surfix,meth.anal,'.RData'))
      }
    }
  }
}

surfix='wikipath'
meth.anals=c('TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend')
res_list.tot<-list()
for(sort.meth in c('logFC','pvalue')){
  res_list=list()
  for(meth.anal in meth.anals){
    meth.anal.use<-c('DESeq2_Cov','edgeR_Cov','voom.limma_Cov','limma_trend_Cov')[meth.anal==c('TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend')]
    result.table<-read.table(paste0(base_dir,'data/analyzed/TCGA-LUAD.rnaseq.processed.htseq_count.',meth.anal.use,'.txt'))
    
    if(sort.meth=='logFC'){
      result.table%<>%dplyr::arrange(desc(abs(logFC)))
      result.table$sort=abs(result.table$logFC)
    }else if(sort.meth=='pvalue'){
      result.table$sort=result.table$pval
      abs.logfc<-abs(result.table$logFC)
      result.table$sort<-sapply(result.table$sort,FUN=function(x)which(x==(result.table$sort%>%unique()%>%sort(decreasing = T))))+1
      
      result.table$sort2=result.table$sort+abs.logfc/ifelse(max(abs.logfc)<1,1,max(abs.logfc))
      
      result.table%<>%dplyr::arrange(desc(sort2))
    }
    ranked_genes<-rownames(result.table)
    res_list[[meth.anal]][['genes']]=ranked_genes
    res_list[[meth.anal]][['pval']]=result.table$pval
    res_list[[meth.anal]][['adj.pval']]=result.table$adj.pval
    res_list[[meth.anal]][['logFC']]=result.table$logFC
  }
  res_list.tot[[sort.meth]]=res_list
}
save(res_list.tot,file=paste0(base_dir,'data/analyzed/TCGA-LUAD.rnaseq.processed.htseq_count.analyzed.RData'))

gene.used='TCGA'
for(sort.meth in c('logFC','pvalue')){
  res_list=res_list.tot[[sort.meth]]
  
  for(meth.anal in meth.anals){
    res.fgsear<-list()
    {
      if(sort.meth=='logFC'){
        res_df<-sapply(names(res_list),FUN = function(x){res_list[[x]][['logFC']][match(res_list[[1]][['genes']],res_list[[x]][['genes']])]})
        res_df<-res_df%>%as.data.frame()%>%dplyr::select(meth.anal)
        rownames(res_df)=res_list[[1]][['genes']]
        rank_mean<-apply(res_df,1,mean)
      }else if(sort.meth=='pvalue'){
        res_df<-sapply(names(res_list),FUN = function(x){match(res_list[[1]][['genes']],res_list[[x]][['genes']])*sign(res_list[[x]][['logFC']])[match(res_list[[1]][['genes']],res_list[[x]][['genes']])]})
        
        res_df<-res_df%>%as.data.frame()%>%dplyr::select(meth.anals)
        rownames(res_df)=res_list[[1]][['genes']]
        for(k in 1:ncol(res_df)){
          k1<-res_df[,k]
          k1.pos<-rank(k1[k1>=0]*(-1))
          k1.neg<-rank(k1[k1<0])*(-1)
          k1[k1>=0]=k1.pos
          k1[k1<0]=k1.neg
          res_df[,k]=k1
        }
        rank_mean<-apply(res_df,1,mean)
      }
      res_df%<>%as.data.frame()
      rank_mean%<>%sort()
      dbs.use<-c('WikiPathway_2021_Human')
      res.temp<-list()
      for(dbtb in dbs.use){
        dbf<-read_lines(file=paste0(base_dir,'data/original/',dbtb,'.txt'))
        dbl<-sapply(1:length(dbf),FUN=function(x){
          rt<-str_split(dbf[x],pattern='\t')[[1]]
        })
        names(dbl)=sapply(dbl,FUN=function(x){x[1]})
        dbl2<-sapply(names(dbl),FUN=function(x){
          setdiff(dbl[[x]],c(x,''))
        })
        res2<-fgsea(dbl2,(rank_mean), nPermSimple=100000,scoreType='std', minSize=15, maxSize=500, nproc=8)
        res2%<>%arrange(pval)
        res.temp[[dbtb]]=res2
      }
      
    }
    res.fgsear[[gene.used]]=res.temp
    
    meta.data=list()
    meta.data[[paste0(gene.used,'_genes')]]=res_list[[1]][['genes']]
    
    save(res.fgsear,meta.data,file=paste0(base_dir,'data/gsea/',"gsea",'_',sort.meth,'_sorted_',gene.used,'_',surfix,meth.anal,'.RData'))
  }
}

surfix='wikipath'
for(gene.used in c('TCGA')){
  for(sort.meth in c('logFC','pvalue')){
    dir.create(paste0(base_dir,'data/gsea/table/',gene.used,'/',sort.meth,'/'),showWarnings = F, recursive = T)
    for(meth.anal in meth.anals){
      load(paste0(base_dir,'data/gsea/',"gsea",'_',sort.meth,'_sorted_',gene.used,'_',surfix,meth.anal,'.RData'))
      for(n in names(res.fgsear)){
        for(dbselect in c('WikiPathway_2021_Human')){
          dt<-res.fgsear[[n]][[dbselect]]
          dt%<>%as.data.frame()%>%arrange(pval)
          dt$leadingEdge<-sapply(dt$leadingEdge,FUN=function(x){paste0(x,collapse ='|')})
          
          write.table(dt, file=paste0(base_dir,'data/gsea/table/',gene.used,'/',sort.meth,'/',n,dbselect,'_',sort.meth,'_sorted_gsea_edge',surfix,meth.anal,'.txt'),sep='\t')
        }
      }
    }
  }
}



meth.comb.selection=list()
meth.comb.selection[['Epithelial cells']]=c('Raw_Wilcox','MAST','DESeq2','ZW_DESeq2_Cov','edgeR_Cov','ZW_edgeR_Cov','limmatrend','MNNCorrect_limmatrend')
surfix='wikipath'
ct='Epithelial cells'

for(filter.rate in c(0.05)){
  for(sort.meth in c('pvalue','logFC')){
    for(gene.used in c('scLUAD-TCGA')){
      ct='Epithelial cells'
      if(gene.used=='scLUAD-TCGA'){
        meth.anals=c(meth.comb.selection[[ct]],'TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend')
      }else{
        meth.anals=meth.comb.selection[[ct]]
      }
      for(meth.anal in meth.anals){
        load(paste0(base_dir,'data/gsea/',filter.rate,"_filtered_gsea",'_',sort.meth,'_sorted_',gene.used,'_',surfix,meth.anal,'.RData'))

        for(n in names(res.fgsear)){
          for(dbselect in c('WikiPathway_2021_Human')){
            dt<-res.fgsear[[n]][[dbselect]]
            dt%<>%as.data.frame()%>%arrange(pval)
            dt$leadingEdge<-sapply(dt$leadingEdge,FUN=function(x){paste0(x,collapse ='|')})
            
            write.table(dt, file=paste0(base_dir,'data/gsea/table/',gene.used,'/',sort.meth,'/',n,dbselect,'_',sort.meth,'_sorted_gsea_edge',surfix,meth.anal,'.txt'),sep='\t')
          }
        }
      }
    }
  }
}
