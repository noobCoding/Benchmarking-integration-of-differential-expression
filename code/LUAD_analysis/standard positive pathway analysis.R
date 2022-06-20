source(paste0(base_dir,'script_bk/visualize.sources.R'))
library(magrittr)
library(tidyverse)
meth.comb.selection=list()
meth.comb.selection[['Epithelial cells']]=c('Raw_Wilcox','MAST','DESeq2','ZW_DESeq2_Cov','edgeR_Cov','ZW_edgeR_Cov','limmatrend','MNNCorrect_limmatrend')
surfix='wikipath'


gold_ref.df<-read.csv(paste0(base_dir,'data/original/Gold standard pathway_epi_analysis.csv'))
gold_ref.df$class%<>%as.numeric()

gold_ref.df$class[which(gold_ref.df$class!=3)]=NA
gold_ref.df$category[which(gold_ref.df$class!=3)]=''



result.list=list()
meth.comb.selection=list()
meth.comb.selection[['Epithelial cells']]=c('Raw_Wilcox','MAST','DESeq2','ZW_DESeq2_Cov','edgeR_Cov','ZW_edgeR_Cov','limmatrend','MNNCorrect_limmatrend')
surfix='wikipath'
ct='Epithelial cells'

for(sort.meth in c('pvalue','logFC')){
  result.list[[sort.meth]]=list()
  for(gene.used in c('scLUAD-TCGA')){
    result.list[[sort.meth]][[gene.used]]=list()
    dir.create(paste0(base_dir,'data/gsea/',gene.used,'/'),recursive = T,showWarnings = F)
    if(gene.used=='scLUAD-TCGA'){
      meth.anals=c(meth.comb.selection[[ct]],'TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend')
      
      filenames.base=paste0(base_dir,'data/gsea/0.05_filtered_gsea_',sort.meth,'_sorted_',gene.used,'_',surfix)
    }else if(gene.used=='scLUAD'){
      meth.anals=c(meth.comb.selection[[ct]])
      filenames.base=paste0(base_dir,'data/gsea/0.05_filtered_gsea_',sort.meth,'_sorted_',gene.used,'_',surfix)
    }else if(gene.used=='TCGA'){
      meth.anals=c('TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend')
      filenames.base=paste0(base_dir,'data/gsea/gsea_',sort.meth,'_sorted_',gene.used,'_',surfix)
    }
    for(meth.anal in meth.anals){
      result.list[[sort.meth]][[gene.used]][[meth.anal]]=list()
      load(paste0(filenames.base,meth.anal,'.RData'))
      for(n in names(res.fgsear)){
        for(dbselect in c('WikiPathway_2021_Human')){
          dt<-res.fgsear[[n]][[dbselect]]
          dt%<>%as.data.frame()%>%arrange(pval)
          # dt<-dt[1:100,c('pathway','pval','padj','ES')]
         
          df_add<-data.frame(category=gold_ref.df$category[match(dt$pathway,gold_ref.df$pathway)],class=gold_ref.df$class[match(dt$pathway,gold_ref.df$pathway)])
          dt_res<-cbind(dt,df_add)
          dt_res<-dt_res[,c('pathway','category','class','padj','pval','ES')]

          write.table(dt_res,paste0(base_dir,'data/gsea/',gene.used,'/',sort.meth,'_sorted_',gene.used,'_',surfix,meth.anal,'_table.txt'),sep='\t')
        }
      }
    }
  }
}



meth.comb.selection=list()
meth.comb.selection[['Epithelial cells']]=c('Raw_Wilcox','MAST','DESeq2','ZW_DESeq2_Cov','edgeR_Cov','ZW_edgeR_Cov','limmatrend','MNNCorrect_limmatrend')
surfix='wikipath'
ct='Epithelial cells'
for(cutoff in c(0.25, 'none')){
  result.list<-list()
  sort.meths=meth.names=N.sig.paths=c()
  for(sort.meth in c('pvalue','logFC')){
    print(sort.meth)
    result.list[[sort.meth]]=list()
    for(gene.used in c('scLUAD-TCGA')){
      result.list[[sort.meth]][[gene.used]]=list()
      if(gene.used=='scLUAD-TCGA'){
        meth.anals=c(meth.comb.selection[[ct]],'TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend')
        filenames.base=paste0(base_dir,'data/gsea/0.05_filtered_gsea_',sort.meth,'_sorted_',gene.used,'_',surfix)
      }else if(gene.used=='scLUAD'){
        meth.anals=c(meth.comb.selection[[ct]])
        filenames.base=paste0(base_dir,'data/gsea/0.05_filtered_gsea_',sort.meth,'_sorted_',gene.used,'_',surfix)
      }else if(gene.used=='TCGA'){
        meth.anals=c('TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend')
        filenames.base=paste0(base_dir,'data/gsea/gsea_',sort.meth,'_sorted_',gene.used,'_',surfix)
      }
      print(gene.used)
      for(meth.anal in meth.anals){
        result.list[[sort.meth]][[gene.used]][[meth.anal]]=list()
        load(paste0(filenames.base,meth.anal,'.RData'))
        for(n in names(res.fgsear)){
          for(dbselect in c('WikiPathway_2021_Human')){
            dt<-res.fgsear[[n]][[dbselect]]
            dt%<>%as.data.frame()%>%arrange(pval)
            if(cutoff=='none'){
            }else{
              dt%<>%dplyr::filter(padj<=cutoff)
            }
            
            if(gene.used!='scLUAD'){
              sort.meths%<>%c(sort.meth)
              if(gene.used=='TCGA'){
                meth.names%<>%c(paste0('whole_',meth.anal))
              }else{
                meth.names%<>%c(meth.anal)
              }
              
              N.sig.paths%<>%c(nrow(dt))
            }
            print(meth.anal)
            print(nrow(dt))
            result.list[[sort.meth]][[gene.used]][[meth.anal]][['pathways']]=dt$pathway
            result.list[[sort.meth]][[gene.used]][[meth.anal]][['pvalue']]=dt$pval
            result.list[[sort.meth]][[gene.used]][[meth.anal]][['qvalue']]=dt$padj
          }
        }
      }
    }
  }
  save(result.list, file=paste0(base_dir,'data/gsea/gsea.list_',cutoff,'cut.',surfix,'.RData'))
  df_nsigpath=data.frame(methods=meth.names,N.sig.paths=N.sig.paths,sort.meth=sort.meths)
  write.table(df_nsigpath%>%dplyr::filter(sort.meth=='pvalue'),file=paste0(base_dir,'data/gsea/N_sig.paths.pvalue_sorted_gsea',cutoff,'cutoff.',surfix,'.txt'),sep='\t')
  write.table(df_nsigpath%>%dplyr::filter(sort.meth=='logFC'),file=paste0(base_dir,'data/gsea/N_sig.paths.logFC_sorted_gsea',cutoff,'cutoff.',surfix,'.txt'),sep='\t')
}


meth.comb.selection2=list()
meth.comb.selection2[['Epithelial cells']]=c('Raw_Wilcox','MAST','DESeq2','ZW_DESeq2_Cov','edgeR_Cov','ZW_edgeR_Cov','limmatrend','MNNCorrect_limmatrend')
surfix='wikipath'




cutoff=0.25
topn_path=50
gold_ref.df$class%<>%as.numeric()
xbase=gold_ref.df%>%nrow()
ybase=sum(gold_ref.df$class==3)
gold_ref.use=gold_ref.df%>%dplyr::filter(class==3)
ref_paths<-gold_ref.use$pathway
for(sort.meth in c('logFC')){
  for(gene.used in c('scLUAD-TCGA')){
    if(gene.used=='scLUAD-TCGA'){
      meth.anals=c(meth.comb.selection[[ct]],'TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend')
      
      filenames.base=paste0(base_dir,'data/gsea/0.05_filtered_gsea_',sort.meth,'_sorted_',gene.used,'_',surfix)
    }else if(gene.used=='scLUAD'){
      meth.anals=c(meth.comb.selection[[ct]])
      filenames.base=paste0(base_dir,'data/gsea/0.05_filtered_gsea_',sort.meth,'_sorted_',gene.used,'_',surfix)
    }
    dir.create(paste0(base_dir,'data/gsea/pathway detection curve/'),recursive = T, showWarnings = F)
    fn=paste0(base_dir,'data/gsea/pathway detection curve/standard positive pathway curves_',gene.used,'_',sort.meth,'sorted gsea_',cutoff,'cutoff')
    fn=paste0(fn,'.',surfix,'.')
    load(paste0(base_dir,'data/gsea/gsea.list_',cutoff,'cut.',surfix,'.RData'))
    draw_plot_paths(fn,gene.used,meth.anal.use=meth.anals,ref_paths,is.TCGA=F,sort.meth, result.list,xmax=topn_path,xbase=xbase,ybase=ybase)
  }
}


std.cancer<-c('cell proliferation/cell fate and differentiation','cell survival','cell metabolism','cell polarity and migration','genomic instability','inflammation','angiogenesis','oxidative stress')
oncosig.cancer<-c('RTK/RAS pathway','WNT pathway','PI3K pathway','HIPPO pathway','NOTCH pathway','TP53 pathway','TGF-beta pathway','MYC pathway')
cancer.term=c('Term tumor/carcinoma/cancer')
meth.anals<-c('Raw_Wilcox','MAST','DESeq2','ZW_DESeq2_Cov','edgeR_Cov','ZW_edgeR_Cov','limmatrend','MNNCorrect_limmatrend','TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend')
rn=c('N_detected',std.cancer, oncosig.cancer,cancer.term)
dm<-matrix(data=NA, nrow=length(rn),ncol=length(meth.anals), dimnames = list(rownames=rn,colnames=meth.anals))


sort.meth="logFC"
surfix='wikipath'
ct='Epithelial cells'
cutoff=0.25
for(topn in c(50)){
  dm<-matrix(data=NA, nrow=length(rn),ncol=length(meth.anals), dimnames = list(rownames=rn,colnames=meth.anals))
  print(sort.meth)
  for(meth.anal in meth.anals){
    print(meth.anal)
    a<-read.table(paste0(base_dir,'data/gsea/',gene.used,'/',sort.meth,'_sorted_',gene.used,'_',surfix,meth.anal,'_table.txt'))
    # a%<>%dplyr::filter(padj<=0.25)
    a%<>%dplyr::arrange(pval)
    if(topn>1){
      a<-a[1:topn,]
    }else{
      a%<>%dplyr::filter(padj<=topn)
    }
    if(cutoff!='none'){
      a%<>%dplyr::filter(padj<=0.25)
    }
    len_a<-nrow(a)
    a%<>%dplyr::filter(class %in% c(3))
    cat_a<-a$category%>%table()%>%sort()
    len_a.1<-sum(cat_a)
    cat_a<-cat_a[cat_a>=1]
    
    divided=cat_a[rn]%>%as.vector()
    divided[is.na(divided)]=0
    divided[1]=sum(divided)
    len_cat_a_1<-gold_ref.df%>%dplyr::filter(class%in%c(3))%>%select(category)%>%table()
    div<-len_cat_a_1[rn]%>%as.vector()
    div[is.na(div)]=0
    div[1]=len_a
    div.vec<-paste0(divided,'/',div)
    
    dm[,meth.anal]=div.vec
    
  }
  write.table(dm,file=paste0(base_dir,'data/gsea/categorized.pathway detection_',sort.meth,'_sorted_top',topn,'_',cutoff,'cutoff_',surfix,'.csv'),sep=',')
}
