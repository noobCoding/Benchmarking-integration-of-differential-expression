library(fgsea)
library(dplyr)
library(magrittr)
library(readr)
library(stringr)
base_dir='set base directory/'##Directory contains scLUAD&TCGA DE analysis results.
setwd(base_dir)
analyzed_data.dir=base_dir ##Can set different directory. 'scLUAD_Epithelial' example folder name. We have performed gsea on epithelial cells of LUAD data only.
standard_positive_dir=base_dir ##Can set different directory. Directory contains standard positive data for real data analysis(Data for selecting known disease genes, survival analysis results for prognostic genes, GO biological pathway data, wikipathway geneset for gsea, scLUAD standard positive pathways ...).

`%ni%`<-Negate(`%in%`)
normalize_p<-function(x){
  #to avoid -Inf value (log2(pvalue))
  x<-x*0.99999+0.00001
  x<-x
  return(x)
}

methods_to_plot<-c('TCGA_DESeq2','TCGA_edgeR','TCGA_limmavoom','TCGA_limmatrend',
                       'DESeq2_FEM','DESeq2', 'ZW_DESeq2','LogN_FEM','edgeR_Cov',
                       'ZW_edgeR_Cov','limmatrend_Cov','RISC_QP','MAST_Cov','Raw_Wilcox')



{
  analyze_with='TCGA'
  fdr_cut='none'
  gsea_dat='wikipath'
  dbs.use<-c('WikiPathway_2021_Human') #select geneset for gsea.
}
for(sort.meth in c('pvalue','logFC')){
  load(paste0(analyzed_data.dir,'/DE_result_analyzed_with_',analyze_with,'_',sort.meth,'_sorted_fdr_',fdr_cut,'_cut.rda'))
  res.fgsear=list()
  geneNames=result.list$gene
  if(sort.meth=='logFC'){
    
    res_df<-sapply(setdiff(names(result.list),'gene'),FUN=function(x){result.list[[x]][match(result.list[['gene']],rownames(result.list[[x]])),'logFC']})
    res_df<-res_df%>%as.data.frame()
    rownames(res_df)=result.list[['gene']]
    if(anyNA(res_df)){
      res_df[which(is.na(res_df),arr.ind = T)]=0
    }
  }else if(sort.meth=='pvalue'){
    res_df<-sapply(setdiff(names(result.list),'gene'),FUN=function(x){
      -log2(normalize_p(result.list[[x]][match(result.list[['gene']],rownames(result.list[[x]])),'pval']))*sign(result.list[[x]][match(result.list[['gene']],rownames(result.list[[x]])),'logFC'])
    })
    res_df<-res_df%>%as.data.frame()
    rownames(res_df)=result.list[['gene']]
    if(anyNA(res_df)){
      res_df[which(is.na(res_df),arr.ind = T)]=0
    }
  }
  for(m in methods_to_plot){
    res_df.temp<-res_df%>%dplyr::select(m)
    rank_mean<-apply(res_df.temp,1,mean)
    rank_mean%<>%sort()
    res.temp<-list()
    for(dbtb in dbs.use){
      dbf<-read_lines(file=paste0(standard_positive_dir,dbtb,'.txt'))
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
    res.fgsear[[m]]=res.temp
  }
  save(res.fgsear,geneNames,file=paste0(analyzed_data.dir,'/gsea','_',sort.meth,'_sorted_',gsea_dat,'.RData'))
}





gold_ref.df<-read.csv(paste0(standard_positive_dir,'LUAD_standard_positive_pathways(wikipathway).csv'))
std.cancer<-c('cell proliferation/cell fate and differentiation','cell survival','cell metabolism','cell polarity and migration','genomic instability','inflammation','angiogenesis','oxidative stress')
oncosig.cancer<-c('RTK/RAS pathway','WNT pathway','PI3K pathway','HIPPO pathway','NOTCH pathway','TP53 pathway','TGF-beta pathway','MYC pathway')
cancer.term=c('Term tumor/carcinoma/cancer')

rn=c('N_detected',std.cancer, oncosig.cancer,cancer.term)
dm<-matrix(data=NA, nrow=length(rn),ncol=length(methods_to_plot), dimnames = list(rownames=rn,colnames=methods_to_plot))

gene.used='1202scLUAD-TCGA'
sort.meth="pvalue"
# sort.meth="logFC"
dbtb='WikiPathway_2021_Human'
cutoff=0.25
for(topn in c(50)){
  dm<-matrix(data=NA, nrow=length(rn),ncol=length(methods_to_plot), dimnames = list(rownames=rn,colnames=methods_to_plot))
  load(paste0(analyzed_data.dir,'/gsea','_',sort.meth,'_sorted_',gsea_dat,meth.anal,'.RData'))
  print(sort.meth)
  
  for(m in methods_to_plot){
    print(m)
    dt<-res.fgsear[[m]][[dbtb]]
    dt%<>%as.data.frame()%>%arrange(pval)
    dt$leadingEdge<-sapply(dt$leadingEdge,FUN=function(x){paste0(x,collapse ='|')})
    category=gold_ref.df$category[match(dt$pathway,gold_ref.df$pathway)]
    dt<-cbind(dt,category)
    dt<-dt[,c('pathway','category','padj','pval','ES')]
    
    dt<-dt[1:topn,]
    if(cutoff!='none'){
      dt%<>%dplyr::filter(padj<=cutoff)
    }
    len_dt<-nrow(dt)
    cat_dt<-dt$category%>%table()%>%sort()
    len_dt.1<-sum(cat_dt)
    cat_dt<-cat_dt[cat_dt>=1]
    
    divided=cat_dt[rn]%>%as.vector()
    divided[is.na(divided)]=0
    divided[1]=sum(divided)
    len_cat_dt_1<-gold_ref.df%>%select(category)%>%table()
    div<-len_cat_dt_1[rn]%>%as.vector()
    div[is.na(div)]=0
    div[1]=len_dt
    div.vec<-paste0(divided,'/',div)
    
    dm[,m]=div.vec
    
  }
  write.table(dm,file=paste0(base_dir,'data/simulation/analyzed/1202gsea/categorized.pathway detection_',sort.meth,'_sorted_top',topn,'_',cutoff,'cutoff_',surfix,'.csv'),sep=',')
}