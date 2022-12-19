# 'TCGA_DESeq2','TCGA_edgeR','TCGA_limmavoom','TCGA_limmatrend', 'GSE43458.limma_Cov','GSE31210.limma_Cov'
# ,'RISC+limmatrend','RISC+Wilcox','RISC+QP',"Scanorama+limmatrend", "scGen+limmatrend", "scVI+limmatrend", "Scanorama+Wilcox", "scGen+Wilcox", "scVI+Wilcox"
setwd('/hdd2/SC lung/run_example/')

{
  library(stringr)
  library(magrittr)
  library(tidyverse)
  library(dplyr)
  #source all run_*.R in ref_script
  dir_refscript='/hdd2/SC lung/ref_script/'
  files.sources = list.files(dir_refscript,full.names = T)
  sapply(files.sources, source)
  
  #generate counts and cellinfo data. If already exist, not necessary.
  load('/hdd2/SC lung/data/filtered/0.05/first_processed/T-NK cells.0.05mincell_filtered_raw_first_processed.RData')
  write.table(first_processed$count,'counts.txt')
  write.table(data.frame(Cell=colnames(first_processed$count),Group=first_processed$group, Batch=first_processed$batch,row.names = colnames(first_processed$count)),'cellinfo.txt')
}
for(i in c('Combat+Wilcox','limma_BEC+Wilcox','MNN+Wilcox','scMerge+Wilcox','Raw_Wilcox','ZW_BEC+Wilcox','Seurat+Wilcox','Pseudobulk+DESeq2','Pseudobulk+edgeR','Pseudobulk+limmavoom','Pseudobulk+limmatrend','MAST','MAST_Cov','DESeq2','DESeq2_Cov','ZW+DESeq2','ZW+DESeq2_Cov','DESeq2_sep+wFisher','edgeR','edgeR_Cov','edgeR_DetRate','edgeR_DetRate_Cov','ZW+edgeR','ZW+edgeR_Cov','edgeR_sep+wFisher','limmavoom','limmavoom_Cov','limmatrend','limmatrend_Cov','Combat+limmatrend','MNN+limmatrend','scMerge+limmatrend','LogN_sep+limmatrend_sep+wFisher','DESeq2_sep+FEM','LogN_sep+ES_sep+FEM','DESeq2_sep+REM','LogN_sep+ES_sep+REM')){
# for(i in c('MNN+Wilcox')){
  
  methods<-str_split(i,pattern='[+]')[[1]]
  print(i)
  former.meth=''
  for(m in methods){
    former.meth=try(do.run(m,former.meth = former.meth,P_ref.ind=NULL))
    print(former.meth)
  }
}

do.run<-function(m,former.meth='',P_ref.ind=NULL){
  if(former.meth==''){
    count=read.table('counts.txt')
    cellinfo=read.table('cellinfo.txt')
    processed<-list()
    for(b in unique(cellinfo$Batch)){
      processed[[b]]=count[,cellinfo$Cell[which(cellinfo$Batch==b)]]%>%as.matrix()
    }
    if(m=='Combat'){
      post.meth=run_combat(count,cellinfo)
    }else if(m=='limma_BEC'){
      post.meth=run_limma_bec(count,cellinfo)
    }else if(m=='MNN'){
      post.meth=run_mnn(count,cellinfo)
    }else if(m=='limma_BEC'){
      post.meth=run_limma_bec(count,cellinfo)
    }else if(m=='scMerge'){
      post.meth=run_scMerge(count,cellinfo)
    }else if(m%in%c('ZW','ZW_BEC')){
      post.meth=run_zinbwave(count,cellinfo)
    }else if(m=='Seurat'){
      post.meth=run_Seurat(count,cellinfo)
    }else if(m=='Pseudobulk'){
      post.meth=run_pseudobulk(count,cellinfo)
    }else if(m=='RISC'){
      post.meth=run_RISC(count,cellinfo,P_ref.ind=P_ref.ind)
    }else if(m=='Raw_Wilcox'){
      post.meth=run_wilcox(processed = count,cellinfo,is.log = F)
    }else if(m%in%c('MAST','MAST_Cov')){
      post.meth=run_MAST(processed = count,cellinfo,cov=str_detect(m,pattern='Cov$'))
    }else if(m%in%c("DESeq2","DESeq2_Cov")){
      post.meth=run_DESeq2(processed = count,cellinfo,cov = str_detect(m,pattern='Cov$'))
    }else if(m%in%c("edgeR","edgeR_Cov","edgeR_DetRate","edgeR_DetRate_Cov")){
      post.meth=run_edgeR(processed = count,cellinfo,cov=str_detect(m,pattern='Cov$'),Det = str_detect(m,pattern='DetRate'))
    }else if(m%in%c("limmavoom","limmavoom_Cov")){
      post.meth=run_limmavoom(processed = count,cellinfo,cov=str_detect(m,pattern='Cov$'))
    }else if(m %in% c('limmatrend','limmatrend_Cov')){
      post.meth=run_limmatrend(processed = count,cellinfo,cov=str_detect(m,pattern='Cov$'))
    }else if(m=='LogN_sep'){
      post.meth=run_LogNormalize(count = count,cellinfo,separate = T)
    }else if(m=='DESeq2_sep'){
      post.meth=run_DESeq2_sep(processed = processed,cellinfo)
    }else if(m=='edgeR_sep'){
      post.meth=run_edgeR_sep(processed = processed,cellinfo)
    }
  }else{
    load(paste0(former.meth,'.rda'))
    if(m=="limmatrend_sep"){
      post.meth=run_limmatrend_sep(processed = processed,cellinfo,former.meth=former.meth)
    }else if(m=="ES_sep"){
      post.meth=run_ES(processed = processed,cellinfo,former.meth=former.meth)
    }else if(m=="limmatrend"){
      if(str_detect(former.meth,pattern = 'pseudobulk')){
        post.meth=run_limmatrend(processed = processed,cellinfo,cov=str_detect(m,pattern='Cov$'),former.meth = former.meth)
      }else{
        post.meth=run_limmatrend_BECdata(processed = processed,cellinfo,former.meth=former.meth)
      }
    }else if(m%in%c('DESeq2','DESeq2_Cov')){
      if(str_detect(former.meth,pattern = 'pseudobulk')){
        post.meth=run_DESeq2(processed = processed,cellinfo,cov = str_detect(m,pattern='Cov$'),former.meth=former.meth)
      }else{
        post.meth=run_DESeq2_zinbwavedata(processed = res,cellinfo,cov = str_detect(m,pattern='Cov$'),former.meth=former.meth)
      }
    }else if(m%in%c('edgeR','edgeR_Cov')){
      if(str_detect(former.meth,pattern = 'pseudobulk')){
        post.meth=run_edgeR(processed = processed,cellinfo,cov=str_detect(m,pattern='Cov$'),Det = str_detect(m,pattern='DetRate'),former.meth=former.meth)
      }else{
        post.meth=run_edgeR_zinbwavedata(processed = res,cellinfo,cov=str_detect(m,pattern='Cov$'),former.meth=former.meth)
      }
    }else if(m=='limmavoom'){
      if(str_detect(former.meth,pattern = 'pseudobulk')){
        post.meth=run_limmavoom(processed = processed,cellinfo,cov=str_detect(m,pattern='Cov$'),former.meth = former.meth)
      }
    }else if(m=='QP'){
      post.meth=run_QP(processed = processed,cellinfo,former.meth = former.meth)
    }else if(m=='REM'){
      post.meth=run_REM(res = res,cellinfo,former.meth = former.meth)
    }else if(m=='FEM'){
      post.meth=run_FEM(res = res,cellinfo,former.meth = former.meth)
    }else if(m=='wFisher'){
      post.meth=run_wFisher(res = res,processed=processed,cellinfo,former.meth = former.meth)
    }else if(m=='Wilcox'){
      post.meth=run_wilcox(processed = processed,cellinfo,is.log = T,former.meth=former.meth)
    }
  }
  return(post.meth)
}

# 
# #example
# run_combat(count,cellinfo)
# load('combat.rda')
# run_wilcox(processed = T,cellinfo=cellinfo,is.log=T,former.meth = 'combat')
# load('combat+wilcox.rda')
# 
# #example2
# run_limmavoom(count,cellinfo,cov=T)
# load('limmavoom.rda')
# #example3
# run_LogNormalize(count,cellinfo,separate = T,former.meth = '')
# load('LogNormalize_sep.rda')
# run_limmatrend_sep(processed=processed,cellinfo=cellinfo,former.meth = 'LogNormalize')
# load('LogNormalize_sep+limmatrend_sep.rda')
# run_wFisher(res,processed=processed,cellinfo=cellinfo,former.meth='LogNormalize_sep+limmatrend')
# load('LogNormalize_sep+limmatrend_sep+wfisher.rda')
# 
# 
# count<-first_processed$count
# cellinfo<-data.frame(Cell=colnames(count),Batch=first_processed$batch,Group=first_processed$group,row.names = colnames(count))
# run_limmavoom(count,cellinfo,cov = T)
