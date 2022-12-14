#raw count as input count=rawcount
run_pseudobulk<-function(count,cellinfo,former.meth=''){
  library(magrittr)
  library(dplyr)
  library(tidyverse)
  library(stringr)
  rownames(cellinfo)=cellinfo$Cell
  cellinfo<-cellinfo[colnames(count),]
  cellinfo$Group%<>%factor()
  cellinfo$Batch%<>%factor()
  
  count_matrix<-count
  Batch_uniq<-unique(cellinfo$Batch)
  ind_list<-list()
  for(b in Batch_uniq){
    cellinfo_b<-cellinfo%>%dplyr::filter(Batch==b)%>%dplyr::filter(Cell %in% colnames(count_matrix))
    for(g in unique(cellinfo$Group)){
      if(length(cellinfo_b$Cell[which(cellinfo_b$Group==g)])>0){
        ind_list[[paste0(g,'_',b)]]=cellinfo_b$Cell[which(cellinfo_b$Group==g)]
      }
    }
  }
  pseudobulk<-matrix(NA,nrow=nrow(count_matrix),ncol=length(ind_list),dimnames = list(row_names=rownames(count_matrix), col_names=ind_list%>%names()%>%sort()))
  for(i in names(ind_list)%>%sort()){
    pseudobulk[,i]=apply(count_matrix[,which(colnames(count_matrix)%in%ind_list[[i]])],1,sum)
  }
  
  cellinfo=data.frame(Cell=colnames(pseudobulk),Group=gsub(colnames(pseudobulk),pattern = '[_].*$',replacement = ''),Batch=rep('Batch0',ncol(pseudobulk)),row.names=colnames(pseudobulk))
  res<-pseudobulk
  processed<-pseudobulk
  
  res_name<-paste0(ifelse(former.meth=='','',paste0(former.meth,'+')),'pseudobulk')
  save(res, processed, cellinfo, file=paste0('./',res_name,'.rda'))
  return(res_name)
}