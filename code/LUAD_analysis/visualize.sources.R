select_bk_df = function(x){
  Methods_bk=c('base_line','TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend','Combat_Wilcox','limma_BEC_Wilcox','MNNCorrect_Wilcox','scMerge_Wilcox','Raw_Wilcox','ZW_Wilcox','Seurat_Wilcox','Pseudobulk_DESeq2','Pseudobulk_edgeR','Pseudobulk_limma','Pseudobulk_limma_trend','MAST','MAST_Cov','DESeq2','DESeq2_Cov','ZW_DESeq2','ZW_DESeq2_Cov','DESeq2_Fisher','DESeq2_wFisher','edgeR','edgeR_Cov','edgeR_DetRate','edgeR_DetRate_Cov','ZW_edgeR','ZW_edgeR_Cov','edgeR_Fisher','edgeR_wFisher','limma','limma_Cov','limmatrend','limmatrend_Cov','Combat_limmatrend','MNNCorrect_limmatrend','scMerge_limmatrend','LogNorm+limmatrend_Fisher','LogNorm+limmatrend_wFisher','DESeq2_FEM','voom_FEM','LogNorm_FEM','DESeq2_REM','voom_REM','LogNorm_REM','voom+modt_Fisher','voom+modt_wFisher', 'GSE43458.limma_Cov','GSE10072.limma_Cov','GSE31210.limma_Cov')
  
  meth_raw<-c('base_line','raw++TCGA_DESeq2_Cov','raw++TCGA_edgeR_Cov','raw++TCGA_voom.limma_Cov','raw++TCGA_limma_trend_Cov','combat++findmarkers','limma_bec++findmarkers','mnn_opt++findmarkers','scMerge++findmarkers','raw++findmarkers','zinbwave++findmarkers', "Seurat++findmarkers", 'raw++pseudo_DESeq2','raw++pseudo_edgeR','raw++pseudo_voom.limma','raw++pseudo_limma_trend','raw++MAST','raw++MAST_Cov','raw++DESeq2','raw++DESeq2_Cov','zinbwave++DESeq2_pseudo','zinbwave++DESeq2_pseudo_Cov','raw+ind.DESeq2.pval+Fisher','raw+ind.DESeq2.pval+wFisher','raw++edgeR','raw++edgeR_Cov','raw++edgeR_DetRate','raw++edgeR_DetRate_Cov','zinbwave++edgeR','zinbwave++edgeR_Cov','raw+ind.edgeR.pval+Fisher','raw+ind.edgeR.pval+wFisher','raw++limma','raw++limma_Cov','raw++limma_trend','raw++limma_trend_Cov','combat++limma_trend','mnn_opt++limma_trend','scMerge++limma_trend','LogNormalize+ind.limma_trend+Fisher','LogNormalize+ind.limma_trend+wFisher','raw+ind.DESeq2.ES+FEM','voom+ind.ES+FEM','LogNormalize+ind.ES+FEM','raw+ind.DESeq2.ES+REM','voom+ind.ES+REM','LogNormalize+ind.ES+REM','voom+ind.modt+Fisher','voom+ind.modt+wFisher', 'raw++LUAD-GSE43458.microarray.limma.cov', 'raw++LUAD-GSE10072.microarray.limma.cov', 'raw++LUAD-GSE31210.microarray.limma.cov')
  
  meth_names<-c('base_line','raw..TCGA_DESeq2_Cov','raw..TCGA_edgeR_Cov','raw..TCGA_voom.limma_Cov','raw..TCGA_limma_trend_Cov','combat..findmarkers','limma_bec..findmarkers','mnn_opt..findmarkers','scMerge..findmarkers','raw..findmarkers','zinbwave..findmarkers',"Seurat..findmarkers",'raw..pseudo_DESeq2','raw..pseudo_edgeR','raw..pseudo_voom.limma','raw..pseudo_limma_trend','raw..MAST','raw..MAST_Cov','raw..DESeq2','raw..DESeq2_Cov','zinbwave..DESeq2_pseudo','zinbwave..DESeq2_pseudo_Cov','raw.ind.DESeq2.pval.Fisher','raw.ind.DESeq2.pval.wFisher','raw..edgeR','raw..edgeR_Cov','raw..edgeR_DetRate','raw..edgeR_DetRate_Cov','zinbwave..edgeR','zinbwave..edgeR_Cov','raw.ind.edgeR.pval.Fisher','raw.ind.edgeR.pval.wFisher','raw..limma','raw..limma_Cov','raw..limma_trend','raw..limma_trend_Cov','combat..limma_trend','mnn_opt..limma_trend','scMerge..limma_trend','LogNormalize.ind.limma_trend.Fisher','LogNormalize.ind.limma_trend.wFisher','raw.ind.DESeq2.ES.FEM','voom.ind.ES.FEM','LogNormalize.ind.ES.FEM','raw.ind.DESeq2.ES.REM','voom.ind.ES.REM','LogNormalize.ind.ES.REM','voom.ind.modt.Fisher','voom.ind.modt.wFisher','raw..LUAD-GSE43458.microarray.limma.cov', 'raw..LUAD-GSE10072.microarray.limma.cov', 'raw..LUAD-GSE31210.microarray.limma.cov')
  df<-data.frame(methods=Methods_bk, raw=meth_raw, raw2=meth_names)
  df_where<-sapply(df, FUN=function(m){str_detect(m,pattern=paste0(x,'$'))})

  if(x%in%c(df$raw,df$raw2, df$methods)){
    return(df$methods[which(x==df,arr.ind = T)[1]])
  }else{
    if(any(df_where)){
      return(df$methods[which(df_where, arr.ind=T)[1,1]])
    }else{
      return(x)
    }
  }
}
get_res_matrix<-function(ct,meth.comb,base_dir,filter.rate=0.01,filt_dir,fc.direction='high'){
  message(meth.comb)
  ind.analysis=F
  first.meth=str_split(meth.comb,pattern='[+]')[[1]][1]
  second.meth=str_split(meth.comb,pattern='[+]')[[1]][2]
  third.meth=str_split(meth.comb,pattern='[+]')[[1]][3]
  if(str_detect(second.meth,pattern='ind')){
    ind.analysis=T
  }
  if(str_detect(pattern='TCGA',third.meth)){
    message(paste0(base_dir,'data/analyzed/TCGA-LUAD.rnaseq.processed.htseq_count.',gsub(pattern='TCGA_',replacement = '',third.meth),'.txt'))
    third_processed<-read.table(paste0(base_dir,'data/analyzed/TCGA-LUAD.rnaseq.processed.htseq_count.',gsub(pattern='TCGA_',replacement = '',third.meth),'.txt'))
  }else if(str_detect(pattern='LUAD-GSE',third.meth)){
    message(paste0(base_dir,'data/analyzed/',third.meth,'.txt'))
    third_processed<-read.table(paste0(base_dir,'data/analyzed/',third.meth,'.txt'))
  }else{
    message(paste0(base_dir,'data/',filt_dir,'/third_processed/',ct,ifelse(is.null(filter.rate),'',paste0('.',filter.rate,'mincell_filtered')),'_',first.meth,'|',second.meth,'|',third.meth,'_third_processed.RData'))
    load(paste0(base_dir,'data/',filt_dir,'/third_processed/',ct,ifelse(is.null(filter.rate),'',paste0('.',filter.rate,'mincell_filtered')),'_',first.meth,'|',second.meth,'|',third.meth,'_third_processed.RData'))
  }
  
  if(third.meth %in% c('wFisher')){
    wfisher.dir<-ifelse(third_processed$two_tailed$direction[,1]=='+',1,-1)
    sign.check<-sign(third_processed$logFC)[names(wfisher.dir)]*wfisher.dir
    log2FC=third_processed$logFC
    minfc<-log2FC[sign.check>=0]%>%abs%>%min
    log2FC[sign.check<0]=sign(log2FC[sign.check<0])*(-1)*ifelse(minfc<1,minfc^2,minfc^(-2))
    res<-data.frame(log2FC=log2FC,pval=third_processed[['two_tailed']]$pval[,1],adj.pval=third_processed[['two_tailed']]$FDR[,1])
    rownames(res)<-rownames(third_processed[['two_tailed']]$FDR)
  }else if(third.meth %in% c('Fisher')){
    res<-data.frame(log2FC=third_processed$logFC[,1],pval=third_processed[['two_tailed']]$pval[,1],adj.pval=third_processed[['two_tailed']]$FDR)
    
  }else if(third.meth %in% c('REM','FEM')){
    # pval.twotailed<-two_tailed_for_ef.meta(third_processed)
    # FDR.twotailed=p.adjust(pval.twotailed,method='BH')
    
    # res<-data.frame(log2FC=third_processed[[fc.direction]]$mu.hat,pval=third_processed[[fc.direction]]$pval, adj.pval=third_processed[[fc.direction]]$FDR)
    # res<-data.frame(log2FC=third_processed[[fc.direction]]$mu.hat,pval=pval.twotailed, adj.pval=FDR.twotailed)
    # rownames(res)<-names(third_processed$high$FDR)    
    third_processed$high$pval
    res<-data.frame(log2FC=third_processed[['abs']]$mu.hat,pval=third_processed[['abs']]$pval, adj.pval=third_processed[['abs']]$FDR)
    rownames(res)<-names(third_processed[['abs']]$FDR)
  }else{
    if(third.meth%in%c('findmarkers','MAST','MAST_Cov')){
      res<-data.frame(log2FC=third_processed$avg_log2FC,pval=third_processed$p_val, adj.pval=third_processed$p_val_adj)
      rownames(res)<-rownames(third_processed)
    }else{
      if(anyNA(third_processed$adjpvalue)){
        message('NA adjpvalues')
        res<-data.frame(log2FC=third_processed$logFC,pval=third_processed$pvalue,adj.pval=p.adjust(third_processed$pvalue,method='BH'))
      }else{
        res<-data.frame(log2FC=third_processed$logFC,pval=third_processed$pvalue,adj.pval=third_processed$adjpvalue)
      }
      rownames(res)<-rownames(third_processed)
    }
  }
  message(nrow(res))
  return(res)
}

get_res_matrix.ind<-function(ct,meth.comb,base_dir,filter.rate=0.01,filt_dir,return.type='logfc'){
  ind.analysis=F
  first.meth=str_split(meth.comb,pattern='[+]')[[1]][1]
  second.meth=str_split(meth.comb,pattern='[+]')[[1]][2]
  third.meth=str_split(meth.comb,pattern='[+]')[[1]][3]
  message(paste0(base_dir,'data/',filt_dir,'/second_processed/',ct,ifelse(is.null(filter.rate),'',paste0('.',filter.rate,'mincell_filtered')),'_',first.meth,'|',second.meth,'_second_processed.RData'))
  load(paste0(base_dir,'data/',filt_dir,'/second_processed/',ct,ifelse(is.null(filter.rate),'',paste0('.',filter.rate,'mincell_filtered')),'_',first.meth,'|',second.meth,'_second_processed.RData'))
  
  if(second.meth %in% c('ind.edgeR.pval',"ind.DESeq2.pval")){
    if(return.type=='logfc'){
      return(second_processed$logFC%<>%as.data.frame())
    }else if(return.type%in%c('pval','qval')){
      two.tailed<-second_processed[['high']]<second_processed[['low']]
      
      for(i in 1:nrow(two.tailed)){
        for(j in 1:ncol(two.tailed)){
          two.tailed[i,j]=ifelse(two.tailed[i,j],second_processed[['high']][i,j],second_processed[['low']][i,j])*2
        }
      }
      if(return.type=='qval'){
        two.tailed<-sapply(as.data.frame(two.tailed),FUN=function(x){p.adjust(x,method='BH')})
      }
      return(two.tailed)
    }else{
      message('set return type')
      break()
    }
  }else{
    message('set other cases')
    break()
  }
}



get_plot_df_saved<-function(saved_dir,ct,meth.comb.tot,ref_genes,ref_weight,base_dir,filter.rate=0.01,filt_dir,adj.p_thr,sub.ct,sort.meth='logFC',up.genes='up'){
  load(paste0(saved_dir,up.genes,' | ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),'gene threshold : expressed in  ',filter.rate*100,'percent cells | adj.p_threshold : ',adj.p_thr,'|',sort.meth, ' sorted.ks.RData'))
  TCGA.exist=any(str_detect(meth.comb.tot,pattern='TCGA'))
  if(adj.p_thr=='none'){
    adj.p_thr=2.2
  }
  res_plots<-data.frame(x=numeric(),y=numeric(),meth.comb=character())
  
  meth.comb.ind<-meth.comb.tot[which(sapply(meth.comb.tot,FUN=function(x)str_split(x,pattern='[+]')[[1]][3])=='ind')]
  meth.comb.tot<-meth.comb.tot[which(sapply(meth.comb.tot,FUN=function(x)str_split(x,pattern='[+]')[[1]][3])!='ind')]
  message(paste0(meth.comb.ind, ' is ind'))
  meth.comb.ind.patients<-c()
  for(inds in meth.comb.ind){
    # print(names(res_list)[str_detect(names(res_list),pattern=gsub(pattern='[+]',replacement='.',inds))])
    # meth.comb.ind.patients%<>%c(names(res_list)[str_detect(names(meth.comb.tot),pattern=gsub(pattern='[+]',replacement='.',inds))])
    meth.comb.ind.patients%<>%c(names(res_list)[str_detect(names(res_list),pattern=gsub(pattern='[+]',replacement='.',inds))])
  }
  
  
  meth.comb.tot%<>%c(meth.comb.ind.patients)
  
  for(meth.comb in meth.comb.tot){
    message(meth.comb)
    second.meth=str_split(meth.comb,pattern='[+]')[[1]][2]
    third.meth=str_split(meth.comb,pattern='[+]')[[1]][3]
    if(adj.p_thr!='none'){adj.p_thr%<>%as.numeric()}
    
    ranked_genes<-res_list[[meth.comb]]$genes
    if(sort.meth=='logFC'){
      message(meth.comb)
      sort2<-abs(res_list[[meth.comb]][['log2FC']])
    }else{
      if(sort.meth=='adj.pvalue'){
        sort<-res_list[[meth.comb]][["adj.pval"]]
        
      }else if(sort.meth=='pvalue'){
        sort<-res_list[[meth.comb]][["pval"]]
      }
      if(length(sort)==0){
        sort=sort2=NULL
      }else{
        abs.logfc<-abs(res_list[[meth.comb]][['log2FC']])
        sort<-sapply(sort,FUN=function(x)which(x==(sort%>%unique()%>%sort(decreasing = T))))+1
        sort2=sort+abs.logfc/ifelse(max(abs.logfc)<1,1,max(abs.logfc))
      }
      
    } 
    if(is.null(sort)){
      si=0
      temp_x<-c(0,0)
      temp_y<-c(0,0)
    }else{
      si<-sapply(c(1:length(sort2)), FUN = function(x){i=x
      if(sum(sort2==sort2[i])>1){i=(which(sort2==sort2[i])%>%max())}
      return(i)})
      si%<>%unique()
      temp_x<-c(0,si)
      temp_y<-c(0,sapply(si,FUN=function(x){ref_weight[ref_genes %in% intersect(ranked_genes[1:x],ref_genes)]%>%sum()}))
    }
    temp_plots<-data.frame(x=temp_x,y=temp_y,methods=rep(meth.comb,length(temp_x)))
    res_plots%<>%rbind(temp_plots)
    
  }
  ed_x<-length(res_list$genes)
  ed_y<-ref_weight[ref_genes %in% intersect(res_list$genes,ref_genes)]%>%sum()
  
  res_plots%<>%rbind(data.frame(x=c(0,ed_x),y=c(0,ed_y),methods=rep(paste0(ed_y, '/', ed_x),2)))
  return(res_plots)
  
}
write.ks.adjp_thr<-function(saved_dir,ct,base_dir,filter.rate=0.01,filt_dir,adj.p_thr,sub.ct,sort.meth='logFC',up.genes='both'){
  load(paste0(saved_dir,up.genes,' | ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),'gene threshold : expressed in  ',filter.rate*100,'percent cells | adj.p_threshold : none|',sort.meth, ' sorted.ks.RData'))
  meth.comb.tot<-setdiff(names(res_list),'genes')
  TCGA.exist=any(str_detect(meth.comb.tot,pattern='TCGA'))
  if(adj.p_thr=='none'){
    adj.p_thr=2.2
  }
  meth.names<-setdiff(names(res_list),'genes')
  
  for(meth.name in meth.names){
    
    temp.res<-res_list[[meth.name]]
    keep.which<-which(temp.res$adj.pval<adj.p_thr)
    temp.res$genes<-temp.res$genes[keep.which]
    temp.res$pval<-temp.res$pval[keep.which]
    temp.res$adj.pval<-temp.res$adj.pval[keep.which]
    temp.res$log2FC<-temp.res$log2FC[keep.which]
    res_list[[meth.name]]=temp.res
  }
  save(res_list,file=paste0(saved_dir,up.genes,' | ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),'gene threshold : expressed in  ',filter.rate*100,'percent cells | adj.p_threshold : ',adj.p_thr,'|',sort.meth, ' sorted.ks.RData'))

}
write_ks_df<-function(ct,meth.comb.tot,ref_df_refined,base_dir,filter.rate=0.01,filt_dir,adj.p_thr,ref,sub.ct,fig_dir,sort.meth='logFC',up.genes='up',use.sort=F,seed='none'){
  adj.p_thr_origin<-adj.p_thr
  if(adj.p_thr=='none'){
    adj.p_thr=2.2
  }
  TCGA.exist=any(str_detect(meth.comb.tot,pattern='TCGA'))
  
  genes<-get_res_matrix(ct,meth.comb='raw++findmarkers',base_dir,filter.rate=filter.rate,filt_dir)%>%rownames()
  if(TCGA.exist){
    genes<-intersect(get_res_matrix(ct,meth.comb='raw++TCGA_DESeq2',base_dir=base_dir)%>%rownames(),genes)
  }
  marray.exist=any(str_detect(meth.comb.tot,pattern='LUAD-GSE'))
  if(marray.exist){
    for(m in meth.comb.tot[str_detect(meth.comb.tot,pattern='LUAD-GSE')]){
      genes<-intersect(get_res_matrix(ct,meth.comb=m,base_dir=base_dir)%>%rownames(),genes)
    }
  }
  # genes<-get_res_matrix(ct,meth.comb='raw++findmarkers',base_dir,filter.rate=filter.rate,filt_dir)%>%rownames()
  res_list<-list()
  res_list[['genes']]=genes
  
  meth.comb.ind<-meth.comb.tot[which(sapply(meth.comb.tot,FUN=function(x)str_split(x,pattern='[+]')[[1]][3])=='ind')]
  meth.comb.tot<-meth.comb.tot[which(sapply(meth.comb.tot,FUN=function(x)str_split(x,pattern='[+]')[[1]][3])!='ind')]
  message(paste0(meth.comb.ind, ' is ind'))
  
  for(meth.comb in meth.comb.tot){
    third.meth=str_split(meth.comb,pattern='[+]')[[1]][3]
    res_matrix<-get_res_matrix(ct,meth.comb,base_dir,filter.rate=filter.rate,filt_dir)
    plot_df<-res_matrix%>%dplyr::filter(adj.pval<adj.p_thr)
    plot_df<-plot_df[which(rownames(plot_df)%in%genes),]
    if(up.genes=='up'){
      plot_df%<>%dplyr::filter(log2FC>0)
    }else if(up.genes=='down'){
      plot_df%<>%dplyr::filter(log2FC<=0)
    }else if(up.genes=='both'){
    }
    if(sort.meth=='logFC'){
      plot_df%<>%dplyr::arrange(desc(abs(log2FC)))
      plot_df$sort=abs(plot_df$log2FC)
    }else if(sort.meth=='adj.pvalue'){
      plot_df%<>%dplyr::arrange(adj.pval)
      plot_df$sort=plot_df$adj.pval
      abs.logfc<-abs(plot_df$log2FC)
      plot_df$sort<-sapply(plot_df$sort,FUN=function(x)which(x==(plot_df$sort%>%unique()%>%sort(decreasing = T))))+1
      
      plot_df$sort2=plot_df$sort+abs.logfc/ifelse(max(abs.logfc)<1,1,max(abs.logfc))
      
      plot_df%<>%dplyr::arrange(desc(sort2))
      
    }else if(sort.meth=='pvalue'){
      plot_df$sort=plot_df$pval
      abs.logfc<-abs(plot_df$log2FC)
      plot_df$sort<-sapply(plot_df$sort,FUN=function(x)which(x==(plot_df$sort%>%unique()%>%sort(decreasing = T))))+1
      
      plot_df$sort2=plot_df$sort+abs.logfc/ifelse(max(abs.logfc)<1,1,max(abs.logfc))
      
      plot_df%<>%dplyr::arrange(desc(sort2))
    }
    ranked_genes<-rownames(plot_df)
    res_list[[meth.comb]][['genes']]=ranked_genes
    res_list[[meth.comb]][['pval']]=plot_df$pval
    res_list[[meth.comb]][['adj.pval']]=plot_df$adj.pval
    res_list[[meth.comb]][['log2FC']]=plot_df$log2FC
  }
  
  for(meth.comb in meth.comb.ind){
    second.meth=str_split(meth.comb,pattern='[+]')[[1]][2]
    third.meth=str_split(meth.comb,pattern='[+]')[[1]][3]
    res_matrix<-get_res_matrix.ind(ct,meth.comb,base_dir,filter.rate=filter.rate,filt_dir,return.type = 'qval')
    res_matrix_pval<-get_res_matrix.ind(ct,meth.comb,base_dir,filter.rate=filter.rate,filt_dir,return.type = 'pval')
    fc_matrix<-get_res_matrix.ind(ct,meth.comb,base_dir,filter.rate=filter.rate,filt_dir,return.type = 'logfc')
    colnames(fc_matrix)<-paste0(colnames(fc_matrix),'_logfc')
    colnames(res_matrix_pval)<-paste0(colnames(res_matrix_pval),'_pval')
    # all(rownames(res_matrix)==rownames(fc_matrix))
    tot_matrix<-as.data.frame(cbind(res_matrix,cbind(res_matrix_pval,fc_matrix)))
    
    for(ind.patient in colnames(res_matrix)){
      
      plot_df<-tot_matrix%>%dplyr::select(all_of(c(ind.patient,paste0(ind.patient,'_pval'),paste0(ind.patient,'_logfc'))))
      plot_df<-plot_df[which(rownames(plot_df)%in%genes),]
      if(up.genes=='up'){
        plot_df%<>%dplyr::filter(get(paste0(ind.patient,'_logfc'))>0)
      }else if(up.genes=='down'){
        plot_df%<>%dplyr::filter(get(paste0(ind.patient,'_logfc'))<=0)
      }else if(up.genes=='both'){
      }
      
      if(sort.meth=='logFC'){
        plot_df%<>%dplyr::arrange(desc(abs(get(paste0(ind.patient,'_logfc')))))
        plot_df$sort=abs(plot_df[[paste0(ind.patient,'_logfc')]])
        plot_df$sort2=plot_df$sort
      }else if(sort.meth=='adj.pvalue'){
        plot_df%<>%dplyr::arrange(get(ind.patient))
        plot_df$sort=plot_df[[ind.patient]]
        abs.logfc<-plot_df[[paste0(ind.patient,'_logfc')]]%>%abs()
        plot_df$sort<-sapply(plot_df$sort,FUN=function(x)which(x==(plot_df$sort%>%unique()%>%sort(decreasing = T))))+1

        plot_df$sort2=plot_df$sort+abs.logfc/ifelse(max(abs.logfc)<1,1,max(abs.logfc))
        
        
        # plot_df$sort2=((1-plot_df[[ind.patient]])*10^8)+abs(plot_df[[paste0(ind.patient,'_logfc')]])
        plot_df%<>%dplyr::arrange(desc(sort2))
      }else if(sort.meth=='pvalue'){
        plot_df%<>%dplyr::arrange(get(paste0(ind.patient,'_pval')))
        plot_df$sort=plot_df[[paste0(ind.patient,'_pval')]]
        abs.logfc<-plot_df[[paste0(ind.patient,'_logfc')]]%>%abs()
        plot_df$sort<-sapply(plot_df$sort,FUN=function(x)which(x==(plot_df$sort%>%unique()%>%sort(decreasing = T))))+1
        
        plot_df$sort2=plot_df$sort+abs.logfc/ifelse(max(abs.logfc)<1,1,max(abs.logfc))
        
        plot_df%<>%dplyr::arrange(desc(sort2))
      }
      ranked_genes<-rownames(plot_df)
      res_list[[paste0(meth.comb,'_',ind.patient)]][['genes']]=ranked_genes
      res_list[[paste0(meth.comb,'_',ind.patient)]][['adj.pval']]=plot_df[[ind.patient]]
      res_list[[paste0(meth.comb,'_',ind.patient)]][['pval']]=plot_df[[paste0(ind.patient,'_pval')]]
      res_list[[paste0(meth.comb,'_',ind.patient)]][['log2FC']]=plot_df[[paste0(ind.patient,'_logfc')]]
    }
  }
  print(paste0(fig_dir,up.genes,' | ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),'gene threshold : expressed in  ',filter.rate*100,'percent cells | adj.p_threshold : ',adj.p_thr_origin, '|',sort.meth, ' sorted',ifelse(seed=='none','',paste0('_seed',seed)),'.ks.RData'))
  save(res_list, file=paste0(fig_dir,up.genes,' | ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),'gene threshold : expressed in  ',filter.rate*100,'percent cells | adj.p_threshold : ',adj.p_thr_origin, '|',sort.meth, ' sorted',ifelse(seed=='none','',paste0('_seed',seed)),'.ks.RData'))
  
}

do.mapping<-function(geo_acc,eset,fset){
  if(!is.null(fset[['Gene Symbol']])){
    rownames(eset)<-fset$`Gene Symbol`[match(rownames(eset),fset$ID)]
  }
  eset<-eset[str_detect(rownames(eset),pattern='[///]', negate = T)%>%which(),]
  mart <- useMart('ENSEMBL_MART_ENSEMBL',host = 'asia.ensembl.org')
  # mart <- useMart('ENSEMBL_MART_ENSEMBL')
  mart <- useDataset("hsapiens_gene_ensembl", mart)
  if(geo_acc=='GSE32863'){
    platform='illumina_humanwg_6_v3'
    filt=attr='external_gene_name'
    
  }else if(geo_acc=='GSE10072'){
    platform='affy_hg_u133a'
    filt=attr='external_gene_name'
  }else if(geo_acc=='GSE43458'){
    platform='affy_hugene_1_0_st_v1'
    filt=attr='affy_hugene_1_0_st_v1'
  }else{
    platform='affy_hg_u133_plus_2'
    filt=attr='external_gene_name'
  }
  
  mrna_attributes <- getBM(mart = mart,
                           # attributes = c(platform,
                           attributes = c(platform,
                                          'ensembl_gene_id',
                                          'gene_biotype',
                                          'external_gene_name'),
                           filter = filt,
                           values = rownames(eset),
                           uniqueRows = TRUE)
  mrna_attributes%<>%dplyr::filter(gene_biotype=='protein_coding')
  eset<-eset[rownames(eset)%in%mrna_attributes[[attr]],]
  rownames(eset)<-mrna_attributes$external_gene_name[match(rownames(eset),mrna_attributes[[attr]])]
  
  res_dt<-eset%>%as.data.frame()%>%dplyr::group_by(GeneID=rownames(eset))%>%summarise_all(mean)%>%as.data.frame()
  rownames(res_dt)<-res_dt$GeneID
  res_dt<-res_dt[,2:ncol(res_dt)]
  res_dt<-res_dt[setdiff(unique(rownames(res_dt))%>%sort(),c("")),]
  return(res_dt)
}

draw_plot_bar<-function(fig_dir,meth.comb.tot,ref,up.genes, filter.rate, ct,sub.ct,seed='none',ngene){
  res.table<-read.table(paste0(fig_dir,ref,'.',up.genes,' | ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),'gene threshold : expressed in  ',filter.rate*100,'percent cells | N_DEG',ifelse(seed=='none','',paste0('_seed',seed)),'.txt'),sep = '\t',check.names=FALSE)
  `%ni%`<-Negate(`%in%`)
  library(reshape2)
  
  raw_ref<-c('base_line','combat++findmarkers','limma_bec++findmarkers','mnn_opt++findmarkers','scMerge++findmarkers','raw++findmarkers','zinbwave++findmarkers', "Seurat++findmarkers",'raw++pseudo_DESeq2','raw++pseudo_edgeR','raw++pseudo_voom.limma','raw++pseudo_limma_trend','raw++MAST','raw++MAST_Cov','raw++DESeq2','raw++DESeq2_Cov','zinbwave++DESeq2_pseudo','zinbwave++DESeq2_pseudo_Cov','raw+ind.DESeq2.pval+Fisher','raw+ind.DESeq2.pval+wFisher','raw++edgeR','raw++edgeR_Cov','raw++edgeR_DetRate','raw++edgeR_DetRate_Cov','zinbwave++edgeR','zinbwave++edgeR_Cov','raw+ind.edgeR.pval+Fisher','raw+ind.edgeR.pval+wFisher','raw++limma','raw++limma_Cov','raw++limma_trend','raw++limma_trend_Cov','combat++limma_trend','mnn_opt++limma_trend','scMerge++limma_trend','LogNormalize+ind.limma_trend+Fisher','LogNormalize+ind.limma_trend+wFisher','raw+ind.DESeq2.ES+FEM','voom+ind.ES+FEM','LogNormalize+ind.ES+FEM','raw+ind.DESeq2.ES+REM','voom+ind.ES+REM','LogNormalize+ind.ES+REM','voom+ind.modt+Fisher','voom+ind.modt+wFisher')
  
  res.table<-res.table[,which(colnames(res.table)%in%raw_ref)]
  col_meths.use=intersect(colnames(res.table), meth.comb.tot)
  res.table<-res.table[,which(colnames(res.table)%in%col_meths.use)]
  colnames(res.table)<-colnames(res.table)%>%sapply(FUN=function(x)return(select_bk_df(x)))%>%as.vector()
  col_meths.use<-col_meths.use%>%sapply(FUN=function(x)return(select_bk_df(x)))%>%as.vector()
  
  default.order=c('base_line','TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend','GSE43458.limma_Cov','GSE31210.limma_Cov','Pseudobulk_DESeq2','Pseudobulk_edgeR','Pseudobulk_limma','Pseudobulk_limma_trend','Combat_Wilcox','limma_BEC_Wilcox','MNNCorrect_Wilcox','scMerge_Wilcox','Seurat_Wilcox','ZINB-WaVE_Wilcox','Raw_Wilcox','MAST','MAST_Cov','DESeq2','DESeq2_Cov','ZINB-WaVE_DESeq2','ZINB-WaVE_DESeq2_Cov','DESeq2_Fisher','DESeq2_wFisher','edgeR_DetRate','edgeR_DetRate_Cov','edgeR','edgeR_Cov','ZINB-WaVE_edgeR','ZINB-WaVE_edgeR_Cov','edgeR_Fisher','edgeR_wFisher','limma','limma_Cov','limmatrend','limmatrend_Cov','Combat_limmatrend','MNNCorrect_limmatrend','scMerge_limmatrend','LogNorm+limmatrend_Fisher','LogNorm+limmatrend_wFisher','DESeq2_FEM','voom_FEM','LogNorm_FEM','DESeq2_REM','voom_REM','LogNorm_REM','voom+modt_Fisher','voom+modt_wFisher')
  
  
  
  library(RColorBrewer)
  library(ggplot2)
  ggplotColours <- function(n = 6, h = c(0, 360) + 15){
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
  }
  `%ni%`<-Negate(`%in%`)
  plot_col<-ggplotColours(n=64)[1:(length(default.order)-11)]
  names(plot_col)<-setdiff(default.order,c('TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend','GSE43458.limma_Cov','GSE31210.limma_Cov','Pseudobulk_DESeq2','Pseudobulk_edgeR','Pseudobulk_limma','Pseudobulk_limma_trend'))
  library(colortools)
  plot_col[['GSE43458.limma_Cov']]=plot_col[['GSE31210.limma_Cov']]=plot_col[['TCGA_limma_trend']]="green"
  plot_col[['TCGA_limma']]=sequential(plot_col[['TCGA_limma_trend']])[[15]]
  plot_col[['TCGA_edgeR']]=sequential(plot_col[['TCGA_limma_trend']])[[10]]
  plot_col[['TCGA_DESeq2']]=sequential(plot_col[['TCGA_limma_trend']])[[5]]
  plot_col[['Pseudobulk_limma_trend']]='blue'
  plot_col[['Pseudobulk_limma']]=sequential(plot_col[['Pseudobulk_limma_trend']])[[15]]
  plot_col[['Pseudobulk_edgeR']]=sequential(plot_col[['Pseudobulk_limma_trend']])[[10]]
  plot_col[['Pseudobulk_DESeq2']]=sequential(plot_col[['Pseudobulk_limma_trend']])[[5]]
  plot_col<-plot_col[order(match(names(plot_col),default.order))]
  
  
  axis_order<-intersect(default.order,col_meths.use%>%sapply(FUN=function(x)return(select_bk_df(x)))%>%as.vector())
  
  meths=cutoff=cutoff.metric=values=cols=c()
  for(thr in c(0.01,0.05)){
    for(metric in c('pval','adj.pval')){
      values%<>%c(res.table[paste0(thr,'_',metric),col_meths.use]%>%as.numeric())
      meths%<>%c(col_meths.use)
      cutoff%<>%c(rep(thr,length(col_meths.use)))
      cutoff.metric%<>%c(rep(metric,length(col_meths.use)))
      cols%<>%c(plot_col[col_meths.use])
    }
  }
  res = data.frame(Methods = meths, Metric=cutoff.metric, Cutoff=cutoff, FP=values, Color = cols)
  res$Methods%<>%factor()
  res$Methods<-factor(res$Methods,levels=axis_order)
  res$Metric%<>%factor()
  res$Metric<-factor(res$Metric,levels=rev(levels(res$Metric)))
  res$Cutoff%<>%factor()
  res$Cutoff<-factor(res$Cutoff,levels=rev(levels(res$Cutoff)))
  res$Color%<>%as.factor()
  res2 = melt(res, measure.vars=c("FP"))
  for(thr in c(0.01,0.05)){
    for(metric in c('pval','adj.pval')){
      pd = position_dodge(width=0.0)
      gbase = ggplot(res2, aes(y=value, x=Methods, fill=Methods))+geom_bar(position=pd, stat = 'identity')+
        # facet_grid(Cutoff~Metric, scales='free')+
        facet_grid(Metric~Cutoff, scales='free')+
        scale_x_discrete(limits = axis_order)+
        theme(axis.text.x=element_text(angle=45, hjust=1), legend.direction = 'horizontal', legend.position = 'bottom')+
        scale_colour_manual(values = plot_col[which(names(plot_col)%ni%c('base_line'))])+
        guides(colour=guide_legend(ncol=2))
      # scale_colour_manual(name = "Methods",
      #                     labels = ts[order(ts)],
      #                     values = tc[order(ts)])
      
      gline = gbase
      print(gline)
      tt = paste(ref,' / ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),' / min.cell = ',filter.rate*100,' percent ',sep="")
      
      print(gline+aes(x=Methods)+labs(x='Methods', y='FP')+ggtitle(tt))
      figurename=paste0(ref,'.',up.genes,' | ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),'gene threshold : expressed in  ',filter.rate*100,'percent cells | N_DEG.pdf')
      ggsave(file=paste(fig_dir,"/",figurename,sep=""),width = 25,height = 8)
      dev.off()
    }
  }
}


draw_plot_box<-function(fig_dir,meth.comb.tot,ref,up.genes, filter.rate, ct,sub.ct,seed='none',ngene){
  thr.use=0.05
  `%ni%`<-Negate(`%in%`)
  
  raw_ref<-c('base_line','combat++findmarkers','limma_bec++findmarkers','mnn_opt++findmarkers','scMerge++findmarkers',"Seurat++findmarkers",'raw++findmarkers','zinbwave++findmarkers', 'raw++pseudo_DESeq2','raw++pseudo_edgeR','raw++pseudo_voom.limma','raw++pseudo_limma_trend','raw++MAST','raw++MAST_Cov','raw++DESeq2','raw++DESeq2_Cov','zinbwave++DESeq2_pseudo','zinbwave++DESeq2_pseudo_Cov','raw+ind.DESeq2.pval+Fisher','raw+ind.DESeq2.pval+wFisher','raw++edgeR','raw++edgeR_Cov','raw++edgeR_DetRate','raw++edgeR_DetRate_Cov','zinbwave++edgeR','zinbwave++edgeR_Cov','raw+ind.edgeR.pval+Fisher','raw+ind.edgeR.pval+wFisher','raw++limma','raw++limma_Cov','raw++limma_trend','raw++limma_trend_Cov','combat++limma_trend','mnn_opt++limma_trend','scMerge++limma_trend','LogNormalize+ind.limma_trend+Fisher','LogNormalize+ind.limma_trend+wFisher','raw+ind.DESeq2.ES+FEM','voom+ind.ES+FEM','LogNormalize+ind.ES+FEM','raw+ind.DESeq2.ES+REM','voom+ind.ES+REM','LogNormalize+ind.ES+REM','voom+ind.modt+Fisher','voom+ind.modt+wFisher')
  

  default.order=c('base_line','TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend','GSE43458.limma_Cov','GSE31210.limma_Cov','Pseudobulk_DESeq2','Pseudobulk_edgeR','Pseudobulk_limma','Pseudobulk_limma_trend','Combat_Wilcox','limma_BEC_Wilcox','MNNCorrect_Wilcox','scMerge_Wilcox','Seurat_Wilcox','ZW_Wilcox','Raw_Wilcox','MAST','MAST_Cov','DESeq2','DESeq2_Cov','ZW_DESeq2','ZW_DESeq2_Cov','DESeq2_Fisher','DESeq2_wFisher','edgeR_DetRate','edgeR_DetRate_Cov','edgeR','edgeR_Cov','ZW_edgeR','ZW_edgeR_Cov','edgeR_Fisher','edgeR_wFisher','limma','limma_Cov','limmatrend','limmatrend_Cov','Combat_limmatrend','MNNCorrect_limmatrend','scMerge_limmatrend','LogNorm+limmatrend_Fisher','LogNorm+limmatrend_wFisher','DESeq2_FEM','voom_FEM','LogNorm_FEM','DESeq2_REM','voom_REM','LogNorm_REM','voom+modt_Fisher','voom+modt_wFisher')
  
  seeds=meths=cutoff=cutoff.metric=values=cols=c()
  for(s in seed){
    print(paste0(fig_dir,ref,'.',up.genes,' | ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),'gene threshold : expressed in  ',filter.rate*100,'percent cells | N_DEG',ifelse(s=='none','',paste0('_seed',s)),'.txt'))
    # res.table<-read.table(paste0(fig_dir,ref,'.',up.genes,' | ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),'gene threshold : expressed in  ',filter.rate*100,'percent cells | N_DEG',ifelse(s=='none','',paste0('_seed',s)),'.txt'),sep = '\t',check.names=FALSE)
    res.table<-read.table(paste0(fig_dir,ref,'.',up.genes,' | ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),'gene threshold : expressed in  ',filter.rate*100,'percent cells | N_DEG',ifelse(s=='none','',paste0('_seed',s)),'.txt'),sep = '\t',check.names=FALSE)
    ngene=res.table['Total genes',1]
    # meth.comb.tot<-colnames(res.table)
    res.table<-res.table[,which(colnames(res.table)%in%raw_ref)]
    col_meths.use=intersect(colnames(res.table), meth.comb.tot)
    res.table<-res.table[,which(colnames(res.table)%in%col_meths.use)]
    colnames(res.table)<-colnames(res.table)%>%sapply(FUN=function(x)return(select_bk_df(x)))%>%as.vector()
    col_meths.use<-col_meths.use%>%sapply(FUN=function(x)return(select_bk_df(x)))%>%as.vector()
    
    
    for(thr in c(thr.use)){
      for(metric in c('pval','adj.pval')){
        values%<>%c(res.table[paste0(thr,'_',metric),col_meths.use]%>%as.numeric())
        meths%<>%c(col_meths.use)
        cutoff%<>%c(rep('',length(col_meths.use)))
        if(metric=='pval'){
          cutoff.metric%<>%c(rep(paste0('p-values < ',thr),length(col_meths.use)))
        }else if(metric=='adj.pval'){
          cutoff.metric%<>%c(rep(paste0('q-values < ',thr),length(col_meths.use)))
        }else{
          stop('metric should be pval or adj.pval')
        }
        
        # cols%<>%c(df_col$color[match(col_meths.use,df_col$name)])
      }
    }
    seeds<-c(seeds,rep(s,(length(values)-length(seeds))))
  }
  library(RColorBrewer)
  library(ggplot2)
  ggplotColours <- function(n = 6, h = c(0, 360) + 15){
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
  }
  `%ni%`<-Negate(`%in%`)
  plot_col<-ggplotColours(n=64)[1:(length(default.order)-11)]
  names(plot_col)<-setdiff(default.order,c('TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend','GSE43458.limma_Cov','GSE31210.limma_Cov','Pseudobulk_DESeq2','Pseudobulk_edgeR','Pseudobulk_limma','Pseudobulk_limma_trend'))
  library(colortools)
  plot_col[['GSE43458.limma_Cov']]=plot_col[['GSE31210.limma_Cov']]==plot_col[['TCGA_limma_trend']]="green"
  plot_col[['TCGA_limma']]=sequential(plot_col[['TCGA_limma_trend']])[[15]]
  plot_col[['TCGA_edgeR']]=sequential(plot_col[['TCGA_limma_trend']])[[10]]
  plot_col[['TCGA_DESeq2']]=sequential(plot_col[['TCGA_limma_trend']])[[5]]
  plot_col[['Pseudobulk_limma_trend']]='blue'
  plot_col[['Pseudobulk_limma']]=sequential(plot_col[['Pseudobulk_limma_trend']])[[15]]
  plot_col[['Pseudobulk_edgeR']]=sequential(plot_col[['Pseudobulk_limma_trend']])[[10]]
  plot_col[['Pseudobulk_DESeq2']]=sequential(plot_col[['Pseudobulk_limma_trend']])[[5]]
  plot_col<-plot_col[order(match(names(plot_col),default.order))]
  
  
  
  axis_order<-intersect(default.order,col_meths.use%>%sapply(FUN=function(x)return(select_bk_df(x)))%>%as.vector())

  
  
  res = data.frame(Methods = meths, Metric=cutoff.metric, Cutoff=cutoff, FP=values, Seed=seeds)
  res$Methods%<>%factor()
  res$Methods<-factor(res$Methods,levels=axis_order)
  res$Metric%<>%factor()
  res$Metric<-factor(res$Metric,levels=sort(unique(res$Metric)))
  res$Cutoff%<>%factor()
  res$Cutoff<-factor(res$Cutoff,levels=rev(sort(unique(res$Cutoff))))

  res2 = melt(res, measure.vars=c("FP"))
  

  # dummy <- data.frame(Cutoff = unique(sort(res2$Cutoff)), Metric=rep('pval',2), hl=(unique(sort(res2$Cutoff))%>%gsub(pattern='cutoff = ', replacement = '')%>%as.numeric())*ngene)
  # dummy <- data.frame(Cutoff = unique(sort(res2$Cutoff)), Metric=rep('pval',1), hl=(unique(sort(res2$Cutoff))%>%gsub(pattern='cutoff = ', replacement = '')%>%as.numeric())*ngene)
  dummy <- data.frame(Cutoff = unique(sort(res2$Cutoff)), Metric=rep(unique(res2$Metric)[str_detect(unique(res2$Metric),pattern='p-values')],1), hl=(thr.use*ngene))
  # dummy <- data.frame(Cutoff = unique(sort(res2$Cutoff)), Metric=rep('pval',length(unique(res2$Cutoff))), hl=(unique(sort(res2$Cutoff))%>%gsub(pattern='cutoff = ', replacement = '')%>%as.numeric())*ngene)
  # dummy <- data.frame(Cutoff = unique(c('cutoff = 0.01', 'cutoff = 0.05')), Metric=rep('pval',2), hl=(unique(sort(c('cutoff = 0.01', 'cutoff = 0.05')))%>%gsub(pattern='cutoff = ', replacement = '')%>%as.numeric())*ngene)
  dummy$Metric%<>%factor(levels=sort(unique(res2$Metric)))
  pd = position_dodge(width=0.0)
  print(sum(res2$value[which(res2$Cutoff=='q-values < 0.05')]>300))
  gbase = ggplot(res2, aes(y=value, x=Methods, fill=Methods))+
    geom_boxplot(position=pd,outlier.shape=1)+
    # facet_grid(Cutoff~Metric, scales='free')+
    facet_grid(Metric~Cutoff, scales='free')+
    # facet_grid(Metric~Cutoff)+
    geom_hline(data = dummy, aes(yintercept=hl), linetype="dashed", color = "red", size=1)+
    # geom_hline(yintercept=ngene*0.01, linetype="dashed", color = "red")+
    scale_x_discrete(limits = axis_order)+
    # scale_y_continuous(limits=c(0,0.25*ngene))+
    scale_y_continuous(limits=c(0,3500))+
    # coord_cartesian(ylim=c(0,2000))+
    scale_fill_manual(values = plot_col[which(names(plot_col)%ni%c('base_line'))])+
    guides(colour=guide_legend(nrow=4, byrow = T))+
    theme(plot.title = element_text(size=40),
          axis.text.x=element_text(angle=60, hjust=1, size = 35), 
          axis.text.y=element_text(size=35), 
          axis.title.x=element_text(size=0),
          axis.title.y=element_text(size=40),
          strip.text=element_text(size=40),
          legend.position = 'none')
  # scale_colour_manual(name = "Methods",
  #                     labels = ts[order(ts)],
  #                     values = tc[order(ts)])
  
  gline = gbase
  print(gline)
  # tt = paste(ref,' / ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),' / min.cell = ',filter.rate*100,' percent ',sep="")
  tt = paste('FP analysis / normal ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),' / ',filter.rate*100,' percent min.cell threshold ',sep="")
  
  print(gline+aes(x=Methods)+labs(x='Methods', y='FP')+ggtitle(tt))
  figurename=paste0(ref,'.',up.genes,' | ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),'gene threshold : expressed in  ',filter.rate*100,'percent cells | N_DEG.pdf')
  ggsave(file=paste(fig_dir,"/",figurename,sep=""),width = 25,height = 16)
  dev.off()
  # for(thr in c(0.01,0.05)){
  #   for(metric in c('pval','adj.pval')){
  #     
  #   }
  # }
}

draw_plot<-function(base_dir, saved_dir,meth.comb.tot,meth.comb.use,ct,filter.rate,adj.p_thr,ref,ref_genes,ref_weight,sub.ct,filt_dir,fig_dir,sort.meth,up.genes='up', plot_topn_genes=c(1000)){
  library(pracma)
  meth.comb.tot<-unique(unlist(meth.comb.use))
  fig_dir.orig<-fig_dir
  TCGA.exist=any(str_detect(pattern = 'TCGA',unique(unlist(meth.comb.use))))

  default.order=c('base_line','TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend','GSE43458.limma_Cov','GSE31210.limma_Cov','Pseudobulk_DESeq2','Pseudobulk_edgeR','Pseudobulk_limma','Pseudobulk_limma_trend','Combat_Wilcox','limma_BEC_Wilcox','MNNCorrect_Wilcox','scMerge_Wilcox','Seurat_Wilcox','ZW_Wilcox','Raw_Wilcox','MAST','MAST_Cov','DESeq2','DESeq2_Cov','ZW_DESeq2','ZW_DESeq2_Cov','DESeq2_Fisher','DESeq2_wFisher','edgeR_DetRate','edgeR_DetRate_Cov','edgeR','edgeR_Cov','ZW_edgeR','ZW_edgeR_Cov','edgeR_Fisher','edgeR_wFisher','limma','limma_Cov','limmatrend','limmatrend_Cov','Combat_limmatrend','MNNCorrect_limmatrend','scMerge_limmatrend','LogNorm+limmatrend_Fisher','LogNorm+limmatrend_wFisher','DESeq2_FEM','voom_FEM','LogNorm_FEM','DESeq2_REM','voom_REM','LogNorm_REM','voom+modt_Fisher','voom+modt_wFisher')
  
  
  plot_df<-get_plot_df_saved(saved_dir,ct=ct,meth.comb.tot=unique(unlist(meth.comb.use)),ref_genes=ref_genes, ref_weight=ref_weight,base_dir=base_dir,filter.rate=filter.rate,filt_dir=filt_dir,adj.p_thr=adj.p_thr,sub.ct=sub.ct,sort.meth=sort.meth,up.genes=up.genes)
  # df_col<-data.frame(name=unique(plot_df$methods), color=c(colorRampPalette(brewer.pal(9,'Set1'))(plot_df$methods%>%unique()%>%length()-1),"#A6D854"))
  
  
  
  
  library(RColorBrewer)
  library(ggplot2)
  ggplotColours <- function(n = 6, h = c(0, 360) + 15){
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
  }
  `%ni%`<-Negate(`%in%`)
  
  plot_col<-ggplotColours(n=64)[1:length(default.order)]
  names(plot_col)<-default.order
  

  plot_df$methods<-plot_df$methods%>%sapply(FUN=function(x)return(select_bk_df(x)))%>%as.vector()
  # plot_col<-ggplotColours(n=length(intersect(plot_df$methods,default.order)))[1:length(intersect(plot_df$methods,default.order))]
  # names(plot_col)<-intersect(default.order,plot_df$methods)
  
  col_meths.use<-unique(plot_df$methods)[!str_detect(unique(plot_df$methods),pattern='ind_P')]
  
  
  
  ind_patients<-unique(plot_df$methods)[str_detect(unique(plot_df$methods),pattern = 'ind_P')]
  # plot_col_sub<-c(rep('grey',length(ind_patients)),"#A6D854")
  plot_col_sub<-c(rep('grey',length(ind_patients)),"black")
  names(plot_col_sub)=c(ind_patients, col_meths.use[str_detect(col_meths.use,pattern='[/]')])
  
  plot_col=c(plot_col,plot_col_sub)
  plot_df%<>%dplyr::filter(methods%in%names(plot_col))
  
  
  
  plot_lty=rep('solid',length(plot_col))
  plot_size=rep(1.8,length(plot_col))
  names(plot_lty)=names(plot_size)=names(plot_col)
  # plot_lty[ind_patients]=1
  plot_lty[col_meths.use[str_detect(col_meths.use,pattern='[/]')]]='longdash'
  plot_size[ind_patients]=0.9
  plot_size[col_meths.use[str_detect(col_meths.use,pattern='[/]')]]=1.5
  
  if(any(str_detect(names(plot_col),pattern='ind_P'))){
    library(colortools)
    
    plot_col[['GSE43458.limma_Cov']]=plot_col[['GSE31210.limma_Cov']]==plot_col[['TCGA_limma_trend']]="green"
    plot_col[['TCGA_limma']]=sequential(plot_col[['TCGA_limma_trend']])[[15]]
    plot_col[['TCGA_edgeR']]=sequential(plot_col[['TCGA_limma_trend']])[[10]]
    plot_col[['TCGA_DESeq2']]=sequential(plot_col[['TCGA_limma_trend']])[[5]]
    
    plot_col[['Pseudobulk_limma_trend']]='blue'
    plot_col[['Pseudobulk_limma']]=sequential(plot_col[['Pseudobulk_limma_trend']])[[15]]
    plot_col[['Pseudobulk_edgeR']]=sequential(plot_col[['Pseudobulk_limma_trend']])[[10]]
    plot_col[['Pseudobulk_DESeq2']]=sequential(plot_col[['Pseudobulk_limma_trend']])[[5]]
    plot_lty[['Pseudobulk_limma']]=plot_lty[['TCGA_limma']]="11"
    plot_lty[['Pseudobulk_edgeR']]=plot_lty[['TCGA_edgeR']]="21"
    plot_lty[['Pseudobulk_DESeq2']]=plot_lty[['TCGA_DESeq2']]="4121"
    
    
    
    plot_col[['DESeq2']]=plot_col[['edgeR']]=plot_col[['limmatrend']]='#CD9600'
    plot_col[['DESeq2_Cov']]=plot_col[['edgeR_Cov']]=plot_col[['limmatrend_Cov']]='#ABA300'
    plot_col[['edgeR_DetRate']]=plot_col[['limma']]='#7CAE00'
    plot_col[['edgeR_DetRate_Cov']]=plot_col[['limma_Cov']]='#0CB702'
    
    # plot_col[str_detect(names(plot_col),pattern='Combat')]='#00BC58'
    # plot_col[str_detect(names(plot_col),pattern='MNNCorrect')]='#00C082'
    # plot_col[str_detect(names(plot_col),pattern='scMerge')]='#00C0BB'
    # plot_col[str_detect(names(plot_col),pattern='limma_BEC')]='#00BADF'
    
    plot_col[str_detect(names(plot_col),pattern='ZW_')]='#00C082' #'#00AEFA'
    plot_col[str_detect(names(plot_col),pattern='ZW.*_Cov')]='#49A0FF'
    plot_col[str_detect(names(plot_col),pattern='_REM')]='#AA88FF'
    plot_col[str_detect(names(plot_col),pattern='_FEM')]='#DD71FA'
    plot_col[str_detect(names(plot_col),pattern='_Fisher')]='#F07E4C' # '#FC61D5'
    plot_col[str_detect(names(plot_col),pattern='_wFisher')]='#FF68A1'
  }else if(any(str_detect(pattern='Main_figure',string=names(meth.comb.use)))){
    library(colortools)
    plot_col[['GSE43458.limma_Cov']]=plot_col[['GSE31210.limma_Cov']]==plot_col[['TCGA_limma_trend']]="green"
    plot_col[['TCGA_limma']]=sequential(plot_col[['TCGA_limma_trend']])[[15]]
    plot_col[['TCGA_edgeR']]=sequential(plot_col[['TCGA_limma_trend']])[[10]]
    plot_col[['TCGA_DESeq2']]=sequential(plot_col[['TCGA_limma_trend']])[[5]]
    
    plot_col[['Pseudobulk_limma_trend']]='blue'
    plot_col[['Pseudobulk_limma']]=sequential(plot_col[['Pseudobulk_limma_trend']])[[15]]
    plot_col[['Pseudobulk_edgeR']]=sequential(plot_col[['Pseudobulk_limma_trend']])[[10]]
    plot_col[['Pseudobulk_DESeq2']]=sequential(plot_col[['Pseudobulk_limma_trend']])[[5]]
    plot_lty[['Pseudobulk_limma']]=plot_lty[['TCGA_limma']]="11"
    plot_lty[['Pseudobulk_edgeR']]=plot_lty[['TCGA_edgeR']]="21"
    plot_lty[['Pseudobulk_DESeq2']]=plot_lty[['TCGA_DESeq2']]="4121"
    
    plot_col[['Seurat_Wilcox']]='#F07E4C'
    plot_col[['Raw_Wilcox']]='#F07E4C'
    # plot_col[['MNNCorrect_Wilcox']]='#CD9600'
    # plot_col[['ZINB-WaVE_Wilcox']]='#ABA300'
    plot_col[['ZW_Wilcox']]='#CD9600'
    plot_col[['MAST']]='#FEBD1A'
    plot_col[['DESeq2']]='#0CB702'
    plot_col[['ZW_DESeq2']]='#7CAE00'
    plot_col[['DESeq2_wFisher']]='#00C0BB'
    plot_col[['DESeq2_FEM']]='#00AEFA'
    plot_col[['limmatrend']]='#AA88FF'
    plot_col[['MNNCorrect_limmatrend']]='#DD71FA'
    plot_col[['edgeR']]='#FC61D5'
    plot_col[['ZW_edgeR']]='#FF68A1'
  }else if(any(str_detect(pattern='findmarker',string=names(meth.comb.use)))){
    library(colortools)
    plot_col[['GSE43458.limma_Cov']]=plot_col[['GSE31210.limma_Cov']]==plot_col[['TCGA_limma_trend']]="green"
    plot_col[['TCGA_limma']]=sequential(plot_col[['TCGA_limma_trend']])[[15]]
    plot_col[['TCGA_edgeR']]=sequential(plot_col[['TCGA_limma_trend']])[[10]]
    plot_col[['TCGA_DESeq2']]=sequential(plot_col[['TCGA_limma_trend']])[[5]]
    
    plot_col[['Pseudobulk_limma_trend']]='blue'
    plot_col[['Pseudobulk_limma']]=sequential(plot_col[['Pseudobulk_limma_trend']])[[15]]
    plot_col[['Pseudobulk_edgeR']]=sequential(plot_col[['Pseudobulk_limma_trend']])[[10]]
    plot_col[['Pseudobulk_DESeq2']]=sequential(plot_col[['Pseudobulk_limma_trend']])[[5]]
    plot_lty[['Pseudobulk_limma']]=plot_lty[['TCGA_limma']]="11"
    plot_lty[['Pseudobulk_edgeR']]=plot_lty[['TCGA_edgeR']]="21"
    plot_lty[['Pseudobulk_DESeq2']]=plot_lty[['TCGA_DESeq2']]="4121"
    

    plot_col[['Raw_Wilcox']]='#F07E4C'
    # plot_col[['MNNCorrect_Wilcox']]='#CD9600'
    # plot_col[['ZINB-WaVE_Wilcox']]='#ABA300'
    plot_col[['ZW_Wilcox']]='#CD9600'
    plot_col[['MNNCorrect_Wilcox']]='#FEBD1A'
    plot_col[['scMerge_Wilcox']]='#0CB702'
    plot_col[['Seurat_d50_Wilcox']]='#7CAE00'
    plot_col[['Seurat_d100_Wilcox']]='#00C0BB'
    plot_col[['Seurat_d150_Wilcox']]='#00AEFA'
    plot_col[['Seurat_d200_Wilcox']]='#AA88FF'
    plot_col[['MNNCorrect_limmatrend']]='#DD71FA'
    plot_col[['edgeR']]='#FC61D5'
    plot_col[['ZINB-WaVE_edgeR']]='#FF68A1'
  }
  
  
  
  
  
  
  for(topn in plot_topn_genes){
    message(topn)
    if(topn!='all'){
      topn%<>%as.numeric()
    }
    # xmax<-plot_df%>%dplyr::filter(methods %in% meth.comb.tot)%>%dplyr::select(x)%>%max()
    # ymax<-plot_df%>%dplyr::filter(methods %in% meth.comb.tot)%>%dplyr::select(y)%>%max()
    xmax<-plot_df%>%dplyr::select(x)%>%max()
    ymax<-plot_df%>%dplyr::select(y)%>%max()
    base_gradient<-plot_df%>%dplyr::filter(methods == col_meths.use[str_detect(col_meths.use,pattern='[/]')])%>%dplyr::select(y)%>%max()/plot_df%>%dplyr::filter(methods == col_meths.use[str_detect(col_meths.use,pattern='[/]')])%>%dplyr::select(x)%>%max()
    
    # if(xmax>ymax/base_gradient){
    #   ymax<-xmax*base_gradient
    # }else{
    #   xmax<-ymax/base_gradient
    # }
    
    xmax_min<-sapply(meth.comb.tot,FUN=function(m){plot_df%>%dplyr::filter(methods == m)%>%dplyr::select(x)%>%max()})%>%as.vector()
    print(meth.comb.tot[which(xmax_min==min(xmax_min))])
    print(min(xmax_min))
    xmax_min%<>%min()
    xmax_min=topn
    # xmax<-150/1.05
    # meth.comb.use<-meth.comb.tot
    # meth.comb.use<-meth.comb.tot[c(1,2,11,12,13,14)]

    
    for(meth.name in names(meth.comb.use)){
      pcol.base<-plot_col[meth.comb.use[[meth.name]]%>%sapply(FUN=function(x)return(select_bk_df(x)))%>%as.vector()]
      pcol.base<-pcol.base[!is.na(pcol.base)]
      if(any(str_detect(meth.comb.use[[meth.name]], pattern = 'ind$'))){
        plot_col.temp<-c(pcol.base,plot_col[names(plot_col)[str_detect(names(plot_col),pattern='ind_P')]],plot_col[names(plot_col)[str_detect(names(plot_col),pattern='[/]')]])
      }else{
        plot_col.temp<-c(pcol.base,plot_col[names(plot_col)[str_detect(names(plot_col),pattern='[/]')]])
      }
      plot_lty.temp=plot_lty[names(plot_col.temp)]
      plot_size.temp=plot_size[names(plot_col.temp)]
      
      # names(plot_col.temp)[str_detect(names(plot_col.temp),pattern='[/]')]='base_line'
      fig_dir<-gsub(pattern='changetomethods',replacement=meth.name,x=fig_dir.orig)
      sted<-plot_df%>%dplyr::filter(str_detect(pattern='[/]',methods))
      sted<-data.frame(x=sort(unique(sted$x)),y=sort(unique(sted$y)))
      res1<-matrix(0,ncol=(length(meth.comb.use[[meth.name]])+1),nrow=1)
      colnames(res1)=c(meth.comb.use[[meth.name]]%>%sapply(FUN=function(x)return(select_bk_df(x)))%>%as.vector(),'base_line')
      res1[,'base_line']=trapz(sted$x,sted$y)
      
      if(topn=='all'){
        message('all')
        res1[,meth.comb.use[[meth.name]]%>%sapply(FUN=function(x)return(select_bk_df(x)))%>%as.vector()]=sapply((meth.comb.use[[meth.name]]%>%sapply(FUN=function(x)return(select_bk_df(x)))%>%as.vector()), FUN=function(x){
          temp1<-plot_df[plot_df$methods==x,]
          if(nrow(temp1)==0){return(0)}else{
            if(max(temp1$x)!=max(sted$x)){
              return(trapz(x=c(temp1$x,max(sted$x)),y=c(temp1$y,max(temp1$y))))
            }else{
              return(trapz(temp1$x,temp1$y))
            }
          }
          
        })
        g1<-ggplot(plot_df%>%dplyr::filter(methods %in% c(meth.comb.use[[meth.name]]%>%sapply(FUN=function(x)return(select_bk_df(x)))%>%as.vector())),aes(x=x,y=y, color=methods,linetype=methods, size=methods))+geom_line()+coord_cartesian(xlim = c(0,xmax*1.05),ylim=c(0,xmax*1.05*base_gradient))
      }else if(topn>1){
        xmax_min%<>%as.numeric()
        message(paste0(topn, ' drawn'))
        res1[,'base_line']<-res1[,'base_line']*((topn/max(sted$x))^2)
        res1[,meth.comb.use[[meth.name]]%>%sapply(FUN=function(x)return(select_bk_df(x)))%>%as.vector()]=sapply((meth.comb.use[[meth.name]]%>%sapply(FUN=function(x)return(select_bk_df(x)))%>%as.vector()), FUN=function(x){
          temp1<-plot_df[plot_df$methods==x,]
          temp2<-temp1%>%dplyr::filter(x>topn)
          temp1%<>%dplyr::filter(x<=topn)
          
          if(nrow(temp1)==0){return(0)}else{
            if(max(temp1$x)!=topn){
              if(nrow(temp2)!=0){
                temp1<-rbind(temp1,c(topn,(((min(temp2$y)-max(temp1$y))/(min(temp2$x)-max(temp1$x)))*(topn-max(temp1$x)))+max(temp1$y),x))
              }else if(nrow(temp2)==0){
                temp1<-rbind(temp1,c(topn,max(temp1$y),x))
              }
            }
            temp1$x%<>%as.numeric()
            temp1$y%<>%as.numeric()
            return(trapz(temp1$x,temp1$y))
          }
          
        })
        g1<-ggplot(plot_df%>%dplyr::filter(methods %in% c(meth.comb.use[[meth.name]]%>%sapply(FUN=function(x)return(select_bk_df(x)))%>%as.vector())),aes(x=x,y=y, color=methods,linetype=methods, size=methods))+geom_line()+coord_cartesian(xlim = c(0,xmax_min),ylim=c(0,xmax_min*2.5*1.05*base_gradient))
      }else if(topn<1){
        xmax_min=round(as.numeric(topn)*xmax)
        message(paste0(xmax_min, ' drawn'))
        res1[,'base_line']<-res1[,'base_line']*((xmax_min/max(sted$x))^2)
        res1[,meth.comb.use[[meth.name]]%>%sapply(FUN=function(x)return(select_bk_df(x)))%>%as.vector()]=sapply((meth.comb.use[[meth.name]]%>%sapply(FUN=function(x)return(select_bk_df(x)))%>%as.vector()), FUN=function(x){
          temp1<-plot_df[plot_df$methods==x,]
          temp2<-temp1%>%dplyr::filter(x>xmax_min)
          temp1%<>%dplyr::filter(x<=xmax_min)
          if(nrow(temp1)==0){return(0)}else{
            if(max(temp1$x)!=xmax_min){
              if(nrow(temp2)!=0){
                temp1<-rbind(temp1,c(xmax_min,(((min(temp2$y)-max(temp1$y))/(min(temp2$x)-max(temp1$x)))*(xmax_min-max(temp1$x)))+max(temp1$y),x))
              }else if(nrow(temp2)==0){
                temp1<-rbind(temp1,c(xmax_min,max(temp1$y),x))
              }
            }
            temp1$x%<>%as.numeric()
            temp1$y%<>%as.numeric()
            return(trapz(temp1$x,temp1$y))
          }
          
        })
        # g1<-ggplot(plot_df%>%dplyr::filter(methods %in% c(meth.comb.use[[meth.name]]%>%sapply(FUN=function(x)return(select_bk_df(x)))%>%as.vector())),aes(x=x,y=y, color=methods))+geom_line(size=1.5)+coord_cartesian(xlim = c(0,xmax_min),ylim=c(0,xmax_min*2.5*1.05*base_gradient))
        # g1<-ggplot(plot_df%>%dplyr::filter(methods %in% c(meth.comb.use[[meth.name]]%>%sapply(FUN=function(x)return(select_bk_df(x)))%>%as.vector())),aes(x=x,y=y, color=methods,linetype=methods, size=methods))+geom_line()+coord_cartesian(xlim = c(0,xmax_min),ylim=c(0,xmax_min*2.5*1.05*base_gradient))
        if(ref=='Known disease genes'){
          g1<-ggplot(plot_df%>%dplyr::filter(methods %in% c(meth.comb.use[[meth.name]]%>%sapply(FUN=function(x)return(select_bk_df(x)))%>%as.vector())),aes(x=x,y=y, color=methods,linetype=methods, size=methods))+geom_line()+coord_cartesian(xlim = c(0,xmax_min),ylim=c(0,25))
        }else if(ref=='Prognostic genes'){
          g1<-ggplot(plot_df%>%dplyr::filter(methods %in% c(meth.comb.use[[meth.name]]%>%sapply(FUN=function(x)return(select_bk_df(x)))%>%as.vector())),aes(x=x,y=y, color=methods,linetype=methods, size=methods))+geom_line()+coord_cartesian(xlim = c(0,xmax_min),ylim=c(0,125))
        }
        
      }
      g1<-g1+geom_line(data=plot_df%>%dplyr::filter(str_detect(pattern='ind_P',methods)),aes(x=x,y=y,color=methods,linetype=methods, size=methods))
      g1<-g1+geom_line(data=plot_df%>%dplyr::filter(str_detect(pattern='[/]',methods)),aes(x=x,y=y,color=methods,linetype=methods, size=methods))
      g1<-g1+scale_linetype_manual(values = plot_lty.temp)+
        scale_colour_manual(values = plot_col.temp)+
        scale_size_manual(values = plot_size.temp)+
        guides(guide_legend(nrow = 3, byrow = TRUE))
      print(g1+labs(x=paste0(sort.meth, ' ranked DE genes'), y='Gold standard genes')+ggtitle(paste0(ref,' | ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),'\n','gene threshold : expressed in  ',filter.rate*100,'% cells | adj.p_threshold : ',adj.p_thr))+theme(axis.text=element_text(size=35), axis.title = element_text(size=40),plot.title = element_text(size=40),legend.title = element_text(size=20),legend.text = element_text(lineheight = .8, size=20),legend.margin=margin(t=1, unit='cm'), legend.direction = 'horizontal', legend.position = 'bottom'))
      if((topn!='all')&(topn<1)){
        dir.create(paste0(fig_dir,'/',topn*100,'percent_genes/'),showWarnings = F)
        print(paste0(fig_dir,'/',topn*100,'percent_genes/',ref,'.',up.genes,' | ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),'gene threshold : expressed in  ',filter.rate*100,'percent cells | adj.p_threshold : ',adj.p_thr, '|',sort.meth, ' sorted.pdf'))
        ggsave(file=paste0(fig_dir,'/',topn*100,'percent_genes/',ref,'.',up.genes,' | ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),'gene threshold : expressed in  ',filter.rate*100,'percent cells | adj.p_threshold : ',adj.p_thr, '|',sort.meth, ' sorted.pdf'), width = 20, height=16)
        write.table(res1,paste0(fig_dir,'/',topn*100,'percent_genes/',ref,'.',up.genes,' | ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),'gene threshold : expressed in  ',filter.rate*100,'percent cells | adj.p_threshold : ',adj.p_thr, '|',sort.meth, ' sorted_pauc.txt'))
      }else{
        dir.create(paste0(fig_dir,'/',topn,'_genes/'),showWarnings = F)
        print(paste0(fig_dir,'/',topn,'_genes/',ref,'.',up.genes,' | ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),'gene threshold : expressed in  ',filter.rate*100,'percent cells | adj.p_threshold : ',adj.p_thr, '|',sort.meth, ' sorted.pdf'))
        ggsave(file=paste0(fig_dir,'/',topn,'_genes/',ref,'.',up.genes,' | ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),'gene threshold : expressed in  ',filter.rate*100,'percent cells | adj.p_threshold : ',adj.p_thr, '|',sort.meth, ' sorted.pdf'), width = 20, height=16)
        write.table(res1,file=paste0(fig_dir,'/',topn,'_genes/',ref,'.',up.genes,' | ',ct,ifelse(is.null(sub.ct),'',paste0('-',sub.ct)),'gene threshold : expressed in  ',filter.rate*100,'percent cells | adj.p_threshold : ',adj.p_thr, '|',sort.meth, ' sorted_pauc.txt'))
      }
      print(paste0(fig_dir,'/',topn,ref))
      dev.off()
    }
  }
}

get_plot_points<-function(paths,pvalues,ref_paths, meth.anal){
  if(anyNA(pvalues)){
    message(paste0(sum(is.na(pvalues)),' NAs in pvalues'))
    pvalues[is.na(pvalues)]=1
  }
  sort=pvalues
  sort<-sapply(sort,FUN=function(x)which(x==(sort%>%unique()%>%sort(decreasing = T))))+1
  si<-sapply(c(1:length(sort)), FUN = function(x){i=x
  if(sum(sort==sort[i])>1){i=(which(sort==sort[i])%>%max())}
  return(i)})
  si%<>%unique()
  temp_x<-c(0,si)
  temp_y<-c(0,sapply(si,FUN=function(x){sum(ref_paths %in% intersect(paths[1:x],ref_paths))}))
  
  temp_plots<-data.frame(x=temp_x,y=temp_y,methods=rep(meth.anal,length(temp_x)))
  return(temp_plots)
}

draw_plot_paths<-function(fn,gene.used,meth.anal.use,ref_paths,is.TCGA=T,sort.meth, result.list,xmax=50,xbase=xbase,ybase=ybase){
  
  library(pracma)
  res_plots=data.frame(x=c(),y=c(),methods=c())
  
  if(is.TCGA){
    for(meth.anal in c('TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend')){
      res_plots%<>%rbind(get_plot_points(paths=result.list[[sort.meth]][['TCGA']][[meth.anal]]$pathways,pvalues=result.list[[sort.meth]][['TCGA']][[meth.anal]]$pvalue,ref_paths, meth.anal=paste0('whole_',meth.anal)))
    }
  }
  for(meth.anal in meth.anal.use){
    res_plots%<>%rbind(get_plot_points(paths=result.list[[sort.meth]][[gene.used]][[meth.anal]]$pathways,pvalues=result.list[[sort.meth]][[gene.used]][[meth.anal]]$pvalue,ref_paths, meth.anal=meth.anal))
  }
  
  default.order=c('Raw_Wilcox','ZW_Wilcox','MAST','MAST_Cov','DESeq2','DESeq2_Cov','ZW_DESeq2','ZW_DESeq2_Cov','DESeq2_Fisher','DESeq2_wFisher','DESeq2_FEM','edgeR','edgeR_Cov','ZW_edgeR','ZW_edgeR_Cov','limmatrend','limmatrend_Cov','MNNCorrect_limmatrend','Combat_limmatrend','LogNorm_FEM','TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend',paste0('whole_',c('TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend')),'base_line')
  default.order2=c('Raw_Wilcox','ZW_Wilcox','limmatrend','limmatrend_Cov','MNNCorrect_limmatrend','Combat_limmatrend','MAST','MAST_Cov','DESeq2','DESeq2_Cov','ZW_DESeq2','ZW_DESeq2_Cov','DESeq2_Fisher','DESeq2_wFisher','DESeq2_FEM','edgeR','edgeR_Cov','ZW_edgeR','ZW_edgeR_Cov','LogNorm_FEM','TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend',paste0('whole_',c('TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend')),'base_line')
  
  library(RColorBrewer)
  library(ggplot2)
  ggplotColours <- function(n = 6, h = c(0, 360) + 15){
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
  }
  `%ni%`<-Negate(`%in%`)
  
  plot_col<-ggplotColours(n=64)[1:length(default.order)]
  names(plot_col)<-default.order
  library(scales)
  pal <- ggplotColours(n=64)
  col_meths.use<-unique(res_plots$methods)
  
  plot_col[['base_line']]='black'
  res_plots%<>%dplyr::filter(methods%in%names(plot_col))
  
  plot_lty=rep('solid',length(plot_col))
  plot_size=rep(2.1,length(plot_col))
  names(plot_lty)=names(plot_size)=names(plot_col)
  
  {
    library(colortools)
    plot_col[['TCGA_limma']]=sequential("green")[[10]]
    plot_col[['TCGA_limma_trend']]=sequential(plot_col[['TCGA_limma']])[[15]]
    plot_col[['TCGA_DESeq2']]=sequential(plot_col[['TCGA_limma']])[[10]]
    plot_col[['TCGA_edgeR']]=sequential(plot_col[['TCGA_limma']])[[5]]
    
    plot_size[['TCGA_limma']]=plot_size[['TCGA_limma_trend']]=plot_size[['TCGA_DESeq2']]=plot_size[['TCGA_edgeR']]=1.6
    
    plot_lty[['base_line']]="42"
    plot_col[['whole_TCGA_DESeq2']]='blue'
    plot_col[['whole_TCGA_edgeR']]=sequential(plot_col[['whole_TCGA_DESeq2']])[[15]]
    plot_col[['whole_TCGA_limma']]=sequential(plot_col[['whole_TCGA_DESeq2']])[[10]]
    plot_col[['whole_TCGA_limma_trend']]=sequential(plot_col[['whole_TCGA_DESeq2']])[[5]]
    plot_lty[['whole_TCGA_edgeR']]=plot_lty[['TCGA_limma']]="11"
    plot_lty[['whole_TCGA_limma']]=plot_lty[['TCGA_DESeq2']]="21"
    plot_lty[['whole_TCGA_limma_trend']]=plot_lty[['TCGA_edgeR']]="4121"
    
    plot_col[['Raw_Wilcox']]='blue'
    plot_col[['MAST']]='#FEBD1A'
    plot_col[['MAST_Cov']]='#FEBD1A'
    
    plot_col[['DESeq2']]='#ABA300'
    plot_col[['DESeq2_Cov']]='#ABA300'
    plot_col[['ZW_DESeq2']]='#964800'
    plot_col[['ZW_DESeq2_Cov']]='#964800'
    plot_col[['DESeq2_Fisher']]='#00C1A6'
    plot_col[['DESeq2_wFisher']]='#00BCD7'
    plot_col[['DESeq2_FEM']]='#00AEFA'
    plot_col[['limmatrend']]="#AA88FF"
    plot_col[['limmatrend_Cov']]="#AA88FF"
    plot_col[['edgeR']]='#FC61D5'
    plot_col[['edgeR_Cov']]='#FC61D5'
    plot_col[['ZW_edgeR']]='#FF68A1'
    plot_col[['ZW_edgeR_Cov']]='red'
  }
  # xmax<-50
  ymax<-res_plots%>%dplyr::filter(x<=xmax)%>%dplyr::select(y)%>%max()
  if(any(str_detect(pattern='[/]',res_plots$methods))){
    sted<-res_plots%>%dplyr::filter(str_detect(pattern='[/]',methods))
    sted<-data.frame(x=sort(unique(sted$x)),y=sort(unique(sted$y)))
    res1<-matrix(0,ncol=(length(unique(res_plots$methods))+1),nrow=1)
    colnames(res1)=c(setdiff(unique(res_plots$methods),'base_line'),'base_line')
    res1[,'base_line']=trapz(sted$x,sted$y)
  }else{
    res1<-matrix(0,ncol=length(unique(res_plots$methods))+1,nrow=1)
    colnames(res1)=c(unique(res_plots$methods),'base_line')
    res1[,'base_line']=trapz(x=c(0,xbase),y=c(0,ybase))*((xmax/xbase)^2)
  }
  res_plots<-rbind(res_plots,data.frame(x=c(0,xmax,xbase),y=c(0,ybase*xmax/xbase,ybase),methods=rep('base_line',3)))
  plot_col<-plot_col[names(plot_col)%in%res_plots$methods]
  plot_lty<-plot_lty[names(plot_lty)%in%res_plots$methods]
  plot_size<-plot_size[names(plot_size)%in%res_plots$methods]
  xmax%<>%as.numeric()
  message(paste0(xmax, ' drawn'))
  res1[,setdiff(unique(res_plots$methods),'base_line')]=sapply(setdiff(unique(res_plots$methods),'base_line'), FUN=function(x){
    temp1<-res_plots[res_plots$methods==x,]
    temp2<-temp1%>%dplyr::filter(x>xmax)
    temp1%<>%dplyr::filter(x<=xmax)
    
    temp1$x%<>%as.numeric()
    temp1$y%<>%as.numeric()
    if(max(temp1$x)<xmax){
      if(nrow(temp2)==0){
        return(trapz(x=c(temp1$x,xmax),y=c(temp1$y,max(temp1$y))))
      }else{
        return(trapz(x=c(temp1$x,xmax)),y=c(temp1$y,((min(temp2$y)-max(temp1$y))/(min(temp2$x)-max(temp1$x)))*(xmax-max(temp1$x))))
      }
    }else{
      return(trapz(temp1$x,temp1$y))
    }
  })
  
  
  res_plots_ed_point<-data.frame(x=c(),y=c(),methods=c())
  for(meth in unique(res_plots$methods)){
    res_plots_temp<-res_plots%>%dplyr::filter(methods==meth)
    res_plots_temp<-res_plots_temp[which(res_plots_temp$x==max(res_plots_temp$x)),]
    res_plots_ed_point%<>%rbind(res_plots_temp)
  }
  
  g1<-ggplot(data=res_plots%>%dplyr::filter(methods %ni% c('TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend','base_line')),aes(x=x,y=y, color=methods,linetype=methods, size=methods))+geom_line()
  g1<-g1+geom_line(data=res_plots%>%dplyr::filter(methods %in% c('TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend')),aes(x=x,y=y, color=methods,linetype=methods, size=methods))
  
  g1<-g1+geom_line(data=res_plots%>%dplyr::filter(methods %in% c('base_line')),aes(x=x,y=y, color=methods,linetype=methods, size=methods))+coord_cartesian(xlim = c(0,xmax),ylim=c(0,ymax))
  
  g1<-g1+scale_linetype_manual(values = plot_lty)+
    scale_colour_manual(values = plot_col)+
    scale_size_manual(values = plot_size)+
    guides(guide_legend(nrow = 3, byrow = TRUE))
  
  
  
  print(g1+labs(x=paste0(sort.meth, ' ranked paths'), y='Gold standard paths')+theme(axis.text=element_text(size=35), axis.title = element_text(size=40),plot.title = element_text(size=40),legend.title = element_text(size=20),legend.text = element_text(lineheight = .8, size=20),legend.margin=margin(t=1, unit='cm'), legend.direction = 'horizontal', legend.position = 'bottom'))
  
  print(paste0(fn,'_plot.pdf'))
  ggsave(file=paste0(fn,'_plot.pdf'), width = 20, height=16)
  write.table(res1,paste0(fn,'_pauc.txt'))
  dev.off()
  
  res_plots2<-data.frame(pauc=res1[1,],methods=colnames(res1))
  gbase = ggplot(res_plots2%>%dplyr::filter(methods %ni% c('base_line')), aes(y=pauc, x=methods, fill=methods))+geom_bar(stat = 'identity')+
    geom_hline(data = res_plots2%>%dplyr::filter(methods %in% c('base_line')), aes(yintercept=pauc), size=2, linetype="dashed", color = "black")+
    scale_x_discrete(limits = intersect(default.order2,setdiff(res_plots2$methods,c('base_line'))))+
    theme(axis.text.x=element_text(angle=60, hjust=1), legend.direction = 'horizontal', legend.position = 'none')+
    scale_fill_manual(values = plot_col)
  gline = gbase
  print(gline+labs(x='Methods', y='pAUC')+theme(axis.text=element_text(size=30), axis.title = element_text(size=30),plot.title = element_text(size=30)))
  ggsave(file=paste0(fn,'_plot_auc.point.pdf'), width = 20, height=16)
  dev.off()
}

get.methcomblist<-function(without=c()){
  tcga<-c('raw++TCGA_DESeq2_Cov','raw++TCGA_edgeR_Cov','raw++TCGA_voom.limma_Cov','raw++TCGA_limma_trend_Cov')
  pseudobulk<-c('raw++pseudo_DESeq2','raw++pseudo_edgeR','raw++pseudo_voom.limma','raw++pseudo_limma_trend')
  meth.comb.list=list()
  # meth.comb.list[['ind']]=c('raw+ind.DESeq2.pval+ind','raw+ind.edgeR.pval+ind')
  # meth.comb.list[['findmarker']]=c('raw++findmarkers','zinbwave++findmarkers','combat++findmarkers', 'scMerge++findmarkers', 'mnn_opt++findmarkers','limma_bec++findmarkers' , "Seurat++findmarkers",tcga,pseudobulk)
  # meth.comb.list[['findmarker+ind']]=c(meth.comb.list[['findmarker']],tcga,pseudobulk,'raw+ind.findmarkers+ind')
  # meth.comb.list[['DESeq2']]=c('raw++DESeq2','raw++DESeq2_Cov','raw+ind.DESeq2.ES+REM','raw+ind.DESeq2.ES+FEM',"raw+ind.DESeq2.pval+Fisher",'raw+ind.DESeq2.pval+wFisher','zinbwave++DESeq2_pseudo','zinbwave++DESeq2_pseudo_Cov')
  # meth.comb.list[['DESeq2+ind']]=c(meth.comb.list[['DESeq2']],tcga,pseudobulk,'raw+ind.DESeq2.pval+ind')
  
  meth.comb.list[['edgeR']]=  c('raw++edgeR','raw++edgeR_Cov', "raw+ind.edgeR.pval+Fisher", 'raw+ind.edgeR.pval+wFisher', 'zinbwave++edgeR', 'zinbwave++edgeR_Cov','raw++edgeR_DetRate', 'raw++edgeR_DetRate_Cov')
  meth.comb.list[['edgeR+ind']]=c(meth.comb.list[['edgeR']],tcga,pseudobulk,'raw+ind.edgeR.pval+ind')
  


  meth.comb.list[['Total']]=c('combat++findmarkers','limma_bec++findmarkers','mnn_opt++findmarkers','scMerge++findmarkers',"Seurat_++findmarkers", 'zinbwave++findmarkers', 'raw++findmarkers','raw++MAST','raw++MAST_Cov','raw++DESeq2','raw++DESeq2_Cov', 'zinbwave++DESeq2_pseudo','zinbwave++DESeq2_pseudo_Cov' ,"raw+ind.DESeq2.pval+Fisher",'raw+ind.DESeq2.pval+wFisher', 'raw++edgeR_DetRate', 'raw++edgeR_DetRate_Cov','raw++edgeR','raw++edgeR_Cov', 'zinbwave++edgeR', 'zinbwave++edgeR_Cov' ,"raw+ind.edgeR.pval+Fisher", 'raw+ind.edgeR.pval+wFisher', "raw++limma", "raw++limma_Cov", 'raw++limma_trend', 'raw++limma_trend_Cov','combat++limma_trend' , 'mnn_opt++limma_trend','scMerge++limma_trend', 'LogNormalize+ind.limma_trend+Fisher', 'LogNormalize+ind.limma_trend+wFisher','raw+ind.DESeq2.ES+FEM', "voom+ind.ES+FEM","LogNormalize+ind.ES+FEM", "raw+ind.DESeq2.ES+REM",  "voom+ind.ES+REM", "LogNormalize+ind.ES+REM", 'voom+ind.modt+Fisher', 'voom+ind.modt+wFisher',pseudobulk)
  meth.comb.list[['Total+TCGA']]=c(meth.comb.list[['Total']],tcga)
  meth.comb.list[['Total+GSE43458']]=c(meth.comb.list[['Total']],'raw++LUAD-GSE43458.microarray.voom.limma.cov')
  meth.comb.list[['Total+GSE31210']]=c(meth.comb.list[['Total']],'raw++LUAD-GSE31210.microarray.voom.limma.cov')
  
  
  meth.comb.list[['Main_figure']]=c('raw++findmarkers', 'zinbwave++findmarkers','raw++MAST','raw++DESeq2', 'zinbwave++DESeq2_pseudo','raw+ind.DESeq2.ES+FEM','raw+ind.DESeq2.pval+wFisher', 'raw++edgeR', 'zinbwave++edgeR', 'raw++limma_trend', 'mnn_opt++limma_trend', tcga,pseudobulk)
  meth.comb.list[['Main_figure++GSE43458']]=c('raw++findmarkers', 'zinbwave++findmarkers','raw++MAST','raw++DESeq2', 'zinbwave++DESeq2_pseudo','raw+ind.DESeq2.ES+FEM','raw+ind.DESeq2.pval+wFisher', 'raw++edgeR', 'zinbwave++edgeR', 'raw++limma_trend', 'mnn_opt++limma_trend','raw++LUAD-GSE43458.microarray.voom.limma.cov',pseudobulk)
  meth.comb.list[['Main_figure++GSE31210']]=c('raw++findmarkers', 'zinbwave++findmarkers','raw++MAST','raw++DESeq2', 'zinbwave++DESeq2_pseudo','raw+ind.DESeq2.ES+FEM','raw+ind.DESeq2.pval+wFisher', 'raw++edgeR', 'zinbwave++edgeR', 'raw++limma_trend', 'mnn_opt++limma_trend','raw++LUAD-GSE31210.microarray.voom.limma.cov',pseudobulk)
  return(meth.comb.list2)
}

get.ctsubct<-function(ctsubct){
  sub.ct=NULL
  if(ctsubct=='M'){
    ct='Myeloid cells'
  }else if(ctsubct=='E'){
    ct="Epithelial cells"
  }else if(ctsubct=='T'){
    ct='T lymphocytes'
  }else if(ctsubct=='T.cd4'){
    ct='T lymphocytes'
    sub.ct='CD4+ Th'
  }else if(ctsubct=='Tnk'){
    ct='T-NK cells'
  }else if(ctsubct=='Tnk.cd4'){
    ct='T-NK cells'
    sub.ct='CD4+ Th'
  }else if(ctsubct=='M.momac'){
    ct='Myeloid cells'
    sub.ct='mo-Mac'
  }else if(ctsubct=='M.mono'){
    ct='Myeloid cells'
    sub.ct='Monocytes'
  }else if(ctsubct=='M.almac'){
    ct='Myeloid cells'
    sub.ct='Alveolar Mac'
  }else if(ctsubct=='M.momac.vs.almac'){
    ct='Myeloid cells'
    sub.ct="Alveolar Mac vs mo-Mac"
  }else if(ctsubct=='M.momac.vs.almac.tlung'){
    ct='Myeloid cells'
    sub.ct="tLung Alveolar Mac vs mo-Mac"
  }else if(ctsubct=='M.momac.vs.almac.nlung'){
    ct='Myeloid cells'
    sub.ct="nLung Alveolar Mac vs mo-Mac"
  }
  return(list(ct=ct,sub.ct=sub.ct))
}

