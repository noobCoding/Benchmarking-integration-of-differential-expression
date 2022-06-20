select_topn_read = function(topn)
{
  answer = switch(topn,
                  "0.1" = paste0(as.numeric(topn)*100,'percent_genes'),
                  "0.2" = paste0(as.numeric(topn)*100,'percent_genes'),
                  'all'= 'all_genes'
  )
  return(answer)
}
select_topn = function(topn)
{
  answer = switch(topn,
                  "500" = paste0('top ',topn, ' genes'),
                  "0.1" = paste0('top ',as.numeric(topn)*100,'%', ' genes'),
                  "0.2" = paste0('top ',as.numeric(topn)*100,'%', ' genes'),
                  'all'= 'all genes'
  )
  return(answer)
}
base_dir='~/'
source(paste0(base_dir,'script_bk/visualize.sources.R'))
topn=c(0.2)
ct.sub.ct=c("Epithelial cells" ,"Myeloid cells","T-NK cells")
filter.rate=c(0.05)
sort.meth=c('pvalue')

ref=c('Known disease genes','Prognostic genes')
adj.p_thr=c('none',0.05)
params.table<-crossing(topn=topn, ct.sub.ct=ct.sub.ct, filter.rate = filter.rate, sort.meth = sort.meth, ref = ref, adj.p_thr=adj.p_thr)
write.table(params.table%>%as.data.frame(),paste0(base_dir,'data/analyzed/params.table.txt'), sep="\t")

# topn=c(0.2)
# ct.sub.ct=c("Epithelial cells" ,"Myeloid cells" ,"T-NK cells")
# filter.rate=c(0.05)
# sort.meth=c('pvalue')
# ref=c('Known disease genes','Prognostic genes')
# 
# adj.p_thr=c('none',0.05)
# params.table<-crossing(topn=topn, ct.sub.ct=ct.sub.ct, filter.rate = filter.rate, sort.meth = sort.meth, ref = ref, adj.p_thr=adj.p_thr)
# write.table(params.table%>%as.data.frame(),'/hdd2/SC lung/data/Tot_result/params.table.pauc.txt', sep="\t")



facet.x=c('topn')
facet.y=c('sort.meth-filter.rate')
# filter.list=list(a=c('ref','CTD+disgenet'),b=c('ct.sub.ct','Epithelial cells'))
# 'prognostic_new.pos'
for(f1 in c('Prognostic genes')){
  # for(f2 in c("Epithelial cells" ,"Myeloid cells" ,"Myeloid cells_Alveolar Mac" ,"Myeloid cells_Alveolar Mac vs mo-Mac" ,"Myeloid cells_mo-Mac" ,"Myeloid cells_Monocytes" ,"Myeloid cells_nLung Alveolar Mac vs mo-Mac" ,"Myeloid cells_tLung Alveolar Mac vs mo-Mac" ,"T-NK cells")){
  width=25
  height=16
  size_element.text=35
    for(f2 in c("Epithelial cells" ,"Myeloid cells" ,"T-NK cells" )){
    for(f3 in c('none',0.05)){
      for(comp.folders in c('Total',"Total+TCGA","Total+GSE43458","Total+GSE31210")){
        filter.list=list(a=c('ref',f1),b=c('ct.sub.ct',f2), d=c('adj.p_thr', f3))
        
        draw_plot_bar_ksp(dat_dir = paste0(base_dir,'data/analyzed/',comp.folders,'/'),fig_dir=paste0(base_dir,'data/analyzed/tks_fig/',comp.folders,'/'), up.genes='both',
                          facet.x=facet.x,facet.y=facet.y, filter.list=filter.list, width=width, height=height)
        
        comp.folders.convert<-c('Total+Hai_comp','Total+TCGA+Hai_comp','Total+Hai_comp_add','Total+Hai_comp_add2','Total+Hai_comp_add3')[which(comp.folders==c('sc','sc+TCGA','sc+add','sc+add2','sc+add3'))]
        draw_plot_bar_pauc(dat_dir = paste0(base_dir,'data/analyzed/gene_detection_curve/',comp.folders.convert,'/'),fig_dir=paste0(base_dir,'data/analyzed/pauc_fig/',comp.folders,'/'), up.genes='both',
                           facet.x=facet.x,facet.y=facet.y, filter.list=filter.list, width=width, height=height)
      }
    }
  }
}

draw_plot_bar_pauc<-function(base_dir,dat_dir,fig_dir,meth.comb.tot,up.genes,facet.x=c('topn'),facet.y=c('sort.meth-filter.rate'), filter.list=NULL, proportional=T, width=25, height=16, size_element.text=35){
  dir.create(dat_dir,showWarnings = F,recursive = T)
  dir.create(fig_dir,showWarnings = F,recursive = T)
  library(reshape2)
  library(magrittr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  source(paste0(base_dir,'script_bk/visualize.sources.R'))
  
  default.order=c('base_line','TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend','GSE43458.limma_Cov','GSE10072.limma_Cov','GSE31210.limma_Cov','Pseudobulk_DESeq2','Pseudobulk_edgeR','Pseudobulk_limma','Pseudobulk_limma_trend','Combat_Wilcox','limma_BEC_Wilcox','MNNCorrect_Wilcox','scMerge_Wilcox','Seurat_Wilcox','ZW_Wilcox','Raw_Wilcox','MAST','MAST_Cov','DESeq2','DESeq2_Cov','ZW_DESeq2','ZW_DESeq2_Cov','DESeq2_Fisher','DESeq2_wFisher','edgeR_DetRate','edgeR_DetRate_Cov','edgeR','edgeR_Cov','ZW_edgeR','ZW_edgeR_Cov','edgeR_Fisher','edgeR_wFisher','limma','limma_Cov','limmatrend','limmatrend_Cov','Combat_limmatrend','MNNCorrect_limmatrend','scMerge_limmatrend','LogNorm+limmatrend_Fisher','LogNorm+limmatrend_wFisher','DESeq2_FEM','voom_FEM','LogNorm_FEM','DESeq2_REM','voom_REM','LogNorm_REM','voom+modt_Fisher','voom+modt_wFisher')
  
  params.table<-read.table(paste0(base_dir,'data/analyzed/params.table.txt'))
  params.table.here<-params.table
  if(!is.null(filter.list)){
    for(f in names(filter.list)){
      f.t<-filter.list[[f]]
      params.table.here<-params.table.here%>%dplyr::filter(get(f.t[1])%in%as.character(f.t[2:length(f.t)]))
    }
  }
  
  if(length(facet.x)>1){
    params.table.here<-params.table.here%>%dplyr::filter(get(facet.x[1])%in%as.numeric(facet.x[2:length(facet.x)]))
  }
  # if(str_detect(facet.y[1], pattern='[-]')){
  #   params.table.here[[facet.y[1]]]=paste0(str_split(facet.y[1], pattern='[-]')[[1]][2],' :',params.table.here[,str_split(facet.y[1], pattern='[-]')[[1]][2]],'_',str_split(facet.y[1], pattern='[-]')[[1]][1],' :',params.table.here[,str_split(facet.y[1], pattern='[-]')[[1]][1]])
  # }
  if(str_detect(facet.y[1], pattern='[-]')){
    params.table.here[[facet.y[1]]]=paste0(params.table.here[,str_split(facet.y[1], pattern='[-]')[[1]][2]],'filt_',params.table.here[,str_split(facet.y[1], pattern='[-]')[[1]][1]],' sorted')
  }
  if(length(facet.y)>1){
    params.table.here<-de.params.table.here%>%dplyr::filter(get(facet.y[1])%in%as.numeric(facet.y[2:length(facet.y)]))
  }
  
  params.table.here[[facet.x[1]]]=factor(params.table.here[[facet.x[1]]], levels=sort(unique(params.table.here[[facet.x[1]]])))
  params.table.here[[facet.y[1]]]=factor(params.table.here[[facet.y[1]]], levels=sort(unique(params.table.here[[facet.y[1]]])))
  

  
  
  for(i in 1:nrow(params.table.here)){
    print(i)
    a<-read.table(paste0(dat_dir, params.table.here$sort.meth[i], '/',select_topn_read(as.character(params.table.here$topn[i])),'/',params.table.here$ref[i],".both | ",gsub(params.table.here$ct.sub.ct[i], pattern='[_]', replacement = '-'),"gene threshold : expressed in  ", params.table.here$filter.rate[i]*100, "percent cells | adj.p_threshold : ",params.table.here$adj.p_thr[i],'|',params.table.here$sort.meth[i], ' sorted_pauc.txt'))
    a<-t(a)
    a<-cbind(a,rownames(a))
    colnames(a)=c('pauc','methods')
    a%<>%as.data.frame()
    a$pauc%<>%as.numeric()
    a%<>%dplyr::arrange(desc(pauc))
    
    a$methods<-sapply(a$methods, FUN=function(x){select_bk_df(x)})%>%unlist()
    
    base_line=a$pauc[which(a$methods=='base_line')]
    
    a%<>%dplyr::filter(methods%in%c('TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend','GSE43458.limma_Cov','GSE31210.limma_Cov','Pseudobulk_DESeq2','Pseudobulk_edgeR','Pseudobulk_limma','Pseudobulk_limma_trend','Combat_Wilcox','limma_BEC_Wilcox','MNNCorrect_Wilcox','scMerge_Wilcox','Seurat_Wilcox','ZW_Wilcox','Raw_Wilcox','MAST','MAST_Cov','DESeq2','DESeq2_Cov','ZW_DESeq2','ZW_DESeq2_Cov','DESeq2_Fisher','DESeq2_wFisher','edgeR_DetRate','edgeR_DetRate_Cov','edgeR','edgeR_Cov','ZW_edgeR','ZW_edgeR_Cov','edgeR_Fisher','edgeR_wFisher','limma','limma_Cov','limmatrend','limmatrend_Cov','Combat_limmatrend','MNNCorrect_limmatrend','scMerge_limmatrend','LogNorm+limmatrend_Fisher','LogNorm+limmatrend_wFisher','DESeq2_FEM','voom_FEM','LogNorm_FEM','DESeq2_REM','voom_REM','LogNorm_REM','voom+modt_Fisher','voom+modt_wFisher'))
    
    rownames(a)<-a$methods
    a%<>%as.data.frame()
    temp1<-a
    temp2<-data.frame(matrix(ncol=ncol(params.table.here), nrow=nrow(temp1)))
    dummy.temp<-data.frame(matrix(ncol=(ncol(params.table.here)+1), nrow=1))
    colnames(temp2)=colnames(params.table.here)
    colnames(dummy.temp)=c('base_line',colnames(params.table.here))
    dummy.temp[,'base_line']=c(base_line)
    for(j in colnames(params.table.here)){
      temp2[[j]]=rep(params.table.here[i,j],nrow(temp1))
      dummy.temp[[j]]=c(params.table.here[i,j])
    }
    temp1<-cbind(temp1,temp2)
    if(proportional==T){
      temp1$pauc<-temp1$pauc/base_line
      if(dummy.temp$topn==0.1){
        temp1$pauc<-temp1$pauc*0.1*0.5
        dummy.temp$base_line=0.1*0.5
      }else if(dummy.temp$topn==0.2){
        temp1$pauc<-temp1$pauc*0.2*0.5
        dummy.temp$base_line=0.2*0.5
      }else if(dummy.temp$topn=='all'){
        temp1$pauc<-temp1$pauc*0.5
        dummy.temp$base_line=1*0.5
      }else{
        message('warning!!!!')
      }
    }
    if(i==1){
      res<-temp1
    }else{
      res<-rbind(res,temp1)
    }
    if(i==1){
      dummy<-dummy.temp
    }else{
      dummy<-rbind(dummy,dummy.temp)
    }
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
  plot_col[['GSE43458.limma_Cov']]=plot_col[['GSE31210.limma_Cov']]=plot_col[['TCGA_limma_trend']]="green"
  plot_col[['TCGA_limma']]=sequential(plot_col[['TCGA_limma_trend']])[[15]]
  plot_col[['TCGA_edgeR']]=sequential(plot_col[['TCGA_limma_trend']])[[10]]
  plot_col[['TCGA_DESeq2']]=sequential(plot_col[['TCGA_limma_trend']])[[5]]
  plot_col[['Pseudobulk_limma_trend']]='blue'
  plot_col[['Pseudobulk_limma']]=sequential(plot_col[['Pseudobulk_limma_trend']])[[15]]
  plot_col[['Pseudobulk_edgeR']]=sequential(plot_col[['Pseudobulk_limma_trend']])[[10]]
  plot_col[['Pseudobulk_DESeq2']]=sequential(plot_col[['Pseudobulk_limma_trend']])[[5]]
  plot_col<-plot_col[order(match(names(plot_col),default.order))]
  axis_order<-intersect(default.order,unique(res$method))
  
  
  res$method%<>%factor(levels=axis_order)
  res$topn<-sapply(res$topn,FUN=function(x){select_topn(as.character(x))})
  res$topn%<>%factor(levels=c('top 10% genes','top 20% genes','all genes'))
  dummy$topn<-sapply(dummy$topn,FUN=function(x){select_topn(as.character(x))})
  dummy$topn%<>%factor(levels=c('top 10% genes','top 20% genes','all genes'))
  for(f in unique(res$`sort.meth-filter.rate`)){
    res1<-res%>%dplyr::filter(topn=='top 10% genes')%>%filter(`sort.meth-filter.rate`==f)
    res1<-res1[match(axis_order,res1$methods),]
    h1=res1$pauc
    res1<-res%>%dplyr::filter(topn=='top 20% genes')%>%filter(`sort.meth-filter.rate`==f)
    res1<-res1[match(axis_order,res1$methods),]
    h2=res1$pauc
    res1<-res%>%dplyr::filter(topn=='all genes')%>%filter(`sort.meth-filter.rate`==f)
    res1<-res1[match(axis_order,res1$methods),]
    h3=res1$pauc
    dt<-data.frame(top10percent=h1, top20percent=h2, all_genes=h3)
    rownames(dt)=axis_order
    write.table(dt,file=paste0(fig_dir,'pAUC_',filter.list[[1]][2],'_',filter.list[[2]][2],f,'.txt'))
  }
  pd = position_dodge(width=0.0)
  gbase = ggplot(res, aes(y=pauc, x=method, fill=method))+geom_bar(position=pd, stat = 'identity')+
    geom_hline(data = dummy, aes(yintercept=base_line), linetype="dashed", color = "black",size=1)+
    # facet_grid(Cutoff~Metric, scales='free')+
    facet_grid(get(facet.y[1])~get(facet.x[1]), scales='free')+
    scale_x_discrete(limits = axis_order)+ylim(0,0.25)+
    theme(axis.text.x=element_text(angle=60, hjust=1), legend.direction = 'horizontal', legend.position = 'none')+
    scale_fill_manual(values = plot_col[which(names(plot_col)%ni%c('base_line'))])

  # tt = paste0(sapply(filter.list, FUN=function(x){paste0(x, collapse=' : ')}), collapse = '\t')
  gline = gbase
  print(gline)
  filter.list.name<-sapply(filter.list,FUN=function(x)x[1])%>%as.vector()
  filter.list.value<-sapply(filter.list,FUN=function(x)x[2])%>%as.vector()
  tt<-paste0(paste0(filter.list.name,' : ', filter.list.value),collapse = ' | ')
  tt<-paste0(tt, '\n x-axis:', facet.x[1], ' y-axis:',facet.y[1])
  print(gline+labs(x='Methods', y='pAUC')+ggtitle(tt)+theme(axis.text=element_text(size=size_element.text), axis.title = element_text(size=40),plot.title = element_text(size=40)))
  figurename=gsub(pattern = "[/]", replacement = "-", x = tt)
  figurename=gsub(pattern = "[|]", replacement = " ", x = figurename)
  figurename=gsub(pattern = "[\n]", replacement = " ", x = figurename)
  figurename=gsub(pattern=" = ",replacement="_",x=figurename,fixed=T)
  figurename=gsub(pattern="%",replacement="percent",x=figurename,fixed=T)
  figurename<-paste0(figurename,'pauc.pdf')

  ggsave(file=paste(fig_dir,figurename,sep=""),width = width,height = height)
  dev.off()
}
draw_plot_bar_ksp<-function(base_dir,dat_dir,fig_dir,meth.comb.tot,up.genes,facet.x=c('topn'),facet.y=c('sort.meth-filter.rate'), filter.list=NULL, width=25, height=16, size_element.text=35){
  dir.create(dat_dir,showWarnings = F,recursive = T)
  dir.create(fig_dir,showWarnings = F,recursive = T)
  library(reshape2)
  library(magrittr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  source(paste0(base_dir,'script_bk/visualize.sources.R'))
  meth.comb.list<-get.methcomblist()

  default.order=c('base_line','TCGA_DESeq2','TCGA_edgeR','TCGA_limma','TCGA_limma_trend','GSE43458.limma_Cov','GSE10072.limma_Cov','GSE31210.limma_Cov','Pseudobulk_DESeq2','Pseudobulk_edgeR','Pseudobulk_limma','Pseudobulk_limma_trend','Combat_Wilcox','limma_BEC_Wilcox','MNNCorrect_Wilcox','scMerge_Wilcox','Seurat_Wilcox','ZW_Wilcox','Raw_Wilcox','MAST','MAST_Cov','DESeq2','DESeq2_Cov','ZW_DESeq2','ZW_DESeq2_Cov','DESeq2_Fisher','DESeq2_wFisher','edgeR_DetRate','edgeR_DetRate_Cov','edgeR','edgeR_Cov','ZW_edgeR','ZW_edgeR_Cov','edgeR_Fisher','edgeR_wFisher','limma','limma_Cov','limmatrend','limmatrend_Cov','Combat_limmatrend','MNNCorrect_limmatrend','scMerge_limmatrend','LogNorm+limmatrend_Fisher','LogNorm+limmatrend_wFisher','DESeq2_FEM','voom_FEM','LogNorm_FEM','DESeq2_REM','voom_REM','LogNorm_REM','voom+modt_Fisher','voom+modt_wFisher')
  
  params.table<-read.table(paste0(base_dir,'data/analyzed/params.table.txt'))
  params.table.here<-params.table
  if(!is.null(filter.list)){
    for(f in names(filter.list)){
      f.t<-filter.list[[f]]
      params.table.here<-params.table.here%>%dplyr::filter(get(f.t[1])%in%as.character(f.t[2:length(f.t)]))
    }
  }
  
  if(length(facet.x)>1){
    params.table.here<-params.table.here%>%dplyr::filter(get(facet.x[1])%in%as.numeric(facet.x[2:length(facet.x)]))
  }
  # if(str_detect(facet.y[1], pattern='[-]')){
  #   params.table.here[[facet.y[1]]]=paste0(str_split(facet.y[1], pattern='[-]')[[1]][2],' :',params.table.here[,str_split(facet.y[1], pattern='[-]')[[1]][2]],'_',str_split(facet.y[1], pattern='[-]')[[1]][1],' :',params.table.here[,str_split(facet.y[1], pattern='[-]')[[1]][1]])
  # }
  if(str_detect(facet.y[1], pattern='[-]')){
    params.table.here[[facet.y[1]]]=paste0(params.table.here[,str_split(facet.y[1], pattern='[-]')[[1]][2]],'filt_',params.table.here[,str_split(facet.y[1], pattern='[-]')[[1]][1]],' sorted')
  }
  if(length(facet.y)>1){
    params.table.here<-de.params.table.here%>%dplyr::filter(get(facet.y[1])%in%as.numeric(facet.y[2:length(facet.y)]))
  }
  
  params.table.here[[facet.x[1]]]=factor(params.table.here[[facet.x[1]]], levels=sort(unique(params.table.here[[facet.x[1]]])))
  params.table.here[[facet.y[1]]]=factor(params.table.here[[facet.y[1]]], levels=sort(unique(params.table.here[[facet.y[1]]])))
  
  # if(length(facet.y)>1){
  #   params.table.here<-de.params.table.here%>%dplyr::filter(get(facet.y[1])%in%as.numeric(facet.y[2:length(facet.y)]))
  # }
  
  for(i in 1:nrow(params.table.here)){

    print(paste0(dat_dir,params.table.here$ct.sub.ct[i],'/',params.table.here$sort.meth[i],'/',params.table.here$ref[i],'.both | ',gsub(params.table.here$ct.sub.ct[i], pattern='[_]', replacement = '-'),'gene threshold : expressed in  ',params.table.here$filter.rate[i]*100,'percent cells | adj.p_threshold : ',params.table.here$adj.p_thr[i],'|',params.table.here$sort.meth[i], ' sorted.',select_topn_read(as.character(params.table.here$topn[i])),'weighted_ks.txt'))
    temp1<-read.table(paste0(dat_dir,params.table.here$ct.sub.ct[i],'/',params.table.here$sort.meth[i],'/',params.table.here$ref[i],'.both | ',gsub(params.table.here$ct.sub.ct[i], pattern='[_]', replacement = '-'),'gene threshold : expressed in  ',params.table.here$filter.rate[i]*100,'percent cells | adj.p_threshold : ',params.table.here$adj.p_thr[i],'|',params.table.here$sort.meth[i], ' sorted.',select_topn_read(as.character(params.table.here$topn[i])),'weighted_ks.txt'))
    
    
    rownames(temp1)<-rownames(temp1)%>%sapply(FUN=function(x)return(select_bk_df(x)))%>%as.vector()
    temp1%<>%dplyr::filter(rownames(temp1)%in%default.order)
    
    temp1<-data.frame(method=rownames(temp1),value=temp1$V1)
    temp1$value<-(-log10(temp1$value))
    
    temp2<-data.frame(matrix(ncol=ncol(params.table.here), nrow=nrow(temp1)))
    colnames(temp2)=colnames(params.table.here)
    for(j in colnames(params.table.here)){
      temp2[[j]]=rep(params.table.here[i,j],nrow(temp1))
    }
    temp1<-cbind(temp1,temp2)
    
    if(i==1){
      res<-temp1
    }else{
      res<-rbind(res,temp1)
    }
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
  plot_col[['GSE43458.limma_Cov']]=plot_col[['GSE31210.limma_Cov']]=plot_col[['TCGA_limma_trend']]="green"
  plot_col[['TCGA_limma']]=sequential(plot_col[['TCGA_limma_trend']])[[15]]
  plot_col[['TCGA_edgeR']]=sequential(plot_col[['TCGA_limma_trend']])[[10]]
  plot_col[['TCGA_DESeq2']]=sequential(plot_col[['TCGA_limma_trend']])[[5]]
  plot_col[['Pseudobulk_limma_trend']]='blue'
  plot_col[['Pseudobulk_limma']]=sequential(plot_col[['Pseudobulk_limma_trend']])[[15]]
  plot_col[['Pseudobulk_edgeR']]=sequential(plot_col[['Pseudobulk_limma_trend']])[[10]]
  plot_col[['Pseudobulk_DESeq2']]=sequential(plot_col[['Pseudobulk_limma_trend']])[[5]]
  
  plot_col<-plot_col[order(match(names(plot_col),default.order))]
  
  axis_order<-intersect(default.order,unique(res$method))
  
  
  
  res$method%<>%factor(levels=axis_order)
  res$topn<-sapply(res$topn,FUN=function(x){select_topn(as.character(x))})
  # res$topn%<>%factor(levels=c('top 500 genes', 'top 10% genes','top 20% genes','all genes'))
  res$topn%<>%factor(levels=c('top 10% genes','top 20% genes','all genes'))
  
  
  
  pd = position_dodge(width=0.0)
  gbase = ggplot(res, aes(y=value, x=method, fill=method))+geom_bar(position=pd, stat = 'identity')+
    geom_hline(yintercept=log10(20), linetype="dashed", color = "red", size=1)+
    geom_hline(yintercept=2, linetype="dashed", color = "black", size=1)+
    # facet_grid(Cutoff~Metric, scales='free')+
    facet_grid(get(facet.y[1])~get(facet.x[1]), scales='free')+
    scale_x_discrete(limits = axis_order)+ylim(0,8)+
    theme(axis.text.x=element_text(angle=60, hjust=1), legend.direction = 'horizontal', legend.position = 'none')+
    scale_fill_manual(values = plot_col[which(names(plot_col)%ni%c('base_line'))])

  # tt = paste0(sapply(filter.list, FUN=function(x){paste0(x, collapse=' : ')}), collapse = '\t')
  gline = gbase
  
  
  print(gline)
  filter.list.name<-sapply(filter.list,FUN=function(x)x[1])%>%as.vector()
  filter.list.value<-sapply(filter.list,FUN=function(x)x[2])%>%as.vector()
  tt<-paste0(paste0(filter.list.name,' : ', filter.list.value),collapse = ' | ')
  tt<-paste0(tt, '\n x-axis:', facet.x[1], ' y-axis:',facet.y[1])
  print(gline+labs(x='Methods', y='-log10(pvalue)')+ggtitle(tt)+theme(axis.text=element_text(size=size_element.text), axis.title = element_text(size=40),plot.title = element_text(size=40)))
  figurename=gsub(pattern = "[/]", replacement = "-", x = tt)
  figurename=gsub(pattern = "[|]", replacement = " ", x = figurename)
  figurename=gsub(pattern = "[\n]", replacement = " ", x = figurename)
  figurename=gsub(pattern=" = ",replacement="_",x=figurename,fixed=T)
  figurename=gsub(pattern="%",replacement="percent",x=figurename,fixed=T)
  figurename<-paste0(figurename,'weightedks.pvalue.pdf')

  print(figurename)
  ggsave(file=paste(fig_dir,figurename,sep=""),width = width,height = height)
  dev.off()
}
