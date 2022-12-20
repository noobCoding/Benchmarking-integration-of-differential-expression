draw_plot<-function(result.list=result.list,methods_to_plot=methods_to_plot,ref=ref,ref_genes=ref_genes,ref_weight=ref_weight,
                    sort.meth=sort.meth,fdr_cut=fdr_cut,plot_topn_genes=c(0.2)){
  library(pracma)
  meth.tot<-unique(unlist(methods_to_plot))
  plot_df<-data.frame(x=numeric(),y=numeric(),methods=character())

  default.order=c('base_line','TCGA_DESeq2','TCGA_edgeR','TCGA_limmavoom','TCGA_limmatrend','GSE43458.limma_Cov','GSE10072.limma_Cov','GSE31210.limma_Cov','Pseudobulk_DESeq2','Pseudobulk_edgeR','Pseudobulk_limmavoom','Pseudobulk_limmatrend','Combat_Wilcox','limma_BEC_Wilcox','MNN_Wilcox','scMerge_Wilcox','Seurat_Wilcox','ZW_BEC_Wilcox','Raw_Wilcox','RISC_Wilcox' , "Scanorama_Wilcox", "scGen_Wilcox", "scVI_Wilcox",'MAST','MAST_Cov','DESeq2','DESeq2_Cov','ZW_DESeq2','ZW_DESeq2_Cov','DESeq2_wFisher','edgeR_DetRate','edgeR_DetRate_Cov','edgeR','edgeR_Cov','ZW_edgeR','ZW_edgeR_Cov','edgeR_wFisher','limmavoom','limmavoom_Cov','limmatrend','limmatrend_Cov','Combat_limmatrend','MNN_limmatrend','scMerge_limmatrend','LogN+limmatrend_wFisher','RISC_limmatrend',"Scanorama_limmatrend", "scGen_limmatrend", "scVI_limmatrend",'DESeq2_FEM','LogN_FEM','DESeq2_REM','LogN_REM','RISC_QP')
  
  
  
  for(meth in meth.tot){
    message(meth)
    ranked_genes<-rownames(result.list[[meth]])
    if(is.null(result.list[[meth]]$ranks)){
      si=0
      temp_x<-c(0,0)
      temp_y<-c(0,0)
    }else{
      si<-unique(result.list[[meth]]$ranks)
      temp_x<-c(0,si)
      temp_y<-c(0,sapply(si,FUN=function(x){ref_weight[ref_genes %in% intersect(ranked_genes[which(result.list[[meth]]$ranks<=x)],ref_genes)]%>%sum()}))
    }
    temp_plots<-data.frame(x=temp_x,y=temp_y,methods=rep(meth,length(temp_x)))
    plot_df%<>%rbind(temp_plots)
  }
  ed_x<-length(result.list$gene)
  ed_y<-ref_weight[ref_genes %in% intersect(result.list$gene,ref_genes)]%>%sum()
  plot_df%<>%rbind(data.frame(x=c(0,ed_x),y=c(0,ed_y),methods=rep(paste0(ed_y, '/', ed_x),2)))
  
  
  
  library(RColorBrewer)
  library(ggplot2)
  ggplotColours <- function(n = 6, h = c(0, 360) + 15){
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
  }
  `%ni%`<-Negate(`%in%`)
  
  plot_col<-ggplotColours(n=64)[1:length(default.order)]
  names(plot_col)<-default.order
  
  # col_meths.use<-unique(plot_df$methods)[!str_detect(unique(plot_df$methods),pattern='[_]ind[.]analysis[_]')]
  ind.analysis<-unique(plot_df$methods)[str_detect(unique(plot_df$methods),pattern = '[_]ind[.]analysis[_]')]
  col_meths.use<-setdiff(unique(plot_df$methods),ind.analysis)
  plot_col_sub<-c(rep('grey',length(ind.analysis)),"black")
  names(plot_col_sub)=c(ind.analysis, col_meths.use[str_detect(col_meths.use,pattern='[/]')])
  # col_meths.use<-setdiff(col_meths.use,col_meths.use[str_detect(col_meths.use,pattern='[/]')])
  

  plot_col=c(plot_col,plot_col_sub)
  plot_df%<>%dplyr::filter(methods%in%names(plot_col))
  
  plot_lty=rep('solid',length(plot_col))
  plot_size=rep(1.8,length(plot_col))
  names(plot_lty)=names(plot_size)=names(plot_col)
  plot_lty[col_meths.use[str_detect(col_meths.use,pattern='[/]')]]='longdash'
  plot_size[ind.analysis]=0.9
  plot_size[col_meths.use[str_detect(col_meths.use,pattern='[/]')]]=1.5
  # print(names(plot_col))
  if(any(str_detect(names(plot_col),pattern='[_]ind[.]analysis[_]'))){
    library(colortools)
    plot_col[['GSE43458.limma_Cov']]=plot_col[['GSE31210.limma_Cov']]=plot_col[['TCGA_limmatrend']]="green"
    plot_col[['TCGA_limmavoom']]=sequential(plot_col[['TCGA_limmatrend']])[[15]]
    plot_col[['TCGA_edgeR']]=sequential(plot_col[['TCGA_limmatrend']])[[10]]
    plot_col[['TCGA_DESeq2']]=sequential(plot_col[['TCGA_limmatrend']])[[5]]
    
    plot_col[['Pseudobulk_limmatrend']]='blue'
    plot_col[['Pseudobulk_limmavoom']]=sequential(plot_col[['Pseudobulk_limmatrend']])[[15]]
    plot_col[['Pseudobulk_edgeR']]=sequential(plot_col[['Pseudobulk_limmatrend']])[[10]]
    plot_col[['Pseudobulk_DESeq2']]=sequential(plot_col[['Pseudobulk_limmatrend']])[[5]]
    plot_lty[['Pseudobulk_limmavoom']]=plot_lty[['TCGA_limmavoom']]="11"
    plot_lty[['Pseudobulk_edgeR']]=plot_lty[['TCGA_edgeR']]="21"
    plot_lty[['Pseudobulk_DESeq2']]=plot_lty[['TCGA_DESeq2']]="4121"
    
    
    plot_col[['DESeq2']]=plot_col[['edgeR']]=plot_col[['limmatrend']]='#125B12'
    plot_col[['DESeq2_Cov']]=plot_col[['edgeR_Cov']]=plot_col[['limmatrend_Cov']]='#7CAE00'
    plot_col[['edgeR_DetRate']]=plot_col[['limmavoom']]='#CD9600'
    plot_col[['edgeR_DetRate_Cov']]=plot_col[['limmavoom_Cov']]='#00AEFA'
    
    
    plot_col[str_detect(names(plot_col),pattern='ZW_')]='#FC61D5' #
    plot_col[str_detect(names(plot_col),pattern='ZW.*_Cov')]='#FF0000'
    plot_col[str_detect(names(plot_col),pattern='_REM')]='#CD9600'
    plot_col[str_detect(names(plot_col),pattern='_FEM')]='#00AEFA'
    plot_col[str_detect(names(plot_col),pattern='_wFisher')]='#FEBD1A'
  }else if(any(str_detect(pattern='Main_figure',string=names(methods_to_plot)))){
    library(colortools)
    plot_col[['GSE43458.limma_Cov']]=plot_col[['GSE31210.limma_Cov']]=plot_col[['TCGA_limmatrend']]="green"
    plot_col[['TCGA_limmavoom']]=sequential(plot_col[['TCGA_limmatrend']])[[15]]
    plot_col[['TCGA_edgeR']]=sequential(plot_col[['TCGA_limmatrend']])[[10]]
    plot_col[['TCGA_DESeq2']]=sequential(plot_col[['TCGA_limmatrend']])[[5]]
    
    plot_col[['Pseudobulk_limmatrend']]='blue'
    plot_col[['Pseudobulk_limmavoom']]=sequential(plot_col[['Pseudobulk_limmatrend']])[[15]]
    plot_col[['Pseudobulk_edgeR']]=sequential(plot_col[['Pseudobulk_limmatrend']])[[10]]
    plot_col[['Pseudobulk_DESeq2']]=sequential(plot_col[['Pseudobulk_limmatrend']])[[5]]
    plot_lty[['Pseudobulk_limmavoom']]=plot_lty[['TCGA_limmavoom']]="11"
    plot_lty[['Pseudobulk_edgeR']]=plot_lty[['TCGA_edgeR']]="21"
    plot_lty[['Pseudobulk_DESeq2']]=plot_lty[['TCGA_DESeq2']]="4121"
    
    plot_col[['Raw_Wilcox']]='#FF0000'
    plot_col[['RISC_Wilcox']]='#CD9600'
    plot_col[['MAST']]=plot_col[['MAST_Cov']]='#FEBD1A'
    plot_col[['DESeq2']]=plot_col[['DESeq2_Cov']]='#125B12'
    plot_col[['ZW_DESeq2']]=plot_col[['ZW_DESeq2_Cov']]='#7CAE00'
    plot_col[['DESeq2_FEM']]='#00AEFA'
    plot_col[['limmatrend']]=plot_col[['limmatrend_Cov']]='#AA88FF'
    plot_col[['edgeR']]=plot_col[['edgeR_Cov']]='#F0724C'
    plot_col[['ZW_edgeR']]=plot_col[['ZW_edgeR_Cov']]='#DD71FA'
    plot_col[["LogN_FEM"]]='#CD9600'
    plot_col[['RISC_QP']]='#00C0BB'
  }
  
  for(topn in plot_topn_genes){
    message(topn)
    if(topn!='all'){
      topn%<>%as.numeric()
    }
    
    xmax<-plot_df%>%dplyr::select(x)%>%max()
    ymax<-plot_df%>%dplyr::select(y)%>%max()
    base_gradient<-plot_df%>%dplyr::filter(methods == col_meths.use[str_detect(col_meths.use,pattern='[/]')])%>%dplyr::select(y)%>%max()/plot_df%>%dplyr::filter(methods == col_meths.use[str_detect(col_meths.use,pattern='[/]')])%>%dplyr::select(x)%>%max()
    
    xmax_min<-sapply(meth.tot,FUN=function(m){plot_df%>%dplyr::filter(methods == m)%>%dplyr::select(x)%>%max()})%>%as.vector()
    print(meth.tot[which(xmax_min==min(xmax_min))])
    print(min(xmax_min))
    xmax_min%<>%min()
    xmax_min=topn
    
    for(meth.title in names(methods_to_plot)){
      pcol.base<-plot_col[methods_to_plot[[meth.title]]]
      pcol.base<-pcol.base[!is.na(pcol.base)]
      if(any(str_detect(methods_to_plot[[meth.title]], pattern = '[_]ind[.]analysis[_]'))){
        plot_col.temp<-c(pcol.base,plot_col[names(plot_col)[str_detect(names(plot_col),pattern='[_]ind[.]analysis[_]')]],plot_col[names(plot_col)[str_detect(names(plot_col),pattern='[/]')]])
      }else{
        plot_col.temp<-c(pcol.base,plot_col[names(plot_col)[str_detect(names(plot_col),pattern='[/]')]])
      }
      plot_lty.temp=plot_lty[names(plot_col.temp)]
      plot_size.temp=plot_size[names(plot_col.temp)]
      
      sted<-plot_df%>%dplyr::filter(str_detect(pattern='[/]',methods))
      sted<-data.frame(x=sort(unique(sted$x)),y=sort(unique(sted$y)))
      res1<-matrix(0,ncol=(length(methods_to_plot[[meth.title]])+1),nrow=1)
      colnames(res1)=c(methods_to_plot[[meth.title]],'base_line')
      res1[,'base_line']=trapz(sted$x,sted$y)
      
      if(topn=='all'){
        message('all')
        res1[,methods_to_plot[[meth.title]]]=sapply(methods_to_plot[[meth.title]], FUN=function(x){
          temp1<-plot_df[plot_df$methods==x,]
          if(nrow(temp1)==0){return(0)}else{
            if(max(temp1$x)!=max(sted$x)){
              return(trapz(x=c(temp1$x,max(sted$x)),y=c(temp1$y,max(temp1$y))))
            }else{
              return(trapz(temp1$x,temp1$y))
            }
          }
          
        })
        g1<-ggplot(plot_df%>%dplyr::filter(methods %in% methods_to_plot[[meth.title]]),aes(x=x,y=y, color=methods,linetype=methods, size=methods))+geom_line()+coord_cartesian(xlim = c(0,xmax*1.05),ylim=c(0,xmax*1.05*base_gradient))
      }else if(topn>1){
        xmax_min%<>%as.numeric()
        message(paste0(topn, ' drawn'))
        res1[,'base_line']<-res1[,'base_line']*((topn/max(sted$x))^2)
        res1[,methods_to_plot[[meth.title]]]=sapply(methods_to_plot[[meth.title]], FUN=function(x){
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
        g1<-ggplot(plot_df%>%dplyr::filter(methods %in% methods_to_plot[[meth.title]]),aes(x=x,y=y, color=methods,linetype=methods, size=methods))+geom_line()+coord_cartesian(xlim = c(0,xmax_min),ylim=c(0,xmax_min*2.5*1.05*base_gradient))
      }else if(topn<1){
        xmax_min=round(as.numeric(topn)*xmax)
        message(paste0(xmax_min, ' drawn'))
        res1[,'base_line']<-res1[,'base_line']*((xmax_min/max(sted$x))^2)
        res1[,methods_to_plot[[meth.title]]]=sapply(methods_to_plot[[meth.title]], FUN=function(x){
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
        g1<-ggplot(plot_df%>%dplyr::filter(methods %in% methods_to_plot[[meth.title]]),aes(x=x,y=y, color=methods))+geom_line(size=1.5)+coord_cartesian(xlim = c(0,xmax_min),ylim=c(0,xmax_min*2.5*1.05*base_gradient))
      }
      
      g1<-g1+geom_line(data=plot_df%>%dplyr::filter(str_detect(pattern='_ind.analysis_',methods)),aes(x=x,y=y,color=methods,linetype=methods, size=methods))
      g1<-g1+geom_line(data=plot_df%>%dplyr::filter(str_detect(pattern='[/]',methods)),aes(x=x,y=y,color=methods,linetype=methods, size=methods))
      g1<-g1+scale_linetype_manual(values = plot_lty.temp)+
        scale_colour_manual(values = plot_col.temp)+
        scale_size_manual(values = plot_size.temp)+
        guides(guide_legend(nrow = 3, byrow = TRUE))
      print(g1+labs(x=paste0(sort.meth, ' ranked DE genes'), y=ref)+ggtitle(paste0(ref,' detection'))+theme(axis.text=element_text(size=35), axis.title = element_text(size=40),plot.title = element_text(size=40),legend.title = element_text(size=20),legend.text = element_text(lineheight = .8, size=20),legend.margin=margin(t=1, unit='cm'), legend.direction = 'horizontal', legend.position = 'bottom'))
      res2<-t(res1)
      res2%<>%as.data.frame()
      if((topn!='all')&(topn<1)){
        ggsave(file=paste0(ref,' detected in top ',topn*100,'percent genes_',sort.meth, 'sorted.pdf'), width = 20, height=16)
        write.table(res2,paste0(ref,' detected in top ',topn*100,'percent genes_',sort.meth, 'sorted_pauc.txt'))
      }else{
        ggsave(file=paste0(ref,' detected in top ',topn,' genes_',sort.meth, 'sorted_fdr_',fdr_cut,'_cut.pdf'), width = 20, height=16)
        write.table(res2,paste0(ref,' detected in top ',topn,' genes_',sort.meth, 'sorted_fdr_',fdr_cut,'_cut_pauc.txt'))
      }
    }
  }
}

draw_plot_bar_pauc<-function(pauc_table,methods_to_plot,ref,sort.meth,top.prop=0.2, width=25, height=16, size_element.text=40){
  library(reshape2)
  library(magrittr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  default.order=c('base_line','TCGA_DESeq2','TCGA_edgeR','TCGA_limmavoom','TCGA_limmatrend','GSE43458.limma_Cov','GSE10072.limma_Cov','GSE31210.limma_Cov','Pseudobulk_DESeq2','Pseudobulk_edgeR','Pseudobulk_limmavoom','Pseudobulk_limmatrend','Combat_Wilcox','limma_BEC_Wilcox','MNN_Wilcox','scMerge_Wilcox','Seurat_Wilcox','ZW_BEC_Wilcox','Raw_Wilcox','RISC_Wilcox' , "Scanorama_Wilcox", "scGen_Wilcox", "scVI_Wilcox",'MAST','MAST_Cov','DESeq2','DESeq2_Cov','ZW_DESeq2','ZW_DESeq2_Cov','DESeq2_wFisher','edgeR_DetRate','edgeR_DetRate_Cov','edgeR','edgeR_Cov','ZW_edgeR','ZW_edgeR_Cov','edgeR_wFisher','limmavoom','limmavoom_Cov','limmatrend','limmatrend_Cov','Combat_limmatrend','MNN_limmatrend','scMerge_limmatrend','LogN+limmatrend_wFisher','RISC_limmatrend',"Scanorama_limmatrend", "scGen_limmatrend", "scVI_limmatrend",'DESeq2_FEM','LogN_FEM','DESeq2_REM','LogN_REM','RISC_QP')
  
  pauc_table<-cbind(pauc_table,rownames(pauc_table))
  colnames(pauc_table)=c('pauc','methods')
  pauc_table%<>%as.data.frame()
  pauc_table$pauc%<>%as.numeric()
  pauc_table%<>%dplyr::arrange(desc(pauc))
  base_line=pauc_table$pauc[which(pauc_table$methods=='base_line')]
  pauc_table%<>%dplyr::filter(methods%in%default.order)
  pauc_table<-cbind(pauc_table,rep(paste0('top',top.prop*100,'percent genes'),nrow(pauc_table)))
  pauc_table<-cbind(pauc_table,rep(ref,nrow(pauc_table)))
  colnames(pauc_table)=c('pauc','methods','ref','top_prop')
  
  # dummy_table<-pauc.table%>%dplyr::filter(methods=='base_line')
  # dummy_table$base_line<-top.prop*0.5
  pauc_table<-pauc.table%>%dplyr::filter(methods!='base_line')
  pauc_table%<>%dplyr::filter(methods%in%methods_to_plot)
  pauc_table$pauc<-pauc_table$pauc/base_line
  pauc_table$pauc<-pauc_table$pauc*top.prop*0.5
  
  res<-pauc_table
  
  
  library(RColorBrewer)
  library(ggplot2)
  ggplotColours <- function(n = 6, h = c(0, 360) + 15){
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
  }
  `%ni%`<-Negate(`%in%`)
  plot_col<-ggplotColours(n=64)[1:(length(default.order)-11)]
  names(plot_col)<-setdiff(default.order,c('base_line','TCGA_DESeq2','TCGA_edgeR','TCGA_limmavoom','TCGA_limmatrend','GSE43458.limma_Cov','GSE31210.limma_Cov','Pseudobulk_DESeq2','Pseudobulk_edgeR','Pseudobulk_limmavoom','Pseudobulk_limmatrend'))
  library(colortools)
  plot_col[['GSE43458.limma_Cov']]=plot_col[['GSE31210.limma_Cov']]=plot_col[['TCGA_limmatrend']]="green"
  plot_col[['TCGA_limmavoom']]=sequential(plot_col[['TCGA_limmatrend']])[[15]]
  plot_col[['TCGA_edgeR']]=sequential(plot_col[['TCGA_limmatrend']])[[10]]
  plot_col[['TCGA_DESeq2']]=sequential(plot_col[['TCGA_limmatrend']])[[5]]
  plot_col[['Pseudobulk_limmatrend']]='blue'
  plot_col[['Pseudobulk_limmavoom']]=sequential(plot_col[['Pseudobulk_limmatrend']])[[15]]
  plot_col[['Pseudobulk_edgeR']]=sequential(plot_col[['Pseudobulk_limmatrend']])[[10]]
  plot_col[['Pseudobulk_DESeq2']]=sequential(plot_col[['Pseudobulk_limmatrend']])[[5]]
  plot_col<-plot_col[order(match(names(plot_col),default.order))]
  axis_order<-intersect(default.order,unique(res$method))
  
  
  res$method%<>%factor(levels=axis_order)
  res$top.prop%<>%factor()
  res$ref%<>%factor()
  write.table(res,file=paste0('pAUC_',ref,'_',detection,'in top',top.prop*100,'percent genes.txt'))
  
  pd = position_dodge(width=0.0)
  gbase = ggplot(res, aes(y=pauc, x=method, fill=method))+geom_bar(position=pd, stat = 'identity')+
    geom_hline(yintercept=top.prop*0.5, linetype="dashed", color = "black", size=1)+
    facet_grid(ref~top_prop, scales='free')+
    scale_x_discrete(limits = axis_order)+ylim(0,0.25)+
    theme(axis.text.x=element_text(angle=60, hjust=1), legend.direction = 'horizontal', legend.position = 'none')+
    scale_fill_manual(values = plot_col[which(names(plot_col)%ni%c('base_line'))])
  gline = gbase
  print(gline)
  
  tt=paste0('pAUC_',ref,'_',detection,'in top',top.prop*100,'percent genes')
  print(gline+labs(x='Methods', y='pAUC')+ggtitle(tt)+theme(axis.text=element_text(size=size_element.text), axis.title = element_text(size=40),plot.title = element_text(size=40)))
  figurename=gsub(pattern = "[/]", replacement = "-", x = tt)
  figurename=gsub(pattern = "[|]", replacement = " ", x = figurename)
  figurename=gsub(pattern = "[\n]", replacement = " ", x = figurename)
  figurename=gsub(pattern=" = ",replacement="_",x=figurename,fixed=T)
  figurename=gsub(pattern="%",replacement="percent",x=figurename,fixed=T)
  figurename<-paste0(figurename,'.pdf')
  ggsave(file=figurename,width = width,height = height)
}

draw_plot_bar_ksp<-function(tks_table,methods_to_plot,ref,sort.meth,top.prop=0.2, width=25, height=16, size_element.text=40){
  
  library(reshape2)
  library(magrittr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  
  default.order=c('TCGA_DESeq2','TCGA_edgeR','TCGA_limmavoom','TCGA_limmatrend','GSE43458.limma_Cov','GSE10072.limma_Cov','GSE31210.limma_Cov','Pseudobulk_DESeq2','Pseudobulk_edgeR','Pseudobulk_limmavoom','Pseudobulk_limmatrend','Combat_Wilcox','limma_BEC_Wilcox','MNN_Wilcox','scMerge_Wilcox','Seurat_Wilcox','ZW_BEC_Wilcox','Raw_Wilcox','RISC_Wilcox' , "Scanorama_Wilcox", "scGen_Wilcox", "scVI_Wilcox",'MAST','MAST_Cov','DESeq2','DESeq2_Cov','ZW_DESeq2','ZW_DESeq2_Cov','DESeq2_wFisher','edgeR_DetRate','edgeR_DetRate_Cov','edgeR','edgeR_Cov','ZW_edgeR','ZW_edgeR_Cov','edgeR_wFisher','limmavoom','limmavoom_Cov','limmatrend','limmatrend_Cov','Combat_limmatrend','MNN_limmatrend','scMerge_limmatrend','LogN+limmatrend_wFisher','RISC_limmatrend',"Scanorama_limmatrend", "scGen_limmatrend", "scVI_limmatrend",'DESeq2_FEM','LogN_FEM','DESeq2_REM','LogN_REM','RISC_QP')
  
  tks_table<-cbind(tks_table,rownames(tks_table))
  colnames(tks_table)=c('tks','methods')
  tks_table%<>%as.data.frame()
  tks_table$tks%<>%as.numeric()
  tks_table%<>%dplyr::arrange(tks)
  
  tks_table<-cbind(tks_table,rep(paste0('top',top.prop*100,'percent genes'),nrow(tks_table)))
  tks_table<-cbind(tks_table,rep(ref,nrow(tks_table)))
  colnames(tks_table)=c('pauc','methods','ref','top_prop')
  
  tks_table%<>%dplyr::filter(methods%in%default.order)
  tks_table%<>%dplyr::filter(methods%in%methods_to_plot)
  tks_table$value=(-log10(tks_table$tks))
  res<-tks_table
  
  
  
  library(RColorBrewer)
  library(ggplot2)
  ggplotColours <- function(n = 6, h = c(0, 360) + 15){
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
  }
  `%ni%`<-Negate(`%in%`)
  plot_col<-ggplotColours(n=64)[1:(length(default.order)-10)]
  names(plot_col)<-setdiff(default.order,c('TCGA_DESeq2','TCGA_edgeR','TCGA_limmavoom','TCGA_limmatrend','GSE43458.limma_Cov','GSE31210.limma_Cov','Pseudobulk_DESeq2','Pseudobulk_edgeR','Pseudobulk_limmavoom','Pseudobulk_limmatrend'))
  library(colortools)
  plot_col[['GSE43458.limma_Cov']]=plot_col[['GSE31210.limma_Cov']]=plot_col[['TCGA_limmatrend']]="green"
  plot_col[['TCGA_limmavoom']]=sequential(plot_col[['TCGA_limmatrend']])[[15]]
  plot_col[['TCGA_edgeR']]=sequential(plot_col[['TCGA_limmatrend']])[[10]]
  plot_col[['TCGA_DESeq2']]=sequential(plot_col[['TCGA_limmatrend']])[[5]]
  plot_col[['Pseudobulk_limmatrend']]='blue'
  plot_col[['Pseudobulk_limmavoom']]=sequential(plot_col[['Pseudobulk_limmatrend']])[[15]]
  plot_col[['Pseudobulk_edgeR']]=sequential(plot_col[['Pseudobulk_limmatrend']])[[10]]
  plot_col[['Pseudobulk_DESeq2']]=sequential(plot_col[['Pseudobulk_limmatrend']])[[5]]
  
  plot_col<-plot_col[order(match(names(plot_col),default.order))]
  
  axis_order<-intersect(default.order,unique(res$method))
  
  
  
  res$method%<>%factor(levels=axis_order)
  res$top.prop%<>%factor()
  res$ref%<>%factor()
  
  
  
  pd = position_dodge(width=0.0)
  gbase = ggplot(res, aes(y=value, x=method, fill=method))+geom_bar(position=pd, stat = 'identity')+
    geom_hline(yintercept=log10(20), linetype="dashed", color = "red", size=1)+
    geom_hline(yintercept=2, linetype="dashed", color = "black", size=1)+
    # facet_grid(Cutoff~Metric, scales='free')+
    facet_grid(ref~top_prop, scales='free')+
    scale_x_discrete(limits = axis_order)+ylim(0,8)+
    theme(axis.text.x=element_text(angle=60, hjust=1), legend.direction = 'horizontal', legend.position = 'none')+
    scale_fill_manual(values = plot_col[which(names(plot_col)%ni%c('base_line'))])
  
  gline = gbase
  print(gline)
  tt=paste0('truncated kolmogorov-smirnov pvalues',ref,'_',detection,'in top',top.prop*100,'percent genes')
  print(gline+labs(x='Methods', y='-log10(pvalue)')+ggtitle(tt)+theme(axis.text=element_text(size=size_element.text), axis.title = element_text(size=40),plot.title = element_text(size=40)))
  figurename=gsub(pattern = "[/]", replacement = "-", x = tt)
  figurename=gsub(pattern = "[|]", replacement = " ", x = figurename)
  figurename=gsub(pattern = "[\n]", replacement = " ", x = figurename)
  figurename=gsub(pattern=" = ",replacement="_",x=figurename,fixed=T)
  figurename=gsub(pattern="%",replacement="percent",x=figurename,fixed=T)
  figurename<-paste0(figurename,'.pdf')
  
  print(figurename)
  ggsave(file=figurename,width = width,height = height)
}