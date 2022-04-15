########## summary confusion matrix
rm(list=ls())

# load libraries
meta_method <- c(  "voom+REM", "voom+FEM",
                   "voom+modt_wFisher_sample+gene"                     ,
                   "edgeR_wFisher_sample+gene"                        ,
                   "deseq2_wFisher_sample+gene"                        ,
                   "LogNormalize+REM"                                 ,
                   "LogNormalize+FEM"                                  ,
                   "LogNormalize+limma.trend_wFisher_sample+gene"     ,
                   "voom+modt_Fisher"              ,  
                   "edgeR_Fisher"                   , 
                   "deseq2_Fisher"                   ,
                   "LogNormalize+limma.trend_Fisher" ,
                   "deseq2+FEM"                                        , "deseq2+REM"                
)

df_med <- c(      "voom+REM", "voom+FEM",
                  "LogNormalize+REM", "LogNormalize+FEM",
                  "deseq2+REM", "deseq2+FEM"
)

meta_med <- c(    "voom+modt_Fisher"              ,  
                  "edgeR_Fisher"                   , 
                  "deseq2_Fisher"                   ,
                  "LogNormalize+limma.trend_Fisher"
)

newmeta_med <- c( "voom+modt_wFisher_sample+gene"                     ,
                  "edgeR_wFisher_sample+gene"                        , 
                  "deseq2_wFisher_sample+gene"                        ,
                  "LogNormalize+limma.trend_wFisher_sample+gene"     
)

all_simu <- list()
vect_simu <- gsub('data/','',list.dirs('data', recursive=FALSE))

for (simu in vect_simu){
  DEG <- read.table(paste0('data/',simu,'/de__genes.txt'), head=T, fill = T)
  up <- read.table(paste0('data/',simu,'/true_up_genes.txt'), head=T, fill = T)
  down <- read.table(paste0('data/',simu,'/true_down_genes.txt'), head=T, fill = T)
  
  hist_001 <- c()
  perc_001 <- c()
  hist_005 <- c()
  perc_005 <- c()
  hist <- c()
  perc <- c()
  
  for (method in meta_method){
    print(method)
   
    count_001 = 0
    count_005 = 0
    count = 0
    
   
      
    if (method %in% c(df_med)){
      
      load(paste0('data/',simu,'/','4batch_righttailed_Meta_Res.RData'))
      tmp <- MetaDE.Res[[method]]
      
      if (is.null(rownames(tmp$FDR))){
        markers <- data.frame(Gene=names(tmp$FDR), pval=tmp$FDR)
      } else {
        markers <- data.frame(Gene=rownames(tmp$FDR), pval=tmp$FDR)
      }
      colnames(markers) <- c('Gene', 'pval') 
      
    }  else if (method %in% c(newmeta_med)){
      
      if (method %in% c('voom+modt_wFisher_sample+gene')){
        load(paste0('data/',simu,'/','4batch_lefttailed_Meta_Res.RData'))
      } else {
        load(paste0('data/',simu,'/','4batch_righttailed_Meta_Res.RData'))
      }
      tmp <- MetaDE.Res[[method]]
      
      if (is.null(rownames(tmp$FDR))){
        markers <- data.frame(Gene=names(tmp$FDR), pval=tmp$FDR)
      } else {
        markers <- data.frame(Gene=rownames(tmp$FDR), pval=tmp$FDR)
      }
      colnames(markers) <- c('Gene', 'pval') 
      
    }  else if (method %in% c(meta_med) ){
      if (method %in% c('voom+modt_Fisher')){
        load(paste0('data/',simu,'/','4batch_lefttailed_Meta_Res.RData'))
      } else {
        
        load(paste0('data/',simu,'/','4batch_righttailed_Meta_Res.RData'))
      }
      tmp <- MetaDE.Res[[method]]
      
      if (is.null(rownames(tmp$meta.analysis$FDR))){
        markers <- data.frame(Gene=names(tmp$meta.analysis$FDR), pval=tmp$meta.analysis$FDR)
      } else {
        markers <- data.frame(Gene=rownames(tmp$meta.analysis$FDR), pval=tmp$meta.analysis$FDR)
      }
      colnames(markers) <- c('Gene', 'pval')
    }
    markers <- markers[rownames(markers) %in% DEG$x, ]
    
    # if (method %in% c(df_med)){
      up_genes <- down$x
     
    # } else {
    #   up_genes <- up$x
    # }
    
    for (g in up$x){
      if (markers[g, 'pval'] < 0.90  && !(g %in% up_genes) ){
        if (markers[g, 'pval'] < 0.01){
          count_001 = count_001 + 1
        }
        
        if (markers[g, 'pval'] < 0.05){
          count_005 = count_005 + 1
        }
        count = count + 1
    }
  }
    
    if (method %in% c(df_med)){
        
        load(paste0('data/',simu,'/','4batch_lefttailed_Meta_Res.RData'))
        tmp <- MetaDE.Res[[method]]
        
        if (is.null(rownames(tmp$FDR))){
          markers <- data.frame(Gene=names(tmp$FDR), pval=tmp$FDR)
        } else {
          markers <- data.frame(Gene=rownames(tmp$FDR), pval=tmp$FDR)
        }
        colnames(markers) <- c('Gene', 'pval') 
        
      }  else if (method %in% c(newmeta_med)){
        
        if (method %in% c('voom+modt_wFisher_sample+gene')){
          load(paste0('data/',simu,'/','4batch_righttailed_Meta_Res.RData'))
        } else {
          load(paste0('data/',simu,'/','4batch_lefttailed_Meta_Res.RData'))
          
        }
        tmp <- MetaDE.Res[[method]]
        
        if (is.null(rownames(tmp$FDR))){
          markers <- data.frame(Gene=names(tmp$FDR), pval=tmp$FDR)
        } else {
          markers <- data.frame(Gene=rownames(tmp$FDR), pval=tmp$FDR)
        }
        colnames(markers) <- c('Gene', 'pval') 
        
      }  else if (method %in% c(meta_med) ){
        if (method %in% c('voom+modt_Fisher')){
          load(paste0('data/',simu,'/','4batch_righttailed_Meta_Res.RData'))
        } else {
          load(paste0('data/',simu,'/','4batch_lefttailed_Meta_Res.RData'))
          
        }
        tmp <- MetaDE.Res[[method]]
        
        if (is.null(rownames(tmp$meta.analysis$FDR))){
          markers <- data.frame(Gene=names(tmp$meta.analysis$FDR), pval=tmp$meta.analysis$FDR)
        } else {
          markers <- data.frame(Gene=rownames(tmp$meta.analysis$FDR), pval=tmp$meta.analysis$FDR)
        }
        colnames(markers) <- c('Gene', 'pval')
        
        markers <- markers[rownames(markers) %in% DEG$x, ]
      }
     
      
      # if (method %in% c(df_med)) {
        down_genes <- up$x
       
      # } else{
      #   down_genes <- down$x
      # }

    for (g in down$x){      
      if (markers[g, 'pval'] < 0.90 && !(g %in% down_genes)){
        
        if (markers[g, 'pval'] < 0.01){
          count_001 = count_001 + 1
        }
        
        if (markers[g, 'pval'] < 0.05){
          count_005 = count_005 + 1
        }
        count = count + 1
      }
    }
    
    hist_001 <- c(hist_001, count_001)
    perc_001 <- c(perc_001, round(100*count_001/length(DEG$x), 2) )
    
    hist_005 <- c(hist_005, count_005)
    perc_005 <- c(perc_005, round(100*count_005/length(DEG$x), 2) )
    
    hist <- c(hist, count)
    perc <- c(perc, round(100*count/length(DEG$x), 2) )
  }
  
  df <- data.frame(percent=perc, percent_005=perc_005, percent_001=perc_001,  
                   count=hist, count_005=hist_005, count_001=hist_001)
  rownames(df) <- meta_method
  all_simu[[simu]] <- df
}


saveRDS(all_simu, file = paste0('4b_meta_deg_sign_twist.rds'))