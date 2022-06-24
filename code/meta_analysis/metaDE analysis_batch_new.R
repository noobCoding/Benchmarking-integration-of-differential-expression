library(MetaDE)
library(metapro)
# Sys.setenv(JAVA_HOME='~/java/jdk1.8.0_261')
# Sys.setenv(CLASSPATH="~/java/jdk1.8.0_261/jre/lib/ext")
library(MetaQC)
library(edgeR)
library(DESeq2)
library(magrittr)
library(tidyverse)
base_dir='~/'
source(paste0(base_dir,'script_simulation_analysis/de.analysis.R'))
HVG=F

# x=ind.Res.temp
# asymptotic = T
# weight = weight


num_batch=2
dir_target_base<-'~/sim_nodown3/'
res_dir<-gsub(pattern='[/]$', replacement='_processed/',x=dir_target_base)
res_dir_targets<-list.files(res_dir,full.names = T)
for(rd in res_dir_targets){
  setwd(rd)
  base_names <- list.dirs('data')
  base_names%<>%setdiff('data')
  for(base_name in base_names){
    if(HVG==T){
      HVG_prefix='HVG_'
    }else{
      HVG_prefix=''
    }
    MetaDE.Res<-list()
    MetaDE.Res.ind<-list()
    load(file=file.path(res_dir,basename(getwd()),'data',basename(base_name),paste0(HVG_prefix,num_batch,'batch_voom.RData')))
    ndata=length(simulated_studies)
    ind.Res.voom.modt<-ind.analysis(simulated_studies,ind.method=rep("modt",ndata),nperm=300,tail='abs')
    
    for(ind.p in c('voom+modt','edgeR','deseq2','LogNormalize+limma_trend')){
      print(ind.p)
      ind.Res.temp<-list()
      if(ind.p=='voom+modt'){
        load(file=file.path(res_dir,basename(getwd()),'data',basename(base_name),paste0(HVG_prefix,num_batch,'batch_voom.RData')))
        weight<-get.weight(simulated_studies,cell.weight = 'sample',weight.gene=T)
        ind.Res.temp<-ind.Res.voom.modt
        ind.Res.temp[['logFC']]=ind.Res.temp$stat
        
        MetaDE.Res.ind[[paste0('ind.',ind.p)]]$p<-ind.Res.temp$p
        MetaDE.Res.ind[[paste0('ind.',ind.p)]]$logFC<-ind.Res.temp[['logFC']]
        MetaDE.Res.ind[[paste0('ind.',ind.p)]]$weight<-weight
      }else if(ind.p=='LogNormalize+limma_trend'){
        load(file=file.path(res_dir,basename(getwd()),'data',basename(base_name),paste0(HVG_prefix,num_batch,'batch_LogNormalize.RData')))
        
        weight<-get.weight(simulated_studies,cell.weight = 'sample',weight.gene=T)
        ind.Res.temp$p=matrix(nrow=nrow(simulated_studies[[1]]$x),ncol=length(simulated_studies))
        rownames(ind.Res.temp$p)=rownames(simulated_studies[[1]]$x)
        colnames(ind.Res.temp$p)=names(simulated_studies)
        ind.Res.temp[['logFC']]=ind.Res.temp$p
        for(ind.study in names(simulated_studies)){
          library(limma)
          library(edgeR)
          count_df<-simulated_studies[[ind.study]]$x
          count_df[is.na(count_df)] = 0.
          cellinfo<-data.frame(group=simulated_studies[[ind.study]]$y%>%factor(), 
                               row.names = colnames(count_df))
          design <- model.matrix(~group, data=cellinfo)
          
          lmfit <- lmFit(count_df, design)
          
          lmfit <- eBayes(lmfit, trend=TRUE, robust = TRUE)
          res <- topTable(lmfit, n = Inf, adjust.method = "BH", coef = 2)
          logFC<-res$logFC
          pval<-res$P.Value
          ind.Res.temp[['logFC']][,ind.study]<-logFC
          ind.Res.temp$p[,ind.study]<-pval
        }
        
        MetaDE.Res.ind[[paste0('ind.',ind.p)]]$p<-ind.Res.temp$p
        MetaDE.Res.ind[[paste0('ind.',ind.p)]]$logFC<-ind.Res.temp[['logFC']]
        MetaDE.Res.ind[[paste0('ind.',ind.p)]]$weight<-weight
      }else if(ind.p=='edgeR'){
        load(file=file.path(res_dir,basename(getwd()),'data',basename(base_name),paste0(num_batch,'batch_count.RData')))
        for(i in names(simulated_studies)){
          simulated_studies[[i]]$count=simulated_studies[[i]]$x
        }
        weight<-get.weight(simulated_studies,cell.weight = 'sample',weight.gene=T)
        ind.Res.temp$p=matrix(nrow=nrow(simulated_studies[[1]]$x),ncol=length(simulated_studies))
        rownames(ind.Res.temp$p)=rownames(simulated_studies[[1]]$x)
        colnames(ind.Res.temp$p)=names(simulated_studies)
        ind.Res.temp[['logFC']]=ind.Res.temp$p
        for(ind.study in names(simulated_studies)){
          count_df<-simulated_studies[[ind.study]]$x
          group=simulated_studies[[ind.study]]$y
          y <- DGEList(counts=count_df, group=factor(group))
          y <- calcNormFactors(y)
          
          
          
          design <- model.matrix(~factor(group))
          rownames(design) <- colnames(y)
          y <- estimateDisp(y, design, robust=TRUE)
          
          fit <- glmQLFit(y, design, robust=TRUE, prior.df = 0)
          qlf <- glmQLFTest(fit)
          logFC<-qlf$table$logFC
          pval<-qlf$table$PValue
          ind.Res.temp[['logFC']][,ind.study]<-logFC
          ind.Res.temp$p[,ind.study]<-pval
        }
        MetaDE.Res.ind[[paste0('ind.',ind.p)]]$p<-ind.Res.temp$p
        MetaDE.Res.ind[[paste0('ind.',ind.p)]]$logFC<-ind.Res.temp[['logFC']]
        MetaDE.Res.ind[[paste0('ind.',ind.p)]]$weight<-weight
        
      }else if(ind.p=='deseq2'){
        load(file=file.path(res_dir,basename(getwd()),'data',basename(base_name),paste0(num_batch,'batch_count.RData')))
        for(i in names(simulated_studies)){
          simulated_studies[[i]]$count=simulated_studies[[i]]$x
        }
        weight<-get.weight(simulated_studies,cell.weight = 'sample',weight.gene=T)
        ind.Res.temp$p=matrix(nrow=nrow(simulated_studies[[1]]$x),ncol=length(simulated_studies))
        rownames(ind.Res.temp$p)=rownames(simulated_studies[[1]]$x)
        colnames(ind.Res.temp$p)=names(simulated_studies)
        ind.Res.temp[['logFC']]=ind.Res.temp$p
        for(ind.study in names(simulated_studies)){
          count_df<-simulated_studies[[ind.study]]$x
          count_df <- round(count_df, 0) + 1 # pseudo count prevents every gene includes at least one 0
          group=simulated_studies[[ind.study]]$y
          dds <- DESeqDataSetFromMatrix(countData = count_df, colData = data.frame(group=factor(group)), design = ~group)
          dds <- DESeq2::DESeq(dds) #, fitType ='mean')
          res<-lfcShrink(dds,coef='group_1_vs_0', type="apeglm", lfcThreshold=0)
          pval=res$pvalue
          logFC=res$log2FoldChange
          ind.Res.temp[['logFC']][,ind.study]<-logFC
          ind.Res.temp$p[,ind.study]<-pval
        }
        MetaDE.Res.ind[[paste0('ind.',ind.p)]]$p<-ind.Res.temp$p
        MetaDE.Res.ind[[paste0('ind.',ind.p)]]$logFC<-ind.Res.temp[['logFC']]
        MetaDE.Res.ind[[paste0('ind.',ind.p)]]$weight<-weight
      }
      # MetaDE.Res[[paste0(ind.p,'_',"wFisher")]]=MetaDE.pvalue.bk(ind.Res.temp,asymptotic = T, weight = weight)
      # MetaDE.Res[[paste0(ind.p,'_',"Fisher")]]=MetaDE.pvalue(ind.Res.temp,meta.method=c("Fisher"),asymptotic = T)
    }
    
    for(ind.p in names(MetaDE.Res.ind)){
      wFisher_input<-list()
      wFisher_input[['weight']]<-MetaDE.Res.ind[[ind.p]]$weight
      wFisher_input[["logFC"]]<-MetaDE.Res.ind[[ind.p]]$logFC
      wFisher_input[['two_tailed']]<-MetaDE.Res.ind[[ind.p]]$p
      
      
      
      Fisher_input<-list()
      p.h=p.l=MetaDE.Res.ind[[ind.p]]$p
      for(i in colnames(p.h)){
        p.h[,i]=unlist(mapply(function(x,y) one_tailed(x=x,logFC=y,ind.tail='high'),MetaDE.Res.ind[[ind.p]]$p[,i],MetaDE.Res.ind[[ind.p]]$logFC[,i]))
        p.l[,i]=unlist(mapply(function(x,y) one_tailed(x=x,logFC=y,ind.tail='low'),MetaDE.Res.ind[[ind.p]]$p[,i],MetaDE.Res.ind[[ind.p]]$logFC[,i]))
      }
      Fisher_input[['high']]=p.h
      Fisher_input[['low']]=p.l
      Fisher_input[['logFC']]=MetaDE.Res.ind[[ind.p]]$logFC
      
      Fisher_output<-do.Meta.pval(Fisher_input,meth='Fisher')
      wFisher_output<-do.Meta.pval(wFisher_input,meth='wFisher')
      MetaDE.Res[[paste0(gsub('ind.',replacement = '',x=ind.p),'_',"Fisher")]]=Fisher_output[['two_tailed']]
      # MetaDE.Res[[paste0(gsub('ind.',replacement = '',x=ind.p),'_',"Fisher")]]$logFC=Fisher_output[['logFC']]
      MetaDE.Res[[paste0(gsub('ind.',replacement = '',x=ind.p),'_',"Fisher")]][['sign']]=sign(Fisher_output[['logFC']])
      MetaDE.Res[[paste0(gsub('ind.',replacement = '',x=ind.p),'_',"Fisher")]]=wFisher_output[['two_tailed']]
      # MetaDE.Res[[paste0(gsub('ind.',replacement = '',x=ind.p),'_',"Fisher")]]$logFC=wFisher_output[['logFC']]
      MetaDE.Res[[paste0(gsub('ind.',replacement = '',x=ind.p),'_',"Fisher")]][['sign']]=sign(Fisher_output[['logFC']])
    }
    
    for(ind.es in c('deseq2','voom','LogNormalize')){
      print(ind.es)
      ind.Res2<-list()
      if(ind.es=='deseq2'){
        load(file=file.path(res_dir,basename(getwd()),'data',basename(base_name),paste0(num_batch,'batch_count.RData')))
        temp_M<-matrix(nrow=nrow(simulated_studies[[1]]$x),ncol=length(simulated_studies))
        rownames(temp_M)=rownames(simulated_studies[[1]]$x)
        colnames(temp_M)=names(simulated_studies)
        ind.Res2$ES=ind.Res2$SE=ind.Res2$Var=temp_M
        for(ind.study in names(simulated_studies)){
          count_df<-simulated_studies[[ind.study]]$x
          count_df <- round(count_df, 0) + 1 # pseudo count prevents every gene includes at least one 0
          group=simulated_studies[[ind.study]]$y
          dds <- DESeqDataSetFromMatrix(countData = count_df, colData = data.frame(group=factor(group)), design = ~group)
          dds <- DESeq2::DESeq(dds) #, fitType ='mean')
          res<-lfcShrink(dds,coef=2, type="apeglm", lfcThreshold=0)
          ind.Res2$ES[,ind.study]=res$log2FoldChange
          ind.Res2$SE[,ind.study]=res$lfcSE
          ind.Res2$Var[,ind.study]=res$lfcSE^2
        }
      }else if(ind.es=='voom'){
        load(file=file.path(res_dir,basename(getwd()),'data',basename(base_name),paste0(HVG_prefix,num_batch,'batch_voom.RData')))
        ndata=length(simulated_studies)
        ind.Res2<-ind.cal.ES.bk(simulated_studies, vs=simulated_studies[[1]][[2]]%>%unique()%>%sort(decreasing = T))
        ind.Res2$SE=sqrt(ind.Res2$Var)
        colnames(ind.Res2$ES)=colnames(ind.Res2$Var)=colnames(ind.Res2$SE)=names(simulated_studies)
      }else if(ind.es=='LogNormalize'){
        load(file=file.path(res_dir,basename(getwd()),'data',basename(base_name),paste0(HVG_prefix,num_batch,'batch_LogNormalize.RData')))
        ndata=length(simulated_studies)
        ind.Res2<-ind.cal.ES.bk(simulated_studies, vs=simulated_studies[[1]][[2]]%>%unique()%>%sort(decreasing = T))
        ind.Res2$SE=sqrt(ind.Res2$Var)
        colnames(ind.Res2$ES)=colnames(ind.Res2$Var)=colnames(ind.Res2$SE)=names(simulated_studies)
      }
      MetaDE.Res.ind[[paste0('ind.ES.',ind.es)]]<-ind.Res2
      
      MetaDE.Res[[paste0(ind.es,"+REM")]]<-do.Meta.ES(ind.Res2,meth='REM')[['abs']]
      MetaDE.Res[[paste0(ind.es,"+REM")]][['sign']]=sign(MetaDE.Res[[paste0(ind.es,"+REM")]][['mu.hat']])
      MetaDE.Res[[paste0(ind.es,"+FEM")]]<-do.Meta.ES(ind.Res2,meth='FEM')[['abs']]
      MetaDE.Res[[paste0(ind.es,"+FEM")]][['sign']]=sign(MetaDE.Res[[paste0(ind.es,"+FEM")]][['mu.hat']])
    }
    save(MetaDE.Res,file=file.path(res_dir,basename(getwd()),'data',basename(base_name),paste0(HVG_prefix, num_batch,'batch_two_tailed_Meta_Res.RData')))
    save(MetaDE.Res.ind,file=file.path(res_dir,basename(getwd()),'data',basename(base_name),paste0(HVG_prefix, num_batch,'batch_two_tailed_Meta_Res_ind.RData')))
  }
}
