library(GEOquery)
{
  `%ni%`<-Negate(`%in%`)
  geo_acc_tot<-c("GSE29013", "GSE30219", "GSE31210", "GSE50081", "GSE37745")
  #GSE29013 85 ADC 25 SQC
  #GSE30219 85 ADC 14 NTL
  #GSE31210 226 ADC 20 NTL
  #GSE50081 127 ADC 42 SQC
  #GSE37745 106 ADC 66 SQC
  library(GEOquery)
  library(stringr)
  library(dplyr)
  library(magrittr)
  library(tidyverse)
  library(survival)
  library(RegParallel)
  library(biomaRt)
  scaled_z_weight<-c()
}
for(i in 1:length(geo_acc_tot)){
  geo_acc=geo_acc_tot[i]
  print(geo_acc)
  gse <- getGEO(geo_acc,GSEMatrix=TRUE)
  gse<-gse[[1]]
  fset <- fData(gse)
  eset <- exprs(gse)
  
  eset<-do.mapping.microarray.survdata(eset,fset)
  
  pset<-pData(gse)
  
  if(geo_acc %in% c("GSE29013")){
    idx <- match(c('death_time:ch1', 'death_event:ch1', 'Stage:ch1',
                   'gender:ch1', 'age.at.surgery:ch1','histology:ch1'),colnames(pset))
    metadata <- data.frame(pset[,idx],
                           row.names = rownames(pset))
    colnames(metadata) <- c('time', 'status', 'stage',
                            'gender', 'age', 'histology')
    metadata <- metadata[!apply(metadata, 1, function(x) any( is.na(x) )),]
    metadata %<>% dplyr::filter(histology == 'Adenocarcinoma')
    scaled_z <- t(scale(t(eset[,colnames(eset)%in%rownames(metadata)])))
  }else if(geo_acc %in% c("GSE30219")){
    idx <- match(c('follow-up time (months):ch1', 'status:ch1', 'pt stage:ch1',
                   'gender:ch1', 'age at surgery:ch1','histology:ch1'),colnames(pset))
    metadata <- data.frame(pset[,idx],
                           row.names = rownames(pset))
    colnames(metadata) <- c('time', 'status', 'stage',
                            'gender', 'age', 'histology')
    
    normal.samples<-rownames(metadata)[metadata$histology=='NTL']
    tumor.samples<-rownames(metadata)[metadata$histology=='ADC']
    metadata$time<-as.numeric(metadata$time)/12
    
    metadata <- metadata[!apply(metadata, 1, function(x) any( is.na(x) )),]
    metadata %<>% dplyr::filter(histology == 'ADC')
    metadata$status<-(metadata$status=='DEAD')%>%as.numeric()
    scaled_z <- scal(eset[,colnames(eset)%in%tumor.samples],eset[,colnames(eset)%in%normal.samples])
  }else if(geo_acc %in% c("GSE31210")){
    pset%>%colnames()
    idx <- match(c("days before death/censor:ch1", "death:ch1", "pathological stage:ch1",
                   'gender:ch1', "age (years):ch1",'description'),colnames(pset))
    metadata <- data.frame(pset[,idx],
                           row.names = rownames(pset))
    colnames(metadata) <- c('time', 'status', 'stage',
                            'gender', 'age', 'histology')
    
    normal.samples<-rownames(metadata)[metadata$histology=='Gene expression data from normal lung.']
    tumor.samples<-rownames(metadata)[metadata$histology=="Gene expression data from primary lung ADC"]
    
    metadata <- metadata[!apply(metadata, 1, function(x) any( is.na(x) )),]
    metadata %<>% dplyr::filter(histology == "Gene expression data from primary lung ADC")
    metadata$status<-(metadata$status=='dead')%>%as.numeric()
    metadata$time<-as.numeric(metadata$time)/12
    scaled_z <- scal(eset[,colnames(eset)%in%tumor.samples],eset[,colnames(eset)%in%normal.samples])
  }else if(geo_acc %in% c("GSE19188")){
    pset%>%colnames()
    idx <- match(c('overall survival:ch1','status:ch1','title',
                   'gender:ch1','title',"cell type:ch1"),colnames(pset))
    
    pset[,c('overall survival:ch1','status:ch1','title','gender:ch1','title',"cell type:ch1")]
    metadata <- data.frame(pset[,idx],
                           row.names = rownames(pset))
    colnames(metadata) <- c('time', 'status', 'stage',
                            'gender', 'age', 'histology')
    metadata[which(metadata=='Not available',arr.ind = T)]=NA
    metadata[which(metadata=='<NA>',arr.ind = T)]=NA
    metadata$age=1
    metadata$stage=1
    
    normal.samples<-rownames(metadata)[metadata$histology=='healthy']
    tumor.samples<-rownames(metadata)[metadata$histology=="ADC"]
    
    metadata <- metadata[!apply(metadata, 1, function(x) any( is.na(x) )),]
    metadata %<>% dplyr::filter(histology == "ADC")
    metadata$status<-(metadata$status=='deceased')%>%as.numeric()
    metadata$time<-as.numeric(metadata$time)/12
    scaled_z <- scal(eset[,colnames(eset)%in%tumor.samples],eset[,colnames(eset)%in%normal.samples])
  }else if(geo_acc %in% c("GSE3141")){
    metadata<-t(data.frame(sapply(pset$characteristics_ch1, FUN=function(x){strsplit(x,split = '[;]')[[1]]})))
    rownames(metadata)<-rownames(pset)
    colnames(metadata)<-c('histology','time','status')
    metadata%<>%as.data.frame()
    metadata$histology%<>%gsub(pattern='.*Cell type[:] ',replacement = '')
    metadata$time%<>%gsub(pattern='.*months)[:] ',replacement = '')
    metadata$time%<>%gsub(pattern='^[>]',replacement = '')
    metadata$status%<>%gsub(pattern='.*[:] ',replacement = '')
    metadata$age=1
    metadata$stage=1
    metadata$gender=1
    
    metadata<-metadata[,c('time', 'status', 'stage',
                          'gender', 'age', 'histology')]
    
    metadata <- metadata[!apply(metadata, 1, function(x) any( is.na(x) )),]
    metadata %<>% dplyr::filter(histology == "A")
    metadata$status%<>%as.numeric()
    
    scaled_z <- t(scale(t(eset[,colnames(eset)%in%rownames(metadata)])))
  }else if(geo_acc %in% c("GSE50081")){
    pset%>%colnames()
    idx <- match(c("survival time:ch1", "status:ch1", "t-stage:ch1",
                   'Sex:ch1', "age:ch1",'histology:ch1'),colnames(pset))
    metadata <- data.frame(pset[,idx],
                           row.names = rownames(pset))
    colnames(metadata) <- c('time', 'status', 'stage',
                            'gender', 'age', 'histology')
    
    
    metadata <- metadata[!apply(metadata, 1, function(x) any( is.na(x) )),]
    metadata %<>% dplyr::filter(histology == "adenocarcinoma")
    metadata$status<-(metadata$status=='dead')%>%as.numeric()
    scaled_z <- t(scale(t(eset[,colnames(eset)%in%rownames(metadata)])))
  }else if(geo_acc %in% c("GSE37745")){
    pset%>%colnames()
    idx <- match(c("days to determined death status:ch1", "dead:ch1", "tumor stage:ch1",
                   'gender:ch1', "age:ch1",'histology:ch1'),colnames(pset))
    metadata <- data.frame(pset[,idx],
                           row.names = rownames(pset))
    colnames(metadata) <- c('time', 'status', 'stage',
                            'gender', 'age', 'histology')
    
    metadata <- metadata[!apply(metadata, 1, function(x) any( is.na(x) )),]
    metadata %<>% dplyr::filter(histology == "adeno")
    metadata$status<-(metadata$status=='yes')%>%as.numeric()
    scaled_z <- t(scale(t(eset[,colnames(eset)%in%rownames(metadata)])))
  }
  
  
  scaled_z <- scaled_z[,which(colnames(scaled_z) %in% rownames(metadata))]
  scaled_z<-scaled_z[,match(rownames(metadata),colnames(scaled_z))]
  coxdata <- data.frame(metadata, t(scaled_z))
  
  scaled_z_weight<-c(scaled_z_weight,ncol(scaled_z))
  # prepare phenotypes
  coxdata$time <- as.numeric(coxdata$time)
  coxdata$status <- as.numeric(coxdata$status)
  coxdata$age <- as.numeric(coxdata$age)
  coxdata$gender <- factor(coxdata$gender, levels = coxdata$gender%>%unique()%>%sort())
  coxdata$stage <- factor(coxdata$stage, levels = coxdata$stage%>%unique()%>%sort())
  
  
  
  library(survival)
  library(RegParallel)
  
  formula.for.use<-paste0('Surv(time, status) ~ [*] ', 
                          ifelse(coxdata$age%>%unique()%>%length()==1,'','+ age'),
                          ifelse(coxdata$gender%>%unique()%>%length()==1,'','+ gender'),
                          ifelse(coxdata$stage%>%unique()%>%length()==1,'','+ stage')
  )
  
  
  # formula.for.use='Surv(time, status) ~ [*] '
  res <- RegParallel(
    data = coxdata,
    formula = formula.for.use,
    # 'Surv(time, status) ~ [*] + age+ gender + stage '
    FUN = function(formula, data)
      coxph(formula = formula,
            data = data,
            ties = 'breslow',
            singular.ok = TRUE),
    FUNtype = 'coxph',
    variables = colnames(coxdata)[7:ncol(coxdata)],
    blocksize = 2000,
    cores = 6,
    nestedParallel = FALSE,
    conflevel = 95)
  res <- res[which(res$Term%in%colnames(coxdata)[7:ncol(coxdata)]),]
  
  if(geo_acc==geo_acc_tot[1]){
    Gene<-res$Term
    HR_table=P_table=matrix(data=NA,nrow=length(Gene),ncol=length(geo_acc_tot),dimnames = list(row=Gene,col=geo_acc_tot))
  }
  if(anyNA(match(rownames(P_table),res$Term))){
    message('NA')
  }
  P_table[,geo_acc]=res$P[match(rownames(P_table),res$Term)]
  HR_table[,geo_acc]=res$HR[match(rownames(HR_table),res$Term)]
}
#Identical weights to all pvalues(for genes) in certain dataset.
weight<-scaled_z_weight/sum(scaled_z_weight)
weight_table<-P_table
for(i in 1:ncol(P_table)){
  weight_table[,i]=rep(weight[i],nrow(P_table))
}

logHR_table<-log2(HR_table)
survival_wFisher<-get.wFisher<-function(p=P_table,weight=weight_table,logfc=logHR_table)
meta.res<-matrix(data=NA,nrow=nrow(P_table),ncol = 5+ncol(HR_table)+ncol(P_table),dimnames = list(row_names=rownames(P_table),col_names=c('stat','pval','FDR','sign',paste0(sapply(colnames(P_table),FUN=function(x)rep(x,2))%>%as.vector(),c('_P','_logHR')))))%>%as.data.frame()
meta.res$stat=survival_wFisher$stat
meta.res$pval=survival_wFisher$pval
meta.res$FDR=survival_wFisher$FDR
meta.res$sign=survival_wFisher$direction
meta.res$sign[meta.res$sign=='+']=1
meta.res$sign[meta.res$sign=='-']=-1

meta.res[,c(5+c(1:ncol(P_table)*2-1))]=P_table
meta.res[,c(5+c(1:ncol(P_table)*2))]=logHR_table
write.table(meta.res,file='survival_analysis_result.txt')

Prognostic_genes<-rownames(meta.res)[meta.res$FDR<0.05]
Prognostic_genes_weight=rep(1,length(Prognostic_genes))



do.mapping.survival<-function(eset,fset){
  eset%<>%as.data.frame()
  eset<-eset[!is.na(fset$`Gene Symbol`),]
  if(any(grep('^AFFX', rownames(eset)))){
    eset <- eset[-grep('^AFFX', rownames(eset)),]
  }
  eset$GeneName<-fset$`Gene Symbol`[match(rownames(eset),fset$ID)]
  # rownames(eset)<-fset$`Gene Symbol`[match(rownames(eset),fset$ID)]
  eset<-eset[str_detect(eset$GeneName,pattern='[///]', negate = T)%>%which(),]
  # eset<-eset[str_detect(rownames(eset),pattern='[///]', negate = T)%>%which(),]
  mart <- useMart('ENSEMBL_MART_ENSEMBL',host = 'asia.ensembl.org')
  mart <- useDataset("hsapiens_gene_ensembl", mart)
  mrna_attributes <- getBM(mart = mart,
                           attributes = c('gene_biotype',
                                          'external_gene_name'),
                           filter = 'external_gene_name',
                           values = eset$GeneName,
                           # values = rownames(eset),
                           uniqueRows = TRUE)
  mrna_attributes%<>%dplyr::filter(gene_biotype=='protein_coding')
  eset<-eset[eset$GeneName%in%mrna_attributes$external_gene_name,]
  # eset<-eset[rownames(eset)%in%mrna_attributes$external_gene_name,]
  
  
  a1<-names(table(eset$GeneName)[table(eset$GeneName)>=2])
  a2<-names(table(eset$GeneName)[table(eset$GeneName)==1])  
  # a1<-names(table(rownames(eset))[table(rownames(eset))>=2])
  # a2<-names(table(rownames(eset))[table(rownames(eset))==1])
  
  
  {
    eset.sub1<-eset[eset$GeneName%in%a1,]
    eset.sub2<-eset[eset$GeneName%in%a2,]
    # eset.sub1<-eset[rownames(eset)%in%a1,]
    # eset.sub2<-eset[rownames(eset)%in%a2,]
    # match(rownames(eset.sub1),mrna_attributes$affy_hg_u133_plus_2)
    if(nrow(eset.sub1)!=0){
      for(g in unique(a1)){
        # print(g)
        c2.t<-eset.sub1[which(eset.sub1$GeneName%in%g),]
        eset.sub2<-rbind(eset.sub2,c2.t[which.max(apply(c2.t[,setdiff(colnames(c2.t),'GeneName')],1,mean)),])
      }
    }
    eset.renew<-eset.sub2
    eset.renew%<>%as.data.frame()
    rownames(eset.renew)<-eset.renew$GeneName
    eset.renew%<>%dplyr::arrange(rownames(eset.renew))
    
    eset.renew<-eset.renew[,setdiff(colnames(eset.renew),'GeneName')]
  }
  
  return(eset.renew)
}

get.wFisher<-function(p,weight,logfc){
  k <- ncol(p)
  direction<-pval <- stat <- rep(NA, nrow(p))
  rnum <- 1:nrow(p)
  pval[rnum]<-unlist(sapply(X=rnum, FUN=function(x){metapro::wFisher(p[x,],weight=ifelse(all(weight[x,]==0),list(rep(1,ncol(weight))),list(weight[x,]))[[1]],is.onetail = F,eff.sign = c(sign(logfc)[x,]))$p}))
  direction[rnum]<-unlist(sapply(X=rnum, FUN=function(x){metapro::wFisher(p[x,],weight=ifelse(all(weight[x,]==0),list(rep(1,ncol(weight))),list(weight[x,]))[[1]],is.onetail = F,eff.sign = c(sign(logfc)[x,]))$overall.eff.direction}))
  qval <- p.adjust(pval, method = "BH")
  res <- list(stat = pval, pval = pval, FDR = qval, direction=direction)
  
  names(res$direction)<-names(res$stat) <- names(res$pval) <- names(res$FDR) <- rownames(p)
  return(res)
}
scal <- function(x,y){
  mean_n <- rowMeans(y)  # mean of normal
  sd_n <- apply(y,1,sd)  # SD of normal
  # z score as (value - mean normal)/SD normal
  res <- matrix(nrow=nrow(x), ncol=ncol(x))
  colnames(res) <- colnames(x)
  rownames(res) <- rownames(x)
  for(i in 1:dim(x)[1]){
    for(j in 1:dim(x)[2]){
      res[i,j] <- (x[i,j]-mean_n[i])/sd_n[i]
    }
  }
  return(res)
}