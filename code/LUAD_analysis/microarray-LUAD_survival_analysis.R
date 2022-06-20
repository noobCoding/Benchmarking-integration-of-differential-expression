library(GEOquery)
base_dir='~/'
surv_dir=paste0(base_dir,'data/analyzed/for surv/')
dir.create(surv_dir,showWarnings = F,recursive=T)
Sys.setenv("VROOM_CONNECTION_SIZE"=131072*10)
{
  `%ni%`<-Negate(`%in%`)
  geo_acc_tot.all<-c("GSE29013", "GSE30219", "GSE31210","GSE19188", "GSE3141", "GSE50081", "GSE37745")
  notinuse<-c('GSE3141','GSE19188')
  #GSE29013 85 ADC 25 SQC
  #GSE30219 85 ADC 14 NTL
  #GSE31210 226 ADC 20 NTL
  #GSE50081 127 ADC 42 SQC
  #GSE37745 106 ADC 66 SQC
  geo_acc_tot<-geo_acc_tot.all[geo_acc_tot.all%ni%notinuse]
  library(GEOquery)
  library(stringr)
  library(dplyr)
  library(magrittr)
  library(tidyverse)
  library(survival)
  library(RegParallel)
  library(biomaRt)
  scaled_z_weight<-c()
  Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)
}
for(i in 1:length(geo_acc_tot)){
  geo_acc=geo_acc_tot[i]
  print(geo_acc)
  gse <- getGEO(geo_acc,GSEMatrix=TRUE)
  gse<-gse[[1]]
  fset <- fData(gse)
  eset <- exprs(gse)
  
  eset<-do.mapping.survival(eset,fset)
  
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
  
  
  

  # idx<-idx[!is.na(idx)]
  
  
  scaled_z <- scaled_z[,which(colnames(scaled_z) %in% rownames(metadata))]
  all((colnames(scaled_z) == rownames(metadata)) == TRUE)
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
    # Surv(Time.RFS, Distant.RFS) ~ [*]
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

write.table(HR_table,paste0(surv_dir,'multi.clinic.HR.filtered.txt'))
write.table(P_table,paste0(surv_dir,'multi.clinic.P.filtered.txt'))
write.table(scaled_z_weight,paste0(surv_dir,'multi.clinic.weight.filtered.txt'))
weight<-scaled_z_weight/sum(scaled_z_weight)
logHR<-apply(log2(HR_table),1,FUN=function(x){sum(x*weight)})
surv_integrated<-get.wFisher.survival(P_table,weight = weight,eff=HR_table,is.onetail=F)
meta.res<-matrix(data=NA,nrow=nrow(P_table),ncol = 5+ncol(HR_table)+ncol(P_table),dimnames = list(row_names=rownames(P_table),col_names=c('stat','pval','FDR','sign','logHR',paste0(sapply(colnames(P_table),FUN=function(x)rep(x,2))%>%as.vector(),c('_P','_HR')))))%>%as.data.frame()
meta.res$stat=surv_integrated$stat
meta.res$pval=surv_integrated$pval
meta.res$FDR=surv_integrated$FDR
if(is.null(surv_integrated$sign)){
  meta.res$sign=NA
}else{
  meta.res$sign=surv_integrated$sign
}
meta.res$logHR=logHR
meta.res[,c(5+c(1:ncol(P_table)*2-1))]=P_table
meta.res[,c(5+c(1:ncol(P_table)*2))]=log2(HR_table)


min_positive_ind<-intersect(which(meta.res$sign=='+'),which(meta.res$logHR<0))
min_negative_ind<-intersect(which(meta.res$sign=='-'),which(meta.res$logHR>0))
meta.res$logHR[min_positive_ind]=min(abs(meta.res$logHR))*0.0001
meta.res$logHR[min_negative_ind]=-min(abs(meta.res$logHR))*0.0001

write.table(meta.res,paste0(surv_dir,'meta.clinic.integrated.filtered.txt'))


