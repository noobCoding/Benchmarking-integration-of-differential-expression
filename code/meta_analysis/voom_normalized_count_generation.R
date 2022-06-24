library(limma)
library(DESeq2)
library(Seurat)
library(stringr)

dir_target_base<-'~/sim_nodown3/'
res_dir<-gsub(pattern='[/]$', replacement='_processed/',x=dir_target_base)
basename(dir_target_base)

dir_targets<-dir_target_base
dir_targets<-list.files(dir_target_base,full.names = T)

for(dir_target in dir_targets){
  setwd(dir_target)
  base_names <- list.dirs('data')
  base_names%<>%setdiff('data')
  for(base_name in base_names){
    print(base_name)
    simulated_studies=list()
    dir.create(file.path(res_dir,basename(getwd()),'data',basename(base_name)), recursive=T,showWarnings = F)
    counts<- read.table(file = paste0(base_name,"/counts.txt"), sep = "\t", fill=T)
    cellinfo<-read.table(file = paste0(base_name,"/cellinfo.txt"), sep = "\t", fill=T)
    geneinfo<-read.table(file = paste0(base_name,"/geneinfo.txt"), sep = "\t", fill=T)
    count_df<-counts
    
    myFilteredData <- count_df
    rv_genes<-which(apply(myFilteredData,1,var)==0)
    rv_genes_names<-rownames(myFilteredData)[rv_genes]
    count_df<-myFilteredData[!(rownames(myFilteredData) %in% rv_genes_names),]
    geneinfo <- geneinfo[!(rownames(geneinfo) %in% rv_genes_names),]
    
    for(b in unique(cellinfo$Batch)){
      if(any(str_detect(cellinfo$Batch, pattern='Batch'))){
        simulated_studies[[b]][['x']]=as.matrix(count_df[,which(cellinfo$Batch==b)])
        b_use=b
      }else{
        simulated_studies[[paste0('Batch_',b)]][['x']]=as.matrix(count_df[,which(cellinfo$Batch==b)])
        b_use=paste0('Batch_',b)
        
      }
      y_temp<-as.vector(cellinfo$Group[which(cellinfo$Batch==b)])
      swch<-sort(unique(y_temp))
      y_temp[which(y_temp==swch[1])]=0
      y_temp[which(y_temp==swch[2])]=1
      simulated_studies[[b_use]][['y']]=as.numeric(y_temp)
    }
    save(simulated_studies, file=file.path(res_dir,basename(getwd()),'data',basename(base_name),paste0(length(unique(cellinfo$Batch)),'batch_count.RData')))
    
    for(k in 1:length(simulated_studies)){
      nf <- edgeR::calcNormFactors(simulated_studies[[k]]$x, method = 'TMM')
      voom.data <- limma::voom(simulated_studies[[k]]$x, design = model.matrix(~factor(simulated_studies[[k]]$y)), lib.size = colSums(simulated_studies[[k]]$x) * nf)
      voom.data$genes <- rownames(simulated_studies[[k]]$x)
      simulated_studies[[k]]$count<-simulated_studies[[k]]$x
      simulated_studies[[k]]$x<-voom.data$E
    }
    save(simulated_studies, file=file.path(res_dir,basename(getwd()),'data',basename(base_name),paste0(length(unique(cellinfo$Batch)),'batch_voom.RData')))
  }
  
  for(base_name in base_names){
    simulated_studies=list()
    dir.create(file.path(res_dir,basename(getwd()),'data',basename(base_name)),showWarnings = F)
    counts<- read.table(file = paste0(base_name,"/counts.txt"), sep = "\t", fill=T)
    cellinfo<-read.table(file = paste0(base_name,"/cellinfo.txt"), sep = "\t", fill=T)
    geneinfo<-read.table(file = paste0(base_name,"/geneinfo.txt"), sep = "\t", fill=T)
    count_df<-counts
    
    myFilteredData <- count_df
    rv_genes<-which(apply(myFilteredData,1,var)==0)
    rv_genes_names<-rownames(myFilteredData)[rv_genes]
    count_df<-myFilteredData[!(rownames(myFilteredData) %in% rv_genes_names),]
    geneinfo <- geneinfo[!(rownames(geneinfo) %in% rv_genes_names),]
  
    for(b in unique(cellinfo$Batch)){
      if(any(str_detect(cellinfo$Batch, pattern='Batch'))){
        simulated_studies[[b]][['x']]=as.matrix(count_df[,which(cellinfo$Batch==b)])
        b_use=b
      }else{
        simulated_studies[[paste0('Batch_',b)]][['x']]=as.matrix(count_df[,which(cellinfo$Batch==b)])
        b_use=paste0('Batch_',b)
        
      }
      y_temp<-as.vector(cellinfo$Group[which(cellinfo$Batch==b)])
      swch<-sort(unique(y_temp))
      y_temp[which(y_temp==swch[1])]=0
      y_temp[which(y_temp==swch[2])]=1
      simulated_studies[[b_use]][['y']]=as.numeric(y_temp)
    }
    
    for(k in 1:length(simulated_studies)){
      simulated_studies[[k]]$count<-simulated_studies[[k]]$x
      simulated_studies[[k]]$x<-NormalizeData(simulated_studies[[k]]$x,normalization.method = "LogNormalize", scale.factor = 10000)%>%as.matrix()
    }
    save(simulated_studies, file=file.path(res_dir,basename(getwd()),'data',basename(base_name),paste0(length(unique(cellinfo$Batch)),'batch_LogNormalize.RData')))
  }
}



