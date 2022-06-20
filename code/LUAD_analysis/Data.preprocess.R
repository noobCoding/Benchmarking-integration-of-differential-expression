library(magrittr)
library(feather)
library(readr)
library(data.table)
`%ni%`=Negate(`%in%`)
base_dir='~/'
dir.create(paste0(base_dir,'data/original'),recursive = T,showWarnings = F)
dir.create(paste0(base_dir,'data/preprocessed'),recursive = T,showWarnings = F)
count<-read.table(paste0(base_dir,'data/original','/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz'),header = T)
annot<-read.delim2('/hdd2/SC lung/data/original/GSE131907_Lung_Cancer_cell_annotation.txt')

#select cells from normal/tumor lung tissue
annot%<>%dplyr::filter(Sample_Origin %in% c('nLung','tLung'))
#assign patient
annot$Patient=annot$Sample
annot$Patient%<>%gsub(pattern='LUNG_N',replacement = 'P00')
annot$Patient%<>%gsub(pattern='LUNG_T',replacement = 'P00')
#select stage1 tumor
annot%<>%dplyr::filter(Patient %ni% c('P0031','P0028',"P0009"))

#filter count data to contain cells from normal/tumor tissue of stage 1 tumor patient
count<-count[,c(1,which(colnames(count)%in%annot$Index))]
colnames(count)[1]='GeneName'
#saved filtered data
write.table(count, file=paste0(base_dir,'data/preprocessed/GSE131907_Lung_Cancer_raw_UMI_matrix.tLung_state1-nLung.txt'))
write.table(annot, file=paste0(base_dir,'data/preprocessed/GSE131907_Lung_Cancer_cell_annotation+patient.txt'))

Cell.types<-c('Myeloid cells','T/NK cells','Epithelial cells')
for(ct in Cell.types){
  if(ct=='T/NK cells'){
    write_feather(count[,c(1,colnames(count)%in%annot$Index[which(annot$Cell_type.refined==ct)]%>%which())],path=paste0(base_dir,'data/analysis','/tLung_state1-nLung.',gsub(pattern='[/]',replacement='-',ct),'.feather'))
  }else{
    write_feather(count[,c(1,colnames(count)%in%annot$Index[which(annot$Cell_type.refined==ct)]%>%which())],path=paste0(base_dir,'data/analysis','/tLung_state1-nLung.',gsub(pattern='[/]',replacement='-',ct),'.feather'))
  }
}
