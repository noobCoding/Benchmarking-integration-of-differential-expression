base_dir='set base directory/'
setwd(base_dir)
library(magrittr)
library(readr)
library(data.table)
`%ni%`=Negate(`%in%`)
# count<-read.table('/hdd2/SC lung/data/original/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz',header = T)
# annot<-read.delim2('/hdd2/SC lung/data/original/GSE131907_Lung_Cancer_cell_annotation.txt')
count<-read.table('GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz',header = T)
annot<-read.delim2('GSE131907_Lung_Cancer_cell_annotation.txt')

#select cells from normal/tumor lung tissue
annot%<>%dplyr::filter(Sample_Origin %in% c('nLung','tLung'))
#assign patient
annot$Patient=annot$Sample
annot$Patient%<>%gsub(pattern='LUNG_N',replacement = 'P00')
annot$Patient%<>%gsub(pattern='LUNG_T',replacement = 'P00')
#Exclude patients not in stage1 tumor
annot%<>%dplyr::filter(Patient %ni% c('P0031','P0028',"P0009"))

rownames(count)<-count$Index
count<-count[,c(which(colnames(count)%in%annot$Index))]

write.table(count, file='count_selected.txt')
write.table(annot, file='annot_selected.txt')

Cell.types<-c('Myeloid cells','T/NK cells','Epithelial cells')
for(ct in Cell.types){
  annot.temp<-annot%>%dplyr::filter(Cell_type.refined==ct)
  count.temp<-count[,which(colnames(count)%in%annot.temp$Index)]
  
  
  dt<-CreateSeuratObject(counts=count.temp, min.cells = 0.05*ncol(count.temp))
  dt[["percent.mt"]] <- PercentageFeatureSet(dt, pattern = "^MT-")
  dt <- subset(dt, subset = percent.mt <= 20)
  g1<-dt@assays$RNA@counts%>%rownames()
  
  #T/NK cells are renamed to T-NK cells since the '/' in filenames can cause errors.
  write.table(count.temp[which(rownames(count.temp) %in% g1),colnames(dt)],'counts.txt')
  annot.temp<-annot.temp[which(annot.temp$Index%in%colnames(dt)),]
  cellinfo.temp<-data.frame(Cell=annot.temp$Index,Group=annot.temp$Sample_Origin,Batch=annot.temp$Patient,row.names = annot.temp$Index)
  cellinfo.temp<-cellinfo.temp[colnames(dt),]
  write.table(cellinfo.temp,'cellinfo.txt')
}