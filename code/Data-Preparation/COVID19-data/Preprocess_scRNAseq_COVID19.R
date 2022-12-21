base_dir='set base directory/'
setwd(base_dir)
{
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(magrittr)
  library(stringr)
  library(readr)
  fts<-read_tsv(paste0(base_dir,'/GSE158055_covid19_features.tsv.gz'))
  annot<-read_csv(paste0(base_dir,'/GSE158055_cell_annotation.csv.gz'))
  mdat<-readxl::read_xlsx(paste0(base_dir,'/GSE158055_sample_metadata.xlsx'))
  # which(is.na(mdat$`# High-throughput sequencing metadata template (version 2.1).`))
  meta.dat<-data.frame(mdat[21:304,2:ncol(mdat)])
  colnames(meta.dat)=mdat[20,2:ncol(mdat)]
  rownames(meta.dat)=mdat$`# High-throughput sequencing metadata template (version 2.1).`[21:304]
}

{
  meta.dat.2<-meta.dat%>%dplyr::filter(`characteristics: Sample time` %in% c('control','progression'))%>%
    dplyr::filter(`characteristics: CoVID-19 severity` %in% c('mild/moderate', 'severe/critical'))
  meta.dat.2%<>%dplyr::filter(`characteristics: Sample type` %in% c('frozen PBMC', 'fresh PBMC'))%>%
    dplyr::filter(`characteristics: Single cell sequencing platform`=="10X 5'")
  annot2<-annot%>%dplyr::filter(sampleID %in% meta.dat.2$title)%>%
    dplyr::filter(majorType %in% c('Mono'))
  
  #Change the file format of GSE158055_covid19_counts.mtx.gz to be read by Read10X as described in https://github.com/satijalab/seurat/issues/4030.
  #Original data is too large to be read by Read10X.
  #Both folder part1 and part contain three files (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz).
  #GSE158055_covid19_barcodes.tsv.gz is renamed to barcodes.tsv.gz of part1/part2.
  #GSE158055_covid19_features.tsv.gz is renamed to features.tsv.gz of part1/part2.
  #GSE158055_covid19_counts.mtx.gz is divided into matrix.mtx.gz of part1/part2.
  c1<-Read10X('part1/',gene.column=1)
  c1.1<-c1[,colnames(c1)%in%annot2$cellName]
  c2<-Read10X('part2/',gene.column=1)
  c2.1<-c2[,colnames(c2)%in%annot2$cellName]
  c3<-c1.1+c2.1
  
  writeMM(obj = c3, file=paste0("COVID19.raw.txt"))
  DN<-dimnames(c3)
  
  count.matrix<-readMM("COVID19.raw.txt")
  dimnames(count.matrix)=DN
  count.matrix%<>%as.matrix%<>%as.data.frame()
  
  
  annot2$Patients=meta.dat.2$Patients[match(annot2$sampleID,meta.dat.2$title)]
  
  #generate cellinfo
  Cell=annot2$cellName
  Cell<-gsub(Cell,pattern='[-]',replacement='.')
  Patients=sapply(annot2$sampleID,FUN=function(x){return(meta.dat.2$`characteristics: Sex`[which(meta.dat.2$title==x)])})%>%as.vector()
  Batch=sapply(annot2$sampleID,FUN=function(x){return(meta.dat.2$`characteristics: Sex`[which(meta.dat.2$title==x)])})%>%as.vector()
  Group=sapply(annot2$sampleID,FUN=function(x){return(meta.dat.2$`characteristics: CoVID-19 severity`[which(meta.dat.2$title==x)])})%>%as.vector()%>%gsub(pattern='[/].*',replacement='')
  cellinfo=data.frame(Cell=Cell,Batch=Batch, Group=Group, Patient=Patients)
  
  #filter COVID19 data
  dt<-CreateSeuratObject(counts=count.matrix, min.cells = 0.05*ncol(count.matrix))
  dt[["percent.mt"]] <- PercentageFeatureSet(dt, pattern = "^MT-")
  dt <- subset(dt, subset = percent.mt <= 10)
  g1<-dt@assays$RNA@counts%>%rownames()
  
  write.table(count.matrix[which(rownames(count.matrix) %in% g1),colnames(dt)],'counts.txt')
  cellinfo%<>%dplyr::filter(Cell%in%gsub(pattern='[-]',replacement='.',colnames(dt)))
  write.table(cellinfo[,c('Cell','Group','Batch')],'cellinfo.txt')
  
  #We aggregate cells for each patients in COVID19 data pseudobulk analysis
  cellinfo$Batch=cellinfo$Patients
  write.table(cellinfo[,c('Cell','Group','Batch')],'cellinfo_pseudobulk.txt')
}