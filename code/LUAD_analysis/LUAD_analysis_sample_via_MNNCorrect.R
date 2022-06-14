library(Seurat)
library(feather)
library(dplyr)
library(magrittr)

count.base='tLung_state1-nLung.'
ct="Epithelial cells"
# ct="Myeloid cells"
# ct="T lymphocytes.CD4+ Th"
# ct="T lymphocytes"

filter.rate=0.05 ###Also tried 0.01.
annot<-read.csv(paste0(base_dir,'/original','/GSE131907_Lung_Cancer_cell_annotation+patient.txt'))
count<-read_feather(paste0(base_dir,'/analysis','/',count.base,ct,'.feather'))

GeneName<-count$GeneName
count_matrix<-as.matrix(count[,2:ncol(count)])
rownames(count_matrix)<-GeneName
dt<-CreateSeuratObject(counts=count_matrix, min.cells = filter.rate*ncol(count_matrix))
g1<-dt@assays$RNA@counts%>%rownames()

##filtered count##
count<-count[which(count$GeneName %in% g1),]
count<-as.matrix(count[,2:ncol(count)])
rownames(count)<-GeneName
annot<-annot[,2:ncol(annot)]
annot%<>%dplyr::filter(Cell_type==ct)%>%dplyr::filter(Index %in% colnames(count))


#######mnn analysis sample#########
library(scran)
library(scales)
require(Rtsne)
library(Seurat)
library(batchelor)

sample_names<-count%>%colnames()
##annot$Patient contains patient number
##annot$Sample_Origin is consist of 'nLung' and 'tLung'
cellinfo<-data.frame(Batch=annot$Patient[match(sample_names,annot$Index)],group=annot$Sample_Origin[match(sample_names,annot$Index)],row.names = sample_names)
dt <- CreateSeuratObject(counts = count, meta.data = cellinfo)
dt <- NormalizeData(dt, normalization.method = "LogNormalize")
# t1 = Sys.time()
out.mnn.total <- batchelor::mnnCorrect(as.matrix(dt@assays$RNA@data),batch=cellinfo$Batch, k=20, sigma=1, cos.norm.in=TRUE, cos.norm.out=TRUE, var.adj=TRUE)
# t2 = Sys.time()
# print(t2-t1)

corre.mnn <- out.mnn.total@assays@data$corrected # @assays[['corrected']]
mnn_processed<-list()
mnn_processed[['count']]=corre.mnn
mnn_processed[['group']]=annot$Sample_Origin[match(corre.mnn%>%colnames(),annot$Index)]
mnn_processed[['batch']]=annot$Patient[match(corre.mnn%>%colnames(),annot$Index)]
mnn_processed
