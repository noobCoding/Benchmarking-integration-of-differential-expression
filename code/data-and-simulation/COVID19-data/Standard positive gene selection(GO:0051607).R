library(readr)
dbf<-read_lines('GO_Biological_Process_2021.txt')
dbl<-sapply(1:length(dbf),FUN=function(x){
  rt<-str_split(dbf[x],pattern='\t')[[1]]
})
names(dbl)=sapply(dbl,FUN=function(x){x[1]})
dbl2<-sapply(names(dbl),FUN=function(x){
  setdiff(dbl[[x]],c(x,''))
})
Standard_Positive_gene_COVID19<-dbl2$`defense response to virus (GO:0051607)`