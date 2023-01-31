#RISC output as input processed.
run_QP<-function(res,cellinfo,former.meth='risc'){
  library(RISC)
  rownames(cellinfo)=cellinfo$Cell
  cellinfo$Group%<>%factor()
  cellinfo$Batch%<>%factor()
  
  obj0<-res
  obj0 = scUMAP(obj0, npc = 3,use = 'PLS')
  obj0 = scCluster(obj0, slot = "cell.umap", k = 3, method = 'louvain')
  
  cell.ctrl = rownames(obj0@coldata)[obj0@coldata$Group == sort(unique(obj0@coldata$Group))[1]]
  cell.sam = rownames(obj0@coldata)[obj0@coldata$Group == sort(unique(obj0@coldata$Group))[2]]
  res = scDEG(obj0, cell.ctrl = cell.ctrl, cell.sam = cell.sam, ncore=1, Padj=1, frac=0, logFC=0,
               min.cells = 0, method = 'QP')
  
  
  res_name<-paste0(ifelse(former.meth=='','',paste0(former.meth,'+')),'QP')
  save(res, cellinfo, file=paste0('./',res_name,'.rda'))
  return(res_name)
}