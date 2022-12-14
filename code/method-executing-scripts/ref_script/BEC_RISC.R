library(RISC)
#Modified InPlot function to provide number of eigen vector. 
n_eigen<-function (object = NULL, var.gene = NULL, Colors = NULL, nPC = 20, 
                   neighbor = 30, method = "louvain", algorithm = "kd_tree", 
                   ncore = 4, minPC = 11, Std.cut = 0.95, bin = 5){
  library(irlba)
  library(doParallel)
  library(igraph)
  set.seed(123)
  if (is.null(objects)) {
    stop("Please keep objects valid")
  }else {
    ncore = as.integer(ncore)
    doParallel::registerDoParallel(ncore)
    npc = as.integer(nPC)
    neighbor = as.integer(neighbor)
    method = as.character(method)
    algorithm = as.character(algorithm)
    nset = dim(summary(object))[1]
    minpc = as.integer(minPC)
    Std.cut = as.numeric(Std.cut)
    bin = as.integer(bin)
    if (minpc >= npc) {
      minpc = npc
    }else {
      minpc = minpc
    }
  }
  gene0 = pc.gene = var0 = PC0 = list()
  vari = PCi = PC1 = PC2 = k0 = ki = ka = kb = kc = Std.Num = Std0 = Cluster = Num = PCs = Group = Group1 = Group2 = mean0 = NULL
  i = j = 1L
  while (i <= nset) {
    gene0[[i]] = object[[i]]@rowdata$Symbol
    var0[[i]] = object[[i]]@vargene
    i = i + 1L
  }
  gene0 = Reduce(intersect, gene0)
  var0 = Reduce(intersect, var0)
  if (is.null(var.gene)) {
    var0 = intersect(gene0, var0)
  }else {
    var0 = intersect(gene0, var.gene)
  }
  PC0 = foreach(i = 1L:nset) %dopar% {
    vari = object[[i]]@assay$logcount[var0, ]
    PCi = irlba(vari, nv = npc, center = T)
    return(PCi)
  }
  PC1 = foreach(i = 1L:nset) %dopar% {
    ka = apply(PC0[[i]]$u, 2, FUN = function(x) {
      ks.test(x, "pnorm", mean = mean(x), sd = sd(x))$statistic
    })
    kb = glm(ka ~ c(1:npc), family = quasipoisson)$deviance
    kc = mean(ka)
    return(list(ka, kb, kc))
  }
  ka = unlist(lapply(PC1, `[[`, 2))
  ka = (1 - ka) * (1/max(1 - ka))
  kc = unlist(lapply(PC1, `[[`, 3))
  PC2 = foreach(i = 1L:nset) %dopar% {
    PCi = PC0[[i]]$d^2/sum(PC0[[i]]$d^2)
    PCi = cumsum(PCi)
    Std0 = c(1:npc)[PCi > Std.cut][1]
    return(list(PCi, Std0))
  }
  kb = sapply(1L:nset, FUN = function(x) {
    glm(PC2[[x]][[1]] ~ c(1L:npc), family = poisson)$coefficients[2]
  })
  kb = 1 - (1/npc - kb)/(1/npc)
  Std.Num = unlist(lapply(PC2, `[[`, 2))
  bin0 = as.integer(seq(minpc, npc, length.out = bin))
  if (method == "louvain") {
    clust0 = foreach(i = rep(1L:nset, length(bin0)), j = rep(bin0, 
                                                             each = nset)) %dopar% {
                                                               k0 = FNN::get.knn(PC0[[i]]$v[, 1L:j], k = neighbor, algorithm = algorithm)
                                                               ki = data.frame(NodStar = rep(1L:nrow(PC0[[i]]$v[, 
                                                                                                                1L:j]), neighbor), NodEnd = as.vector(k0$nn.index), 
                                                                               stringsAsFactors = FALSE)
                                                               ki = igraph::graph_from_data_frame(ki, directed = FALSE)
                                                               ki = simplify(ki)
                                                               ki = cluster_louvain(ki, weights = 1/(1 + as.vector(k0$nn.dist)))
                                                               k0 = length(unique(ki$membership))
                                                               return(as.integer(k0))
                                                             }
  }else {
    stop("A new method later")
  }
  PC3 = data.frame(PC1 = unlist(lapply(PC1, `[[`, 1)), PC2 = unlist(lapply(PC2, 
                                                                           `[[`, 1)), Group = rep(paste0("Set-", 1L:nset), each = npc), 
                   Num = rep(1L:npc, nset), stringsAsFactors = F)
  PC3$Group = factor(PC3$Group, levels = paste0("Set-", 1L:nset))
  PC3$mean0 = rep(round(kc, 2), each = npc)
  PC3$Group1 = factor(PC3$Group, labels = paste0(levels(PC3$Group), 
                                                 ": ", round(ka, 2)))
  PC3$Group2 = factor(PC3$Group, labels = paste0(levels(PC3$Group), 
                                                 ": ", round(as.numeric(kb), 2)))
  PC3.use<-PC3%>%dplyr::filter(PC2<Std.cut)
  #If number of eigen vectors are too small, allow maximum value of recommended Std.cut range.
  if(max(PC3.use$Num)<5){
    PC3.use<-PC3%>%dplyr::filter(PC2<c(0.9,0.95,0.98)[min(which(sapply(c(0.9,0.95,0.98),FUN=function(x){return(PC3%>%dplyr::filter(PC2<x)%>%dplyr::select(Num)%>%max())})>=5))])
  }
  #Use the number of PC below std.cut as eigenvector to use.
  return(PC3.use$Num%>%max())
}
process0 <- function(obj0,large=F){
  # Filter cells and genes
  
  obj0 = scFilter(obj0, min.UMI = 0, max.UMI = Inf, min.gene = 0, min.cell = 0, mitochon
                  = 0, gene.ratio = 0,is.filter = F)
  # Normalize the raw counts
  obj0 = scNormalize(obj0, ncore = 4,large=large)
  # Find highly variable genes
  obj0 = scDisperse(obj0)
  # print(length(obj0@vargene))
  return(obj0)
}



#raw count as input count=rawcount
run_RISC<-function(count,cellinfo,P_ref.ind=NULL){
  rownames(cellinfo)=cellinfo$Cell
  cellinfo<-cellinfo[colnames(count),]
  cellinfo$Group%<>%factor()
  cellinfo$Batch%<>%factor()
  mat0 = as(count, 'dgCMatrix')
  large=(ncol(mat0)>50000)
  coldata0 = cellinfo[,c('Cell','Batch','Group')]
  coldata0%<>%dplyr::filter(Cell %in% colnames(mat0))
  data0.origin<-list()
  for(i in sort(unique(coldata0$Batch))){
    cell1 = coldata0$Cell[coldata0$Batch == i]
    coldata1 = coldata0[coldata0$Cell %in% cell1,]
    rownames(coldata1) = coldata1$Cell
    mat1 = mat0[,rownames(coldata1)]
    rowdata0 = data.frame(Symbol = rownames(mat1), row.names = rownames(mat1))
    ## Create RISC object
    # make sure the row-names of mat1 equal to row-names of rowdata0,
    # and the col-names of mat1 equal to row-names of coldata1
    dat1 = readsc(mat1, coldata1, rowdata0, is.filter = F)
    data0.origin[[i]]=process0(dat1,large)
  }
  
  var0 = Reduce(
    intersect, lapply(data0.origin,FUN=function(x){x@rowdata$Symbol})
  )
  data0 = data0.origin
  #"Std.cut" can be 0.85~0.9 for small size datasets (like total cells < 10,000)
  # and be 0.90~0.98 for large size datasets (total cells > 10,000).
  if(ncol(mat0)<10000){
    Std.cut<-((log10(ncol(mat0))-2)*0.05/2)+0.85
  }else if(ncol(mat0)<100000){
    Std.cut<-((log10(ncol(mat0))-4)*0.08/2)+0.9
  }else if(ncol(mat0)<1000000){
    Std.cut<-((log10(ncol(mat0))-5)*0.08/2)+0.94
  }
  if(is.null(P_ref.ind)){
    #Select the reference dataset from InPlot function.
    pdf('InPlot_output_RISC.pdf')
    InPlot(data0, var.gene = var0, Std.cut = Std.cut, ncore = 4, minPC = 16, nPC = 40)
    dev.off()
    return(NULL)
  }
  #If reference dataset is already set as input, we can do the further integration
  #Select the number of eigen_vector for further integration
  num_eigen=n_eigen(data0, var.gene = var0, Std.cut = Std.cut, ncore = 4, minPC = 16, nPC = 40)
  data0.integrate=data0.origin
  P_ref=names(data0.origin)[P_ref.ind]
  
  #Take reference data to the first place of list and name it 0Batch
  dat_p_ref=data0.integrate[[P_ref]]
  data0.integrate[[P_ref]]=data0.integrate[[1]]
  data0.integrate[[1]]=dat_p_ref
  n_dat_integrate<-names(data0.integrate)
  n_dat_integrate[n_dat_integrate==P_ref]=n_dat_integrate[1]
  n_dat_integrate[1]='0Batch'
  names(data0.integrate)=n_dat_integrate
  
  data0 = scMultiIntegrate(
    objects = data0.integrate, eigens = num_eigen, add.Id = NULL, var.gene = var0,
    align = 'OLS', npc = 50, adjust = TRUE, ncore = 4
  )
  res<-data0
  
  processed <- do.call(cbind, data0@assay$logcount)
  processed%<>%as.matrix()
  colnames(processed)<-gsub(pattern='^Set._',replace='',colnames(processed))
  
  res_name<-paste0(ifelse(former.meth=='','',paste0(former.meth,'+')),'risc')
  save(res, processed, cellinfo, file=paste0('./',res_name,'.rda'))
  return(res_name)
}
