#calculate effectsize from modified function of MetaDE
#List of Log transformed count as input indexed by each batch input processed=list(batch1=batch1logcount, batch2=batch2logcount, ...)
run_ES<-function(processed,cellinfo,former.meth){
  library(MetaDE)
  library(DESeq2)
  rownames(cellinfo)=cellinfo$Cell
  cellinfo$Group%<>%factor()
  cellinfo$Batch%<>%factor()
  
  processed2<-lapply(processed,FUN=function(x){y<-list()
  y$count<-x
  y$group<-factor(cellinfo[colnames(x),'Group'])
  return(y)})
  
  res<-ind.cal.ES.fix(first_processed_pseudo, vs=processed$group%>%unique()%>%sort(decreasing = T))
  colnames(res$ES)=colnames(res$Var)=names(processed2)
  
  save(res, cellinfo, file=paste0('./',former.meth,'es.rda'))
}

#ind.cal.ES, cal.ES, get.ES functions from MetaDE(version 1.0.5) are modified for this analysis
ind.cal.ES.fix<-function(x, vs){
  K <- length(x)
  res <- get.ES.fix(x, vs)
  if (is.null(names(x))) {
    colnames(res$ES) <- colnames(res$Var) <- paste("dataset", 
                                                   1:K, sep = "")
  }else {
    colnames(res$ES) <- colnames(res$Var) <- names(x)
  }
  result <- list(ES = res$ES, Var = res$Var)
  return(result)  
}
get.ES.fix<-function(x, vs){
  K <- length(x)
  ES.m <- Var.m <- N <- n <- NULL
  for (k in 1:K) {
    y <- x[[k]][[1]]
    l <- x[[k]][[2]]
    temp <- cal.ES.fix(y, l, vs=vs)
    ES.m <- cbind(ES.m, temp[, "dprime"])
    Var.m <- cbind(Var.m, temp[, "vardprime"])
    N <- c(N, length(l))
    n <- c(n, table(l))
  }
  rownames(ES.m) <- rownames(y)
  rownames(Var.m) <- rownames(y)
  colnames(ES.m) <- names(x)
  colnames(Var.m) <- names(x)
  res <- list(ES = ES.m, Var = Var.m)
  return(res)
}


cal.ES.fix<-function (y, l, vs){
  library(gmp)
  l <- unclass(factor(l))
  n <- table(factor(l))
  ind <- diag(rep(1, length(n)))[l, ]
  ym <- y %*% ind %*% diag(1/n)
  ntilde <- 1/sum(1/n)
  m = sum(n) - 2
  # cm = gamma(m/2)/(sqrt(m/2) * gamma((m - 1)/2))
  cm = as.numeric(factorialZ((m/2)-1)/(sqrt(m/2) * factorialZ(((m - 1)/2)-1)))
  ind <- diag(rep(1, length(n)))[l, ]
  ym <- y %*% ind %*% diag(1/n)
  
  s <- sqrt((1/(sum(n)- 2)) * ((((y-ym[,l])^2)%*%ind)%*%diag(1/(n - 1))%*%(n - 1)))
  
  # s <- sqrt((1/(sum(n) - 2) * ((y^2 %*% ind) %*% diag(1/(n - 1)) - ym^2 %*% diag(n/(n - 1))) %*% (n - 1)))
  # d <- (ym[, 2] - ym[, 1])/s
  d <- (ym[, match(vs,levels(l))[1]] - ym[, match(vs,levels(l))[2]])/s
  
  d[which(s==0)]=Inf
  d[which((ym[, match(vs,levels(l))[1]] - ym[, match(vs,levels(l))[2]])==0)]=0
  
  dprime = d - 3 * d/(4 * (sum(n) - 2) - 1)
  terme1 = 1/ntilde
  vard = terme1 + d^2 * (terme1 * ntilde - 1/cm^2)
  vardprime = sum(1/n) + dprime^2/(2 * sum(n))
  result = cbind(dprime, vardprime)
  colnames(result) = c("dprime", "vardprime")
  rownames(result) <- rownames(y)
  return(result)
}