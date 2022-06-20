library(DESeq2)
library(edgeR)
library(limma)


run_TCGA.deg<-function(dat, metadata, meth='deseq2', min.gene.filter=5, cov=T){
  count_df<-dat[which(apply(dat,1,FUN = function(x){as.numeric(x)%>%mean()})>=min.gene.filter),]
  # rowN=rownames(count_df)
  # count_df<-count_df%>%sapply(as.numeric)%>%as.data.frame()
  # rownames(count_df)<-rowN
  library(DESeq2)
  library(edgeR)
  library(limma)
  
  cellinfo<-metadata
  for(i in colnames(cellinfo)){
    if(i=='age'){
      cellinfo[[i]]%<>%as.numeric()
    }else{
      cellinfo[[i]]%<>%factor()
    }
    
  }
  if(cov==T){
    design<-model.matrix(formula(paste(c("~ histology", setdiff(colnames(cellinfo),c('histology'))), collapse = '+')), data=cellinfo)
  }else{
    design<-model.matrix(formula(paste0("~ histology")), data=cellinfo)
  }
  
  if(meth=='DESeq2'){
    count_df <- round(dat, 0) + 1
    if(cov){
      dds <- DESeqDataSetFromMatrix(countData = count_df, colData = cellinfo, design = formula(paste(c("~ histology", setdiff(colnames(cellinfo),c('histology'))), collapse = '+')))
      dds <- DESeq2::DESeq(dds) #, fitType ='mean')
      res.origin<-results(dds, name="histology_Tumor_vs_Normal", cooksCutoff = F, independentFiltering=F)
      res <- lfcShrink(dds, coef =2 , res=res.origin, type="apeglm", lfcThreshold=0)
      # res <- lfcShrink(dds, coef =2 , type="apeglm", lfcThreshold=0)
      result.table <- data.frame('pvalue' = res$pvalue, 'adjpvalue' = res$padj, 'logFC' = res$log2FoldChange)
      
      rownames(result.table) <- rownames(dds)
    }else{
      dds <- DESeqDataSetFromMatrix(countData = count_df, colData = cellinfo, design = formula(paste0("~ histology")))
      dds <- DESeq2::DESeq(dds) #, fitType ='mean')
      res.origin<-results(dds, name="histology_Tumor_vs_Normal", cooksCutoff = F, independentFiltering=F)
      res <- lfcShrink(dds, coef =2 , res=res.origin, type="apeglm", lfcThreshold=0)
      # res <- lfcShrink(dds, coef=2, type="apeglm", lfcThreshold=0)
      result.table <- data.frame('pvalue' = res$pvalue, 'adjpvalue' = res$padj, 'logFC' = res$log2FoldChange)
      rownames(result.table) <- rownames(dds)
    }
  }else if(meth=='edgeR'){
    library(Seurat)
    library(edgeR)
    count_df<-dat
    y <- DGEList(counts=count_df, group=cellinfo[['histology']])
    y <- calcNormFactors(y)
    if(cov){
      y <- estimateDisp(y, design, robust=TRUE)
      fit <- glmQLFit(y, design, robust=TRUE, prior.df = 0)
      qlf <- glmQLFTest(fit, coef=2)
      FDR<-p.adjust(qlf$table$PValue,method = "BH")
      qlf$table$FDR <- FDR
      result.table <- data.frame('pvalue' = qlf$table$PValue, 'adjpvalue' = qlf$table$FDR, 'logFC' = qlf$table$logFC)
      rownames(result.table) <- rownames(qlf)
    }else{
      y <- estimateDisp(y, design, robust=TRUE)
      fit <- glmQLFit(y, design, robust=TRUE, prior.df = 0)
      qlf <- glmQLFTest(fit, coef=2)
      FDR<-p.adjust(qlf$table$PValue,method = "BH")
      qlf$table$FDR <- FDR
      result.table <- data.frame('pvalue' = qlf$table$PValue, 'adjpvalue' = qlf$table$FDR, 'logFC' = qlf$table$logFC)
      rownames(result.table) <- rownames(qlf)
    }
  }else if(meth=='voom.limma'){
    count_df<-dat
    nf <- edgeR::calcNormFactors(count_df, method = 'TMM')
    voom.data <- limma::voom(count_df, design = design, lib.size = colSums(count_df) * nf)
    voom.data$genes <- rownames(count_df)
    voom.fitlimma <- limma::lmFit(voom.data, design = design)
    voom.fitbayes <- limma::eBayes(voom.fitlimma)
    voom.pvalues <- voom.fitbayes$p.value[, 2]
    voom.adjpvalues <- p.adjust(voom.pvalues, method = 'BH')
    voom.logFC <- voom.fitbayes$coefficients[, 2]
    voom.score <- 1 - voom.pvalues
    result.table <- data.frame('pvalue' = voom.pvalues, 'adjpvalue' = voom.adjpvalues, 'logFC' = voom.logFC, 'score' = voom.score)
    rownames(result.table) <- rownames(count_df)
    return(result.table)
  }else if(meth=='limma_trend'){
    count_df<-dat
    count_df[is.na(count_df)] = 0.
    # count_df <- count_df + abs(min(count_df))
    
    ## Convert to an edgeR object
    dgeObj <- DGEList(count_df)
    ## Perform TMM normalisation
    dgeObj <- calcNormFactors(dgeObj)
    logCPM <- cpm(dgeObj, log=TRUE, prior.count=3)
    if(cov){
      lmfit <- lmFit(logCPM, design)
      lmfit <- eBayes(lmfit, trend=TRUE, robust = TRUE)
      res <- topTable(lmfit, n = Inf, adjust.method = "BH", coef = 2)
    }else{
      lmfit <- lmFit(logCPM, design)
      lmfit <- eBayes(lmfit, trend=TRUE, robust = TRUE)
      res <- topTable(lmfit, n = Inf, adjust.method = "BH", coef = 2)
    }
    result.table <- data.frame('pvalue' = res$P.Value, 
                               'adjpvalue' = res$adj.P.Val, 
                               'logFC' = res$logFC,
                               row.names = rownames(res))
  }
  return(result.table)
}
run_limma_cov.microarray<-function(dat, metadata, min.gene.filter=5, log.transformed=T){
  count_df<-dat[which(apply(dat,1,FUN = function(x){as.numeric(x)%>%mean()})>=ifelse(log.transformed,log2(min.gene.filter),min.gene.filter)),]
  
  library(edgeR)
  library(limma)
  
  cellinfo<-metadata
  for(i in colnames(cellinfo)){
    if(i=='age'){
      cellinfo[[i]]%<>%as.numeric()
    }else{
      cellinfo[[i]]%<>%factor()
    }
  }
  
  expset<-ExpressionSet(assayData=count_df%>%as.matrix())
  design<-model.matrix(formula(paste(c("~ histology", setdiff(colnames(cellinfo),c('histology'))), collapse = '+')), data=cellinfo)
  fit <- lmFit(expset, design)
  fitbayes <- limma::eBayes(fit)
  pvalues <- fitbayes$p.value[, 2]
  adjpvalues <- p.adjust(pvalues, method = 'BH')
  logFC <- fitbayes$coefficients[, 2]
  score <- 1 - pvalues
  result.table <- data.frame('pvalue' = pvalues, 'adjpvalue' = adjpvalues, 'logFC' = logFC, 'score' = score)
  rownames(result.table) <- rownames(count_df)
  return(result.table)
}
do.mapping.microarray<-function(geo_acc,eset,fset){
  mart <- useMart('ENSEMBL_MART_ENSEMBL',host = 'useast.ensembl.org')
  mart <- useDataset("hsapiens_gene_ensembl", mart)
  if(geo_acc=='GSE32863'){
    platform='illumina_humanwg_6_v3'
    filt=attr='external_gene_name'
    
  }else if(geo_acc=='GSE10072'){
    platform='affy_hg_u133a'
    attr=filt='affy_hg_u133a'
  }else if(geo_acc=='GSE43458'){
    platform='affy_hugene_1_0_st_v1'
    filt=attr='affy_hugene_1_0_st_v1'
  }else if(geo_acc=="GSE31210"){
    platform='affy_hg_u133_plus_2'
    attr=filt='affy_hg_u133_plus_2'
  }else{
    platform='affy_hg_u133_plus_2'
    filt=attr='external_gene_name'
  }
  mrna_attributes <- getBM(mart = mart,
                           # attributes = c(platform,
                           attributes = c(platform,
                                          'ensembl_gene_id',
                                          'gene_biotype',
                                          'external_gene_name'),
                           filter = filt,
                           values = rownames(eset),
                           uniqueRows = TRUE)
  mrna_attributes%<>%dplyr::filter(gene_biotype=='protein_coding')
  eset<-eset[rownames(eset)%in%mrna_attributes[[attr]],]
  rownames(eset)<-mrna_attributes$external_gene_name[match(rownames(eset),mrna_attributes[[attr]])]
  eset<-eset[str_detect(rownames(eset),pattern='[///]', negate = T)%>%which(),]
  {
    eset.1<-eset[!duplicated(rownames(eset)),]
    eset.2<-eset[duplicated(rownames(eset)),]
    probavg<-apply(eset.2, 1, mean)%>%as.vector()
    
    maxprob<-sapply(unique(rownames(eset.2)),FUN=function(x){
      m1<-probavg[x==rownames(eset.2)]
      return(which(x==rownames(eset.2))[which(m1==max(m1))])
    })%>%unlist()%>%as.vector()
    res_dt<-rbind(eset.1,eset.2[maxprob,])
    
    res_dt%<>%as.data.frame()
    res_dt<-res_dt[setdiff(unique(rownames(res_dt))%>%sort(),c("")),]
  }
  return(res_dt)
}

scal <- function(x,y){
  mean_n <- rowMeans(y)  # mean of normal
  sd_n <- apply(y,1,sd)  # SD of normal
  # z score as (value - mean normal)/SD normal
  res <- matrix(nrow=nrow(x), ncol=ncol(x))
  colnames(res) <- colnames(x)
  rownames(res) <- rownames(x)
  for(i in 1:dim(x)[1]){
    for(j in 1:dim(x)[2]){
      res[i,j] <- (x[i,j]-mean_n[i])/sd_n[i]
    }
  }
  return(res)
}
do.mapping.survival<-function(eset,fset){
  eset<-eset[!is.na(fset$`Gene Symbol`),]
  if(any(grep('^AFFX', rownames(eset)))){
    eset <- eset[-grep('^AFFX', rownames(eset)),]
  }
  rownames(eset)<-fset$`Gene Symbol`[match(rownames(eset),fset$ID)]
  eset<-eset[str_detect(rownames(eset),pattern='[///]', negate = T)%>%which(),]
  mart <- useMart('ENSEMBL_MART_ENSEMBL',host = 'asia.ensembl.org')
  mart <- useDataset("hsapiens_gene_ensembl", mart)
  mrna_attributes <- getBM(mart = mart,
                           attributes = c('affy_hg_u133_plus_2',
                                          'ensembl_gene_id',
                                          'gene_biotype',
                                          'external_gene_name'),
                           filter = 'external_gene_name',
                           values = rownames(eset),
                           uniqueRows = TRUE)
  mrna_attributes%<>%dplyr::filter(gene_biotype=='protein_coding')
  eset<-eset[rownames(eset)%in%mrna_attributes$external_gene_name,]
  
  
  a1<-names(table(mrna_attributes$external_gene_name)[table(mrna_attributes$external_gene_name)>=2])
  a2<-names(table(mrna_attributes$external_gene_name)[table(mrna_attributes$external_gene_name)==1])
  
  
  {
    eset.sub1<-eset[rownames(eset)%in%a1,]
    eset.sub2<-eset[rownames(eset)%in%a2,]
    # match(rownames(eset.sub1),mrna_attributes$affy_hg_u133_plus_2)
    sum.dup<-apply(eset.sub1,1,sum)
    eset.sub1.max<-sapply(a1,FUN=function(x){
      which(names(sum.dup)==x)[which(sum.dup[which(names(sum.dup)==x)]==max(sum.dup[which(names(sum.dup)==x)]))]
      # names(sum.dup[mrna_attributes$affy_hg_u133_plus_2[mrna_attributes$external_gene_name==x]]%>%sort(decreasing = T))[1]
    })%>%as.vector()
    eset.sub1<-eset.sub1[eset.sub1.max,]
    eset.renew<-rbind(eset.sub1,eset.sub2)
    # rownames(eset.renew)<-mrna_attributes$external_gene_name[match(rownames(eset.renew),mrna_attributes$affy_hg_u133_plus_2)]
    eset.renew%<>%as.data.frame()
    eset.renew%<>%dplyr::arrange(rownames(eset.renew))
  }
  
  return(eset.renew)
}
get.wFisher.survival<-function(p,weight,eff=NULL,is.onetail=T){
  library(metapro)
  if(!is.null(eff)){
    eff.sign=log2(eff)%>%as.data.frame()
    eff.sign2<- lapply(eff.sign, function(x) ifelse(x>0,1,-1))%>%as.data.frame(row.names = rownames(HR_table))
    message('eff exist, do one tail false')
  }else{
    eff.sign2=NULL
    message('eff null, do one tail')
  }
  
  k <- ncol(p)
  pval <- stat <- rep(NA, nrow(p))
  rnum <- 1:nrow(p)
  res.proto<-unlist(sapply(X=rnum, FUN=function(x){metapro::wFisher(p[x,]%>%as.numeric(),weight=weight,is.onetail = is.onetail,eff.sign = ifelse(is.null(eff),NULL,eff.sign2[x,]))}))
  
  if(is.onetail==T){
    pval=res.proto
    sign=NULL
  }else{
    pval=res.proto[c(1:(length(res.proto)/2))*2-1]%>%as.numeric()
    sign=res.proto[c(1:(length(res.proto)/2))*2]
  }
  qval <- p.adjust(pval, method = "BH")
  res <- list(stat = pval, pval = pval, FDR = qval, sign=sign)
  
  names(res$stat) <- names(res$pval) <- names(res$FDR) <- rownames(p)
  return(res)
}