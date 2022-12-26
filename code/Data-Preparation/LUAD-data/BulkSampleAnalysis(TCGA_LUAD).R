`%ni%`<-Negate(`%in%`)
data.type='htseq_count'
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl"
)
# setwd('~/)
for(tumor.type in c('LUAD')){
  count<-read.delim(paste0('TCGA-',tumor.type,'.htseq_counts.tsv.gz'))
  count[,2:ncol(count)]<-(2^count[,2:ncol(count)])-1
  samples<-colnames(count)[2:ncol(count)]
  probmap<-read.table('gencode.v22.annotation.gene.probeMap',header = T)
  count$GeneName<-probmap$gene[match(count$Ensembl_ID,probmap$id)]
  # count$Ensembl_ID.original<-count$Ensembl_ID
  # count$Ensembl_ID%<>%gsub(pattern = '[.].*',replacement = '')
  
  mrna_attributes <- getBM(attributes=c("external_gene_name",
                                        "ensembl_gene_id",
                                        "gene_biotype"),
                           filters = c("external_gene_name"),
                           values = count$GeneName,
                           mart = mart)
  mrna_attributes<-mrna_attributes[which(mrna_attributes$gene_biotype=='protein_coding'),]  
  count<-count[which(count$GeneName%in%mrna_attributes$external_gene_name),]
  c1<-count[which(count$GeneName%ni%names(table(count$GeneName)[table(count$GeneName)==1])),]
  if(any(table(GeneName)>=2)){
    for(g in names(table(count$GeneName))[table(count$GeneName)>1]){
      c2.t<-count%>%dplyr::filter(GeneName==g)
      c1<-rbind(c1,c2.t[which.max(apply(c2.t[,2:(ncol(c2.t)-1)],1,mean)),])
    }
  }
  count<-c1
  rownames(count)<-count$GeneName
  count<-count[,setdiff(colnames(count),c('Ensembl_ID','GeneName'))]
  write.table(count, file=paste0('TCGA-',tumor.type,'.rnaseq.mapped',ifelse(is.null(data.type),'',paste0('.',data.type)),'.txt'))
}

for(tumor.type in c('LUAD')){
  rnaseq2<-read.table(file=paste0('TCGA-',tumor.type,'.rnaseq.mapped',ifelse(is.null(data.type),'',paste0('.',data.type)),'.txt'))
  if(tumor.type == "LUAD"){
    pset<-read.delim(paste0('TCGA-',tumor.type,'.GDC_phenotype.tsv.gz'))
    pset[["submitter_id.samples"]]%<>%gsub(pattern='[-]',replacement='.')
    pset%<>%dplyr::filter(submitter_id.samples %in% colnames(rnaseq2))
    
    rownames(pset)=pset[['submitter_id.samples']]
    
    idx<-match(c("gender.demographic",'age_at_initial_pathologic_diagnosis','sample_type.samples','tobacco_smoking_history'),colnames(pset))
    metadata <- data.frame(pset[,idx[!is.na(idx)]],
                           row.names = rownames(pset))
    
    colnames(metadata) <- c('gender', 'age', 'histology',
                            'smoking')[!is.na(idx)]
    metadata$smoking[metadata$smoking%in%c(1)]='Non-smoker'
    metadata$smoking[metadata$smoking%in%c(2,3,4,5)]='Smoker'
    metadata$histology[metadata$histology!="Solid Tissue Normal"]='Tumor'
    metadata$histology[metadata$histology=="Solid Tissue Normal"]='Normal'
    metadata <- metadata[!apply(metadata, 1, function(x) any( is.na(x) )),]
    
    rnaseq2<-rnaseq2[,which(colnames(rnaseq2)%in%rownames(metadata))]
    metadata<-metadata[match(colnames(rnaseq2),rownames(metadata)),]
    
    
  }
  for(meth in c('DESeq2','edgeR','limmatrend','limmavoom')){
    for(cov in c(T)){
      result.table<-run_TCGA.deg(dat=rnaseq2,metadata=metadata, meth=meth, min.gene.filter=5, cov=cov)
      write.table(result.table,file=paste0('TCGA-',tumor.type,'.rnaseq.processed',ifelse(is.null(data.type),'',paste0('.',data.type)),'.',meth,ifelse(cov,'_Cov',''),'.txt'))
    }
  }
}





run_TCGA.deg<-function(dat, metadata, meth='DESeq2', min.gene.filter=5, cov=T){
  library(DESeq2)
  library(edgeR)
  library(limma)
  count_df<-dat[which(apply(dat,1,FUN = function(x){as.numeric(x)%>%mean()})>=min.gene.filter),]
  
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
    count_df <- round(count_df, 0) + 1
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
    # count_df<-dat
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
  }else if(meth=='limmavoom'){
    # count_df<-dat
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
  }else if(meth=='limmatrend'){
    # count_df<-dat
    count_df[is.na(count_df)] = 0.
    
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