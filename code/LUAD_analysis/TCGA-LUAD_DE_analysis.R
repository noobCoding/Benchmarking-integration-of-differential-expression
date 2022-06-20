source('/hdd2/SC lung/script_bk/Bulk RNA-seq_analysis.sources.R')
`%ni%`<-Negate(`%in%`)
data.type='htseq_count'
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl"
                # , host="grch37.ensembl.org"
)
base_dir='~/'
# surv_dir=paste0(base_dir,'data/for surv/')
# dir.create(surv_dir,showWarnings = F,recursive=T)

for(tumor.type in c('LUAD')){
  count<-read.delim(paste0(base_dir,'data/original/TCGA-',tumor.type,'.htseq_counts.tsv.gz'))
  count[,2:ncol(count)]<-(2^count[,2:ncol(count)])-1
  
  
  # count$Ensembl_ID.original<-count$Ensembl_ID
  count$Ensembl_ID%<>%gsub(pattern = '[.].*',replacement = '')
  mrna_attributes <- getBM(attributes=c("external_gene_name",
                                        "ensembl_gene_id",
                                        "gene_biotype"),
                           filters = c("ensembl_gene_id"),
                           values = count$Ensembl_ID,
                           mart = mart)
  mrna_attributes<-mrna_attributes[which(mrna_attributes$gene_biotype=='protein_coding'),]  
  count<-count[which(count$Ensembl_ID %in% mrna_attributes$ensembl_gene_id),]
  GeneName=mrna_attributes$external_gene_name[match(mrna_attributes$ensembl_gene_id,count$Ensembl_ID)]
  colnames(count)[1]="GeneName"
  count$GeneName = GeneName
  count<-count[which(count$GeneName%ni%names(table(GeneName)[table(GeneName)>=2])),]
  rownames(count)<-count$GeneName
  count<-count[,2:ncol(count)]
  write.table(count, file=paste0(base_dir,'data/preprocessed/TCGA-',tumor.type,'.rnaseq.processed',ifelse(is.null(data.type),'',paste0('.',data.type)),'.txt'))
}

for(tumor.type in c('LUAD')){
  rnaseq2<-read.table(file=paste0(base_dir,'data/preprocessed/TCGA-',tumor.type,'.rnaseq.processed',ifelse(is.null(data.type),'',paste0('.',data.type)),'.txt'))
  if(tumor.type == "LUAD"){
    pset<-read.delim(paste0(base_dir,'data/original/TCGA-',tumor.type,'.GDC_phenotype.tsv.gz'))
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
  for(meth in c('DESeq2','edgeR','limma_trend','voom.limma')){
    for(cov in c(T,F)){
      result.table<-run_TCGA.deg(dat=rnaseq2,metadata=metadata, meth=meth, min.gene.filter=5, cov=cov)
      
      write.table(result.table,file=paste0(base_dir,'data/analyzed/TCGA-',tumor.type,'.rnaseq.processed',ifelse(is.null(data.type),'',paste0('.',data.type)),'.',meth,ifelse(cov,'_Cov',''),'.txt'))
    }
  }
}
