ref_df<-read.delim('C0152013_disease_gda_summary.tsv')
ref_df%<>%dplyr::filter(Score_gda>=0.3)
ref_df%<>%dplyr::arrange(desc(Score_gda))
ref_genes.1<-ref_df$Gene
ref_weight.1<-ref_df$Score_gda
ref_df<-read.delim('CTD_D000077192_genes_20210713220641.tsv')
ref_df%<>%dplyr::filter(Direct.Evidence=="marker/mechanism")
ref_genes.2<-ref_df$Gene.Symbol
both.in<-intersect(ref_genes.1,ref_genes.2)
Known_disease_genes<-c(ref_genes.1,setdiff(ref_genes.2,ref_genes.1))
Known_disease_genes_weight=c(ref_weight.1, rep(median(ref_weight.1[ref_genes.1%in%both.in]), setdiff(ref_genes.2,ref_genes.1)%>%length()))
