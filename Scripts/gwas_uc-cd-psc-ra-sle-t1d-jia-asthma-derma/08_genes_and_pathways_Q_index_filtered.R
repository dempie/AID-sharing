library(data.table)
library(ChIPpeakAnno)
library(stringr)
library(gprofiler2)
library(GenomicRanges)
library(EnsDb.Hsapiens.v75)
library(biomaRt)
library(ComplexHeatmap)
library(clusterProfiler)
library(enrichplot)
library(GOSemSim)
library(gtools)
library(RColorBrewer)


#------- add the Q heterogeneity index to the factro loci list -----------------

factor_loci <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/factor_loci_moloc_closest_gene_info.txt', data.table = F) 

factor_loci[,'Q_index_pvalue'] <- rep(0, nrow(factor_loci))
gwas_list <- list()

for(i in 1:3){
  #load the gwas 
  tt <- c('f1', 'f2', 'f3')[i]
  gwas_list[[tt]] <- fread(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/factor',i,'_gwas_final_withQindex.txt'), data.table = F)
  #add the qpval   
  for(k in 1:nrow(factor_loci[factor_loci$trait==tt, ])){
    to_take <- factor_loci[factor_loci$trait==tt,]$SNP[k]
    factor_loci[factor_loci$trait==tt & factor_loci$SNP==to_take,]$Q_index_pvalue <- gwas_list[[tt]][gwas_list[[tt]]$SNP==to_take, ]$Q_chisq_pval
  }
  
}


fwrite(factor_loci, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/factor_loci_with_qindex.csv')
factor_loci <-fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/factor_loci_with_qindex.csv', data.table = F) 

#make a list of genes for all the factors with ensemble_gened and gene symbol
mart <- useDataset("hsapiens_gene_ensembl", useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")) #select the database to convert and maake sure it is build 37 
f_list_genes <- list()
f_list_genes_Q <- list()
for(i in 1:3){
  tt <- c('f1', 'f2', 'f3')[i]
  f_list_genes[[tt]] <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values= factor_loci[factor_loci$trait==tt & factor_loci$Q_index_pvalue>5e-8, ]$closest_gene ,mart= mart)
  f_list_genes_Q[[tt]] <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values= factor_loci[factor_loci$trait==tt & factor_loci$Q_index_pvalue<5e-8, ]$closest_gene ,mart= mart)
}


#----- plot the upset plot of the Q heterogeneity index genes ------------------
q <- make_comb_mat(list(f1=f_list_genes$f1$ensembl_gene_id, f2=f_list_genes$f2$ensembl_gene_id, f3=f_list_genes$f3$ensembl_gene_id))
pdf(width = 10, height = 5, file = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/with_Q_index/intersect_genes_Qindex_no_significant.pdf')
UpSet(q, set_order = c("f1", "f2", "f3"), comb_order = order(comb_size(q), decreasing = T),
      comb_col = c(rcartocolor::carto_pal(12, 'Safe')[c(1,2,3)][comb_degree(q)]),
      top_annotation = upset_top_annotation(q, add_numbers = TRUE, height = unit(6, "cm")),
      right_annotation = upset_right_annotation(q, add_numbers = TRUE, width = unit(5,'cm') ),
      row_title = "Factor", 
      column_title = "Intersection of genes with no significant Q index ")

dev.off()


#
pdf(width = 10, height = 5, file = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/with_Q_index/intersect_genes_Qindex_significant.pdf')
u <- make_comb_mat(list(f1=f_list_genes_Q$f1$ensembl_gene_id, f2=f_list_genes_Q$f2$ensembl_gene_id, f3=f_list_genes_Q$f3$ensembl_gene_id))
UpSet(u, set_order = c("f1", "f2", "f3"), comb_order = order(  comb_size(u), decreasing =F),
      comb_col = c(rcartocolor::carto_pal(12, 'Safe')[c(1,2,3)][comb_degree(u)]),
      top_annotation = upset_top_annotation(u, add_numbers = TRUE, height = unit(6, "cm")),
      right_annotation = upset_right_annotation(u, add_numbers = TRUE, width = unit(5,'cm') ),
      row_title = "Factor",
      column_title = "Intersection of genes significant for Q index ")
dev.off()



















