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


# plot an upset plot of the q index genes
pdf(width = 10, height = 5, file = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/with_Q_index/intersect_genes_Qindex_significant.pdf')
u <- make_comb_mat(list(f1=f_list_genes_Q$f1$ensembl_gene_id, f2=f_list_genes_Q$f2$ensembl_gene_id, f3=f_list_genes_Q$f3$ensembl_gene_id))
UpSet(u, set_order = c("f1", "f2", "f3"), comb_order = order(  comb_size(u), decreasing =F),
      comb_col = c(rcartocolor::carto_pal(12, 'Safe')[c(1,2,3)][comb_degree(u)]),
      top_annotation = upset_top_annotation(u, add_numbers = TRUE, height = unit(6, "cm")),
      right_annotation = upset_right_annotation(u, add_numbers = TRUE, width = unit(5,'cm') ),
      row_title = "Factor",
      column_title = "Intersection of genes significant for Q index ")
dev.off()


#-------heatmap ----------------------------------------------------------------

gost_no_q <- gost(query = list(f1=f_list_genes$f1$ensembl_gene_id, f2=f_list_genes$f2$ensembl_gene_id, f3=f_list_genes$f3$ensembl_gene_id),
                           sources = c('KEGG'), significant = T, evcodes = T)
saveRDS(gost_no_q, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/with_Q_index/kegg_pathways_without_Q_index_significant.RDS')

he <- gost_no_q$result
genes <- unique(strsplit(paste0(he$intersection, collapse = ','), split = ',')[[1]])
genes_symbol <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values=  genes,mart=mart )

t_t <- data.frame(genes) 
for(i in 1:length(he$term_name)){
  t_t[, paste0(paste(strsplit(he$term_name[i], split=' ')[[1]], collapse = '_'),'_' , he$query[i])] <- as.numeric(genes %in% strsplit(he[he$term_name==he$term_name[i] & he$query==he$query[i], 'intersection'], split=',' )[[1]])
  
}

#add gene symbol and remove ENSEMBLE
rownames(t_t) <- genes_symbol$hgnc_symbol[match(t_t$genes, genes_symbol$ensembl_gene_id)]   
t_t$genes <- NULL

t_t['trait',] <- c(he$query) 
t_t['p_value',] <- c(he$p_value)
tt <- t(t_t)

#color selection
library(rcartocolor)
display_carto_all(colorblind_friendly = TRUE)

colori<- ComplexHeatmap:::default_col(x = tt)
colori[c('f1', 'f3')] <-  rcartocolor::carto_pal(12, 'Safe')[c(4,10)]#colorblind safe
colori[c('0','1')] <- c('#a6cee3','#1f78b4')
c('#a6cee3','#1f78b4','#b2df8a')
#column split
sep <- c(rep('a', 25), rep('b', 2))
names(sep)<- colnames(tt)

pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/with_Q_index/heatmap_gene_and_pathways_no_significant_Q_genes.pdf', height = 8, width = 16)
Heatmap(tt, column_title = "Genes and pathways no significant Q",
        rect_gp = gpar(col = "white", lwd = 0.25), 
        column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        left_annotation = rowAnnotation( '-log10(Padj)' = anno_barplot(   axis_param = list(direction = "reverse"),-log10(as.numeric(tt[,'p_value'])), width = unit(2, "cm"))),
        row_split = tt[,'trait'],
        column_split = sep,
        column_gap = unit(3, "mm"),
        row_gap = unit(5, "mm"), 
        row_names_gp = gpar(col =  rcartocolor::carto_pal(12, 'Safe')[c(4,10)], fontsize=10),
        column_labels = c(colnames(tt)[1:25], 'Factor', '') ,
        column_names_side = 'bottom',
        show_heatmap_legend = F,
        row_title = "Pathways",
        na_col = 'black', 
        col = colori, 
        heatmap_width = unit(28, 'cm'), heatmap_height = unit(15, 'cm')
)
dev.off()

#----- are Q index genes significant for some KEGG pathway?---------------------

gost_q <- gost(query = list(f1=f_list_genes_Q$f1$ensembl_gene_id, f2=f_list_genes_Q$f2$ensembl_gene_id, f3=f_list_genes_Q$f3$ensembl_gene_id),
                  sources = c('KEGG'), significant = F, evcodes = T)
saveRDS(gost_q, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/with_Q_index/kegg_pathways_only_Q_index_significant.RDS')
gost_q$result




