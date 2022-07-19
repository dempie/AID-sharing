library(data.table)
library(ComplexHeatmap)
library(RColorBrewer)
library(qgraph)

#------ load the dataset -------------------------------------------------------
loci.table <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/loci_table_names_nicola_genes.csv', data.table = F)


#----- make an upset plot-------------------------------------------------------
q <- make_comb_mat(list(f1=loci.table[loci.table$trait=='f1',]$final.locus, f2=loci.table[loci.table$trait=='f2',]$final.locus, f3=loci.table[loci.table$trait=='f3',]$final.locus), mode = 'distinct')

pdf(width = 10, height = 4, file = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure_coloc/upset_plot_final_loci_nicola.pdf')
UpSet(q, set_order = c("f1", "f2", "f3"), 
      comb_order = c(5,6,7,2,3,4,1),
      comb_col = c(brewer.pal(12, 'Paired')[c(12,12,12,12,2,4,6)]),
      top_annotation = upset_top_annotation(q, add_numbers = TRUE, height = unit(6, "cm")),
      right_annotation = upset_right_annotation(q, add_numbers = TRUE, width = unit(5,'cm'),gp = gpar(fill = brewer.pal(5, 'Greys')[1])  ),
      row_title = "Factor", 
      column_title = "Intersection of all loci")
dev.off()


head(loci.table)
coloc_tab <-  fread('Nicola/loci_definitions/colocalization.table.H4.tsv', data.table = F)

topl<- loci.table[loci.table$pan.locus==23, c('trait', 'SNP')]


topl$trait <- toupper(topl$trait)

pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure_coloc/qgraph_locus_23_factors.pdf', width = 16, height = 16)
qgraph(topl[, c(2,1)], directed=F, layout='spring', 
       edge.labels=F, 
       groups=list(c(1,7), c(2,3,8), c(4,5,6,9)), 
       color=brewer.pal(5, 'Paired')[c(1,3,5)],
       edge.color=brewer.pal(5, 'Paired')[c(1,3,3,5,5,5)])
dev.off()




















