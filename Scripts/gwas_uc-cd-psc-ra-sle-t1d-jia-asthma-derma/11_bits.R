library(data.table)

loci.table <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/loci_table_names_nicola_genes.csv', data.table = F) #ta


loci.table$f1 <- rep(NA, nrow(loci.table))
loci.table$f2 <- rep(NA, nrow(loci.table))
loci.table$f3 <- rep(NA, nrow(loci.table))
loci.table$multi <- rep(NA, nrow(loci.table))

loci.table$final.locus[i]

unique(loci.table[loci.table$final.locus==loci.table$final.locus[i], ]$trait)

for(i in 1:nrow(loci.table)){
loci.table[loci.table$final.locus==loci.table$final.locus[i],]$multi <-  paste0(unique(loci.table[loci.table$final.locus==loci.table$final.locus[i], ]$trait), collapse = '-')
}

loci.table[!loci.table$multi %in% c('f1', 'f2', 'f3', 'f1-f3', 'f1-f2', 'f2-f3' ),c('trait','SNP','Chr',,'closest_gene_name')]