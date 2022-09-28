library(data.table)
library(dplyr)

#----Supplementary table 2, table of conditional loci --------------------------

loci.table <-  fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/10_genes_and_pathways/loci_table_names_nicola_genes.csv', data.table = F) 
loci.table$final.locus=paste0(loci.table$Chr,":",loci.table$start,"-",loci.table$end,"_",loci.table$sub_locus)
head(loci.table)

#Add the column for independent signals from the same trait ended up colocalising in the same locus
loci.table$indip_coloc <- rep(NA, nrow(loci.table))

for(i in 1:3){
  tt <- c('f1', 'f2', 'f3')[i]
  dup <- loci.table[loci.table$trait== tt, ]$final.locus[duplicated(loci.table[loci.table$trait==tt, ]$final.locus)]
  loci.table[loci.table$trait== tt & loci.table$final.locus %in% dup, ]$indip_coloc <- 'multiple singal colocalising in the same locus'
  loci.table[loci.table$trait== tt & (!loci.table$final.locus %in% dup), ]$indip_coloc <- '-'
}


#change the factor names
loci.table$trait[loci.table$trait=='f1'] <- 'Fgut'
loci.table$trait[loci.table$trait=='f2'] <- 'Faid'
loci.table$trait[loci.table$trait=='f3'] <- 'Falrg'

#renmae the columns
loci <- loci.table %>% rename(Factor=trait, 'ENSEMBLEID_CLOSEST_PROTEIN_CODING_GENE'=closest_gene, 
                              'GENE_NAME_CLOSEST_PROTEIN_CODING_GENE'=closest_gene_name, 'LOCUS_NAME(CHR_START_END)'= final.locus, 'NOTE'=indip_coloc  ) %>% 
  select(-c(pan.locus, sub_locus) )
colnames(loci) <- toupper(colnames(loci))

#save
fwrite(loci,'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/tables/Supplementary_table_conditional_loci.csv', sep = ',', col.names = T, quote = F)

