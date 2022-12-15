library(stringr)
library(data.table)
library(dplyr)


#mr table
mr_list <- list.files('/project/aid_sharing/sc_eqtl/')[str_ends(list.files('/project/aid_sharing/sc_eqtl/'), 'MR.tsv')]

setwd('/project/aid_sharing/sc_eqtl/')

mr_l <- list()
mr_ok <- list()

for(i in mr_list){
  
  mr_l[[i]] <- try(fread(i))
  
  if(!is(mr_l[[i]], 'try-error')){
    mr_ok[[i]] <- mr_l[[i]]
  }
  
}


mr_res <- Reduce(rbind, mr_ok)
mr_res <- mr_res[mr_res$trait %in% c('f1', 'f2', 'f3'),]

#add the locus name to the table
loci <- fread('/project/aid_sharing/AID_sharing/outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/09_loci_definitions_Nicola/loci_definitions/final_locus_table.tsv', data.table = F)
loci$final.locus=paste0(loci$Chr,":",loci$start,"-",loci$end,"_",loci$sub_locus)
loci$tmp <- paste0(loci$pan.locus, '_',loci$sub_locus)
head(loci)

mr_res$locus_name <- rep('', nrow(mr_res))
for(i in 1:nrow(mr_res)){
  mr_res[i, ]$locus_name <- loci[loci$tmp==mr_res[i, ]$locus,]$final.locus[1]
  
}

mr_final <- mr_res %>% rename( 'locus_name_(chr:start_end)'=locus_name, 'cell_type'=cell.type) %>% select(-c('locus'))

mr_final$trait[mr_final$trait=='f1'] <- 'Fgut'
mr_final$trait[mr_final$trait=='f2'] <- 'Faid'
mr_final$trait[mr_final$trait=='f3'] <- 'Falrg'
mr_final$trait[mr_final$trait=='asthma'] <- 'eczema'

names(mr_final) <- toupper(names(mr_final))

#remove duplicated rows for multiple locus associated with the same eQTL
mr_final <- mr_final[!is.na(mr_final$B),]

#save
fwrite(mr_final,'/project/aid_sharing/AID_sharing/outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/tables/Supplementary_table_MR_results.csv', sep = ',', col.names = T, quote = F)



#--------eQTL table-------------------------------------------------------------
qtl_names <- list.files()[str_ends(list.files(), 'coloc.tsv')]
list.files()

qtl_l <- list()
qtl_ok <- list()

for(i in qtl_names){
  
  qtl_l[[i]] <- try(fread(i))
  
  if(c(!is(qtl_l[[i]], 'try-error') & !ncol(qtl_l[[i]])==1)){
    qtl_ok[[i]] <- qtl_l[[i]]
  }
  
}


qtl_res<- Reduce(rbind, qtl_ok)
qtl_factors <- qtl_res[qtl_res$t1 %in% c('f1', 'f2', 'f3'),]


qtl_final <- qtl_factors %>% rename('hit_trait_1'=hit1, 'hit_trait_2'=hit2, 'cell_type'=cell_type, 'trait_1'= t1, 'trait_2'=t2  )

qtl_final$trait_1[qtl_final$trait_1=='f1'] <- 'Fgut'
qtl_final$trait_1[qtl_final$trait_1=='f2'] <- 'Faid'
qtl_final$trait_1[qtl_final$trait_1=='f3'] <- 'Falrg'

colnames(qtl_final) <- toupper(colnames(qtl_final))

head(qtl_final)

#save
fwrite(qtl_final,'/project/aid_sharing/AID_sharing/outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/tables/Supplementary_table_eQTL_results.csv', sep = ',', col.names = T, quote = F)


#---- colocalization table -----------------------------------------------------

posterior <- fread('/project/aid_sharing/AID_sharing/outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/09_loci_definitions_Nicola/loci_definitions/colocalization.table.all.tsv', data.table = F)
regions <- fread('/project/aid_sharing/AID_sharing/outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/06_genomic_regions_all_traits_munged/loci_all_traits.txt')

posterior$genomic_region <- rep('-', nrow(posterior))

for(i in 1:nrow(posterior)){
  
  posterior[i,]$genomic_region <- unique(regions[regions$pan_locus== posterior[i,]$locus, ]$pan_locus_name)

  }


final_post <- posterior %>% select(-c('locus')) %>% rename('genomic_region(CHR_START_END)'=genomic_region, 'trait_1'=t1, 'trait_2'=t2, 'leadSNP_locus_trait_1'=hit1, 'leadSNP_locus_trait_2'=hit2)

final_post$trait_1[final_post$trait_1=='f1'] <- 'Fgut'
final_post$trait_2[final_post$trait_2=='f1'] <- 'Fgut'

final_post$trait_1[final_post$trait_1=='f2'] <- 'Faid'
final_post$trait_2[final_post$trait_2=='f2'] <- 'Faid'

final_post$trait_1[final_post$trait_1=='f3'] <- 'Falrg'
final_post$trait_2[final_post$trait_2=='f3'] <- 'Falrg'

final_post$trait_1[final_post$trait_1=='derma'] <- 'eczema'
final_post$trait_2[final_post$trait_2=='derma'] <- 'eczema'


fwrite(final_post,'/project/aid_sharing/AID_sharing/outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/tables/Supplementary_table_coloc_all_loci_results.csv', sep = ',', col.names = T, quote = F)








