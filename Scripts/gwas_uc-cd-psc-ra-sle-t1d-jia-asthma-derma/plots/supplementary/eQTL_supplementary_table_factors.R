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


mr_res<- Reduce(rbind, mr_ok)
mr_factors <- mr_res[mr_res$trait%in% c('f1', 'f2', 'f3')]
#add the locus name to the table
loci <- fread('/project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/loci_table_names_nicola_genes.csv', data.table = F)
loci$tmp <- paste0(loci$pan.locus, '_',loci$sub_locus)
head(loci)

mr_factors$locus_name <- rep('', nrow(mr_factors))
for(i in 1:nrow(mr_factors)){
  mr_factors[i, ]$locus_name <- loci[loci$tmp==mr_factors[i, ]$locus,]$final.locus[1]
  
}

mr_final <- mr_factors %>% rename( 'locus_name_(chr:start_end)'=locus_name, 'cell_type'=cell.type) %>% select(-c('locus'))

mr_final$trait[mr_final$trait=='f1'] <- 'Fgut'
mr_final$trait[mr_final$trait=='f2'] <- 'Faid'
mr_final$trait[mr_final$trait=='f3'] <- 'Falrg'

names(mr_final) <- toupper(names(mr_final))

#save
fwrite(mr_final,'/project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/supplementary/supplementary_tables/Supplementary_table_3_MR_results.csv', sep = ',', col.names = T, quote = F)



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

qtl_final$TRAIT_1[qtl_final$TRAIT_1=='f1'] <- 'Fgut'
qtl_final$TRAIT_1[qtl_final$TRAIT_1=='f2'] <- 'Faid'
qtl_final$TRAIT_1[qtl_final$TRAIT_1=='f3'] <- 'Falrg'

colnames(qtl_final) <- toupper(colnames(qtl_final))

head(qtl_final)

#save
fwrite(qtl_final,'/project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/supplementary/supplementary_tables/Supplementary_table_4_eQTL_results.csv', sep = ',', col.names = T, quote = F)
# 
# final_table<- do.call(rbind, chunk_ok)
# only_f <- final_table[final_table$V9 %in% c('f1', 'f2', 'f3'),]
# 
# only_f_H4 <- only_f[as.numeric(only_f$V6)>0.90, ]
# 
# unique(only_f_H4[only_f_H4$V9=='f1',]$V10)
# 
# f1_h4<- only_f_H4[only_f_H4$V9=='f1',]
# 
# f1_h4[duplicated(f1_h4$V10), ]
# 
# length(unique(f1_h4$V10))
# 
# 
# rownames(f1_h4) <- NULL
# dim(f1_h4[, c(9,10,11)][!duplicated(f1_h4[, c(9,10,11)]),])
# 
# f1_h4
# 
# 
# f2_h4<- only_f_H4[only_f_H4$V9=='f2',]
# 
# dim(f2_h4)
# 
# dim(f1_h4[!duplicated(f1_h4[, c(9,10,11)]),])
# 
# 
# f3_h4<- only_f_H4[only_f_H4$V9=='f3',]
# 
# 
# 
# dim(f2_h4[!duplicated(f2_h4[, c(9,10,11)]),])
# dim(f3_h4[!duplicated(f3_h4[, c(9,10,11)]),])
# 
# 
# dim(final_table[round(as.numeric(final_table$V6) ,2)>0.9 & final_table$V9%in%c("f1","f2","f3"),])
# 
# 
# factors <- final_table[round(as.numeric(final_table$V6) ,2)>0.9 & final_table$V9%in%c("f1","f2","f3"),]
# 
# dim(factors[factors$V9=='f1',][!duplicated(factors[factors$V9=='f1', c(9,10,11)]),])
# dim(factors[factors$V9=='f2',][!duplicated(factors[factors$V9=='f2', c(9,10,11)]),])
# dim(factors[factors$V9=='f3',][!duplicated(factors[factors$V9=='f3', c(9,10,11)]),])
# 
# 



