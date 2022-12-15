#blagoje grant
library(data.table)
coclo <- fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/tables/Supplementary_table_eQTL_results.csv')


bcells <- unique(coclo[coclo$PP.H4.ABF>0.90 & coclo$CELL_TYPE %in% c("BIN", "BMem", 'Plasma'), c('TRAIT_1', 'TRAIT_2' ,'CELL_TYPE')])
tcells <- unique(coclo[coclo$PP.H4.ABF>0.90 & coclo$CELL_TYPE %in% c('CD4NC', "CD8NC", "CD8S100B", 'CD4ET', 'CD8ET', 'CD4SOX4'), c('TRAIT_1', 'TRAIT_2' ,'CELL_TYPE')])

fwrite(bcells, 'outputs/tests/bcells_0.90_eqtls.csv', sep=',', col.names = T, row.names = F)
fwrite(tcells, 'outputs/tests/tcells_0.90_eqtls.csv', sep=',', col.names = T, row.names = F)

#-------------
library(stringr)
library(data.table)
library(dplyr)


setwd('/project/aid_sharing/sc_eqtl/')

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


dim(qtl_res)


traits <- unique(qtl_res[qtl_res$PP.H4.abf>0.90, ]$t1)
coclo <- qtl_res
bcells <- list()
tcells <- list()
for(tt in traits){
bcells[[tt]] <- length(unique(coclo[coclo$PP.H4.abf>0.90 & coclo$cell_type %in% c("BIN", "BMem", 'Plasma') & coclo$t1==tt, ]$t2))
tcells[[tt]] <- length(unique(coclo[coclo$PP.H4.abf>0.90 & coclo$cell_type %in% c('CD4NC', "CD8NC", "CD8S100B", 'CD4ET', 'CD8ET', 'CD4SOX4') & coclo$t1==tt ,]$t2))
}

bcells
tcells


bcells_t <- list()
tcells_t <- list()
all_traits <- 
factors <- c('f1', 'f2', 'f3')
c(all_traits, factors)


        bcells_t <- length(unique(coclo[coclo$PP.H4.abf>0.90 & coclo$cell_type %in% c("BIN", "BMem", 'Plasma') & coclo$t1 %in%traits[!traits %in% c('f1', 'f2', 'f3')], ]$t2))
        tcells_t <- length(unique(coclo[coclo$PP.H4.abf>0.90 & coclo$cell_type %in% c('CD4NC', "CD8NC", "CD8S100B", 'CD4ET', 'CD8ET', 'CD4SOX4') & coclo$t1 %in% traits[!traits %in% c('f1', 'f2', 'f3')] ,]$t2))

        
        bcells_f <- length(unique(coclo[coclo$PP.H4.abf>0.90 & coclo$cell_type %in% c("BIN", "BMem", 'Plasma') & coclo$t1 %in% c('f1', 'f2', 'f3'), ]$t2))
        tcells_f <- length(unique(coclo[coclo$PP.H4.abf>0.90 & coclo$cell_type %in% c('CD4NC', "CD8NC", "CD8S100B", 'CD4ET', 'CD8ET', 'CD4SOX4') & coclo$t1  %in% c('f1', 'f2', 'f3') ,]$t2))
        
        
        
        bcells_f_not_t <- length(setdiff(unique(coclo[coclo$PP.H4.abf>0.90 & coclo$cell_type %in% c("BIN", "BMem", 'Plasma') & coclo$t1 %in% c('f1', 'f2', 'f3'), ]$t2),
                                  unique(coclo[coclo$PP.H4.abf>0.90 & coclo$cell_type %in% c("BIN", "BMem", 'Plasma') & coclo$t1 %in%traits[!traits %in% c('f1', 'f2', 'f3')], ]$t2) ))
        
        tcells_f_not_t <- length(setdiff(unique(coclo[coclo$PP.H4.abf>0.90 & coclo$cell_type %in% c('CD4NC', "CD8NC", "CD8S100B", 'CD4ET', 'CD8ET', 'CD4SOX4') & coclo$t1  %in% c('f1', 'f2', 'f3') ,]$t2),
                                  unique(coclo[coclo$PP.H4.abf>0.90 & coclo$cell_type %in% c('CD4NC', "CD8NC", "CD8S100B", 'CD4ET', 'CD8ET', 'CD4SOX4') & coclo$t1 %in% traits[!traits %in% c('f1', 'f2', 'f3')] ,]$t2)
                                  
                        ))
        
print(paste0('eGENEs in B cells in all traits:  ' , bcells_t))
print(paste0('eGENEs in B cells in factors  ' , bcells_f))
print(paste0('eGENEs in T cells in all traits:  ' , tcells_t))
print(paste0('eGENEs in T cells in factors  ' , tcells_f))
print(paste0('eGENEs in B cells in factors but not in traits:  ' , bcells_f_not_t))
print(paste0('eGENEs in T cells in factors but not in traits:  ' , tcells_f_not_t))




unique(coclo[coclo$PP.H4.abf>0.90 & coclo$cell_type %in% c('CD4NC', "CD8NC", "CD8S100B", 'CD4ET', 'CD8ET', 'CD4SOX4') & coclo$t1  %in% c('f1', 'f2', 'f3') ,]$t2)








