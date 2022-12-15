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

#------- all conditional loci --------------------------------------------------

all_loci <-  fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/09_loci_definitions_Nicola/loci_definitions/final_locus_table.tsv', data.table = F) 

all_loci$final.locus=paste0(all_loci$Chr,":",all_loci$start,"-",all_loci$end,"_",all_loci$sub_locus)

all_loci <- select(all_loci, -c('pan.locus', 'sub_locus')) %>% rename('LOCUS_NAME(CHR_START_END)'= final.locus)

# head(all_loci)
# temp <- all_loci[, c('Chr',  'start', 'end', 'SNP', 'bp' ,'refA', 'othA', 'trait', 'p', 'pJ', 'b', 'LD_r', 'final.locus', 'pan.locus')]
# temp[is.na(temp$othA),]
# all_loci[all_loci$pan.locus==238,]

colnames(all_loci) <- toupper(colnames(all_loci))


all_loci$TRAIT[all_loci$TRAIT=='f1'] <- 'Fgut'
all_loci$TRAIT[all_loci$TRAIT=='f2'] <- 'Faid'
all_loci$TRAIT[all_loci$TRAIT=='f3'] <- 'Falrg'


head(all_loci)
table(all_loci$TRAIT)





#use all_genomic regions table to add the alternate allele


asthma <- fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/munged/asthma_munged_build37.txt')
sle <- fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/munged/sle_munged_build37.txt')
uc <- fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/munged/uc_munged_build37.txt')
cd <- fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/munged/cd_munged_build37.txt')
t1d <- fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/munged/t1d_munged_build37.txt')





nam <- c('asthma', 'sle', 'uc', 'cd', 't1d')
files <- list(asthma, sle, uc, cd, t1d)
final <- list()

#cycle through the gwas 
for(i in 1:5){
 #select the active gwas select the nas
       
       act <- files[[i]] 
       qq <- act[act$SNP %in% temp[is.na(temp$othA) & temp$trait==nam[i],]$SNP,]
       ttmpo <- temp[is.na(temp$othA) & temp$trait==nam[i],]
      
              jj<- inner_join(qq, ttmpo, 'SNP')
      
              for( k in 1:nrow(jj)){
                
                if(jj$refA[k]==jj$A1[k]){ 
                  jj$othA[k] = jj$A2[k] 
                    
                    } else if( jj$refA[k]==jj$A2[k] ) {
                      jj$othA[k] = jj$A1[k] 
                    }
                
              }
              
        final[[nam[i]]] <- jj[, c('SNP', 'CHR', 'trait', 'refA', 'othA')]     
      

}


ref<-do.call(rbind, final)


for(k in 1:nrow(all_loci[is.na(all_loci$OTHA), ])){
  trait <-  all_loci[is.na(all_loci$OTHA), ][i,]$TRAIT
  snp <- all_loci[is.na(all_loci$OTHA), ][i,]$SNP
   all_loci[is.na(all_loci$OTHA), ][i,]$OTHA <- ref[ref$trait==trait & ref$SNP==snp, ]$othA
   
}



all_loci_ok <- all_loci
fwrite(all_loci_ok , 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/tables/Supplementary_all_conditional_loci.csv', sep = ',', col.names = T, quote = F)



# head(all_loci_ok)
# 
# 
# table(all_loci$othA == all_loci_ok$OTHA)
# 
# 
# 
# all_loci
# table(is.na(all_loci$othA == all_loci_ok$OTHA))








