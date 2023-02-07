library(data.table)
library(dplyr)

all_loci <-  fread('outputs/rev_1/coloc_eqtl_mr/loci_definitions/final_locus_table.tsv', data.table = F) 

all_loci$final.locus=paste0(all_loci$Chr,":",all_loci$start,"-",all_loci$end,"_",all_loci$sub_locus)

all_loci <- select(all_loci, -c('pan.locus', 'sub_locus')) %>% rename('LOCUS_NAME(CHR_START_END)'= final.locus)

head(all_loci)
 temp <- all_loci
temp[is.na(temp$othA),]


colnames(all_loci) <- toupper(colnames(all_loci))


all_loci$TRAIT[which(all_loci$TRAIT=='f1')] <- 'Fgut'
all_loci$TRAIT[which(all_loci$TRAIT=='f2')] <- 'Faid'
all_loci$TRAIT[which(all_loci$TRAIT=='f3')] <- 'Falrg'


head(all_loci)
table(all_loci$TRAIT)
table(temp[is.na(temp$othA),]$trait)




#use all_genomic regions table to add the alternate allele


asthma <- fread('outputs/rev_1/05_munge/munged/asthma_munged_build37.txt')
sle <- fread('outputs/rev_1/05_munge/munged/sle_munged_build37.txt')
cd <- fread('outputs/rev_1/05_munge/munged/cd_munged_build37.txt')
t1d <- fread('outputs/rev_1/05_munge/munged/t1d_munged_build37.txt')





nam <- c('asthma', 'sle',  'cd', 't1d')
files <- list(asthma, sle, cd, t1d)
final <- list()

#cycle through the gwas 
for(i in 1:4){
 #select the active gwas select the nas
       
       act <- files[[i]] 
       qq <- act[which(act$SNP %in% temp[which(is.na(temp$othA) & temp$trait==nam[i]),]$SNP),]
       ttmpo <- temp[which(is.na(temp$othA) & temp$trait==nam[i]),]
      
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
colnames(ref) <- toupper(colnames(ref))

qq <- inner_join(all_loci,ref, c('SNP', 'TRAIT'), all.x=T)

qq <- select(qq, -c('OTHA.x', 'CHR.y', 'REFA.y'))

qq <- rename(qq, 'REFA'=REFA.x, 'OTHA'=OTHA.y, 'CHR'=CHR.x)
head(qq)

all_loci_noNA <- all_loci[which(!is.na(all_loci$OTHA)),]

dim(all_loci_noNA)


all_loci_miss <- rbind(all_loci_noNA, qq)
dim(all_loci_miss)

fwrite(all_loci_miss , 'outputs/rev_1/Supplementary_tables/8_Supplementary_all_conditional_loci.csv', sep = ',', col.names = T, quote = F)



# head(all_loci_ok)
# 
# 
# table(all_loci$othA == all_loci_ok$OTHA)
# 
# 
# 
# all_loci
# table(is.na(all_loci$othA == all_loci_ok$OTHA))








