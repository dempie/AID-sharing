# in this script all the necessary steps for the colocalizzation will be performed 

library(GenomicSEM)
library(data.table)
library(dplyr)
library(MungeSumstats)

#the idea is to use the sumstat function to flip the alleles and returned a partially munged sumstat list
#however munge will do something to the standard errors, so the standard errors will be added by me after munge has done

#all the sumstats in the gwas have Beta, except for psc and asthma that has OR, just had the beta column

#psc
psc <-  fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/psc_ji-2016.txt', data.table = F)
head(psc)
psc$beta <- log(psc$OR)
psc <- rename(psc, the_effect_in_OR=OR )
fwrite(psc,'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/summary_stats_with_beta/psc_ji-2016_withBeta.txt',
       sep = '\t', col.names = T, row.names = F, quote = F)

#asthma
asthma <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/asthma_ban-2020.txt', data.table = F)
head(asthma)
asthma$beta <- log(asthma$OR)
asthma <- rename(asthma, the_effect_in_OR=OR)
fwrite(psc,'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/summary_stats_with_beta/asthma_han-2020_withBeta.txt',
       sep = '\t', col.names = T, row.names = F, quote = F)


#prepare the  factor GWAS
f1 <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/factor1_gwas_final_withQindex.txt', data.table = F)
f1 <- rename(f1, c(beta='est', p='Pval_Estimate')) %>% select(c('beta', 'p', 'SE', 'SNP', 'CHR', 'BP'))
f2 <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/factor2_gwas_final_withQindex.txt', data.table = F)
f2 <- rename(f2, c(beta='est', p='Pval_Estimate')) %>% select(c('beta', 'p', 'SE', 'SNP', 'CHR', 'BP'))
f3 <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/factor3_gwas_final_withQindex.txt', data.table = F)
f3 <- rename(f3, c(beta='est', p='Pval_Estimate')) %>% select(c('beta', 'p', 'SE', 'SNP', 'CHR', 'BP'))

f<- list(f1,f2,f3) 
for(i in c(1:3)){
  fwrite(f[i], paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/summary_stats_with_beta/f',i,'_summary_stats_for_munge.txt'), 
         sep = '\t', col.names = T, row.names = F, quote = F)
}

#-------------------------------------------------------------------------------


file_names_ok <- c('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/summary_stats_with_beta/f1_summary_stats_for_munge.txt', 
                   'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/summary_stats_with_beta/f2_summary_stats_for_munge.txt', 
                   'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/summary_stats_with_beta/f3_summary_stats_for_munge.txt', 
                    'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/t1d_chiou-2021.txt', 
                   'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/crohn_delange-2017.txt', 
                   'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/uc_delange-2017.txt',
                   'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/summary_stats_with_beta/psc_ji-2016_withBeta.txt',
                   'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/jia_beta_lopezisac-2020.txt',
                   'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/sle_beta_bentham-2015.txt',
                   'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/ra_eu_okada-2014.txt',
                   'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/summary_stats_with_beta/asthma_han-2020_withBeta.txt',
                   'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/derma_sliz-2021.txt'
)


res <- sumstats(files = file_names_ok , ref = 'SNP/reference.1000G.maf.0.005.txt.gz',
         trait.names = c('f1','f2','f3','t1d','crohn', 'uc', 'psc', 'jia', 'sle', 'ra',  'asthma', 'derma') ,
         se.logit =  c(F,F,F,F,F,F,F,F,F),  
         OLS = c(T,T,T,T,T,T,T,T,T), 
         linprob = c(F,F,F,F,F,F,F,F,F), 
         N = c(NA, NA, NA, NA, NA, NA, NA, NA, NA), 
         parallel = T,
         keep.indel= F 
)

saveRDS(res, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/sumstats_output/intermediate_sumstats.RDS')
res <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/sumstats_output/intermediate_sumstats.RDS')
system('mv *log /project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/sumstats_output')





#----------------check if the betas are identical and add the standard errors from the original GWAS

f1 <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/summary_stats_with_beta/f1_summary_stats_for_munge.txt', data.table = F)
f2 <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/summary_stats_with_beta/f2_summary_stats_for_munge.txt', data.table = F)
f3 <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/summary_stats_with_beta/f3_summary_stats_for_munge.txt', data.table = F)
uc <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/uc_delange-2017.txt', data.table = F)
cd <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/crohn_delange-2017.txt', data.table = F)
psc <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/psc_ji-2016.txt', data.table = F)
jia <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/jia_beta_lopezisac-2020.txt', data.table = F)
sle <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/sle_beta_bentham-2015.txt', data.table = F)
t1d <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/t1d_chiou-2021_build37.txt', data.table=F) 
asthma <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/asthma_ban-2020.txt', data.table = F)
ra <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/ra_eu_okada-2014_chr_bp.txt', data.table = F)
derma <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/derma_sliz-2021_build37.txt', data.table = F)

#check if the betas in res are equal to the original gwas is abs() values as they should have been flipped and munged

dim(res)#3344679      24

head(res)



psc <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/summary_stats_with_beta/psc_ji-2016_withBeta.txt', data.table = F)





for(i in c(1:7)){
  
trait <- gwas[i]
colnames(trait) <- toupper(colnames(trait))

trait <- (trait[(trait$SNP %in% res$SNP),  ])

res <- res[order(as.numeric(res$CHR), as.numeric(res$BP)),]
trait <- trait[order(as.numeric(unlist(select(trait, 'CHR')), as.numeric(unlist(select(trait,'BP'))))) ,]
num_row[i] <- nrow(trait)==nrow(res)
 
 
beta_identical[i] <- mean(abs(res$beta.psc)==abs(select(trait, 'BETA')))
summary[[i]][[1]] <- summary(abs(res$beta.psc))
summary[[i]][[2]] <- summary(abs(select(trait, 'BETA')))
res$se.psc <-select(trait, 'SE')
 
se_identical[[i]] <- mean(res$se.psc == select(trait, 'SE') )

}

 head(res)
 head(trait)


 
 
 
 
 
 
 
 
 










