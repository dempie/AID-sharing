library(GenomicSEM)
library(data.table)

aid_sumstats <- readRDS('/project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/04_sumstat_outputs/sumstats_output/sumstats_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS')
ldsc_model <- readRDS('/project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/04_sumstat_outputs/ldsc_output/ldsc_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS')
loci.table <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/loci_table_names_nicola_genes.csv', data.table = F) 


loci_snps  <- aid_sumstats[aid_sumstats$SNP %in% loci.table$SNP,]

aid_model_f1 <- 'F1 =~ NA*crohn + uc  + psc  
F2 =~ NA*jia + sle + ra+ t1d 
F3 =~ NA*asthma+ derma 

F1~~F2
F1~~F3
F2~~F3
F1 ~~ 1*F1
F2 ~~ 1*F2
F3 ~~ 1*F3



F2 ~ SNP
F3 ~ SNP

crohn + uc  + psc ~ SNP

derma~~a*derma
a>0.001'


aid_model_f2 <- 'F1 =~ NA*crohn + uc  + psc  
F2 =~ NA*jia + sle + ra+ t1d 
F3 =~ NA*asthma+ derma 

F1~~F2
F1~~F3
F2~~F3
F1 ~~ 1*F1
F2 ~~ 1*F2
F3 ~~ 1*F3



F1 ~ SNP
F3 ~ SNP

jia + sle + ra+ t1d ~ SNP

derma~~a*derma
a>0.001'


aid_model_f3 <- 'F1 =~ NA*crohn + uc  + psc  
F2 =~ NA*jia + sle + ra+ t1d 
F3 =~ NA*asthma+ derma 

F1~~F2
F1~~F3
F2~~F3
F1 ~~ 1*F1
F2 ~~ 1*F2
F3 ~~ 1*F3



F1 ~ SNP
F2 ~ SNP

asthma + derma~ SNP

derma~~a*derma
a>0.001'


output <- list()

for(i in 1:3){
  tt <-  c('f1', 'f2', 'f3')[i]
  output[[tt]] <- userGWAS(covstruc = ldsc_model, 
                           SNPs = loci_snps, 
                           model = c(aid_model_f1, aid_model_f2, aid_model_f3)[i], 
                           sub=c('F1~~F2'),
                           parallel = T,cores = 2 )
}

saveRDS(output, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/10_Q_index_estimation/output_model_estimation.RDS')






