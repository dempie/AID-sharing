library(GenomicSEM)
library(data.table)
library(ggplot2)

aid_sumstats <- readRDS('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/02_sumstats_function/sumstats_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS')
ldsc_model <- readRDS('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/02_sumstats_function/ldsc_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS')
factor_loci <-fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/06_genomic_regions_all_traits_munged/factor.regions.txt', data.table = F)


loci_snps  <- aid_sumstats[aid_sumstats$SNP %in% factor_loci$SNP,]

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


crohn~~a*crohn
uc~~b*uc
psc~~c*psc
jia~~d*jia
sle~~e*sle
ra~~f*ra
t1d~~g*t1d
asthma~~h*asthma
derma~~i*derma

a>0.001
b>0.001
c>0.001
d>0.001
e>0.001
f>0.001
g>0.001
h>0.001
i>0.001
'


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

crohn~~a*crohn
uc~~b*uc
psc~~c*psc
jia~~d*jia
sle~~e*sle
ra~~f*ra
t1d~~g*t1d
asthma~~h*asthma
derma~~i*derma

a>0.001
b>0.001
c>0.001
d>0.001
e>0.001
f>0.001
g>0.001
h>0.001
i>0.001'


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

crohn~~a*crohn
uc~~b*uc
psc~~c*psc
jia~~d*jia
sle~~e*sle
ra~~f*ra
t1d~~g*t1d
asthma~~h*asthma
derma~~i*derma

a>0.001
b>0.001
c>0.001
d>0.001
e>0.001
f>0.001
g>0.001
h>0.001
i>0.001'


output <- list()

for(i in 1:3){
  tt <-  c('f1', 'f2', 'f3')[i]
  output[[tt]] <- userGWAS(covstruc = ldsc_model, 
                           SNPs = loci_snps, 
                           model = c(aid_model_f1, aid_model_f2, aid_model_f3)[i], 
                           sub=c('F1~~F2'),
                           parallel = T,cores = 2 )
}

saveRDS(output, 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_qsnp/output_model_estimation.RDS')


# #-------- calculation of Q index for the lead SNP ------------------------------
q_ind <- readRDS('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_qsnp/output_model_estimation.RDS')

fac <- readRDS('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/03_gwas_estimation/factors_summary_stats.RDS')

head(fac)
factor_loci$Q_chisq_pval <- rep(NA, nrow(factor_loci))



for(i in 1:3){
  tt <- c('f1', 'f2', 'f3')[i]
  q_act <- q_ind[[tt]][[1]]
  q_act_2 <- q_act[q_act$SNP %in% factor_loci[factor_loci$trait==tt,]$SNP,]
  gwas <- fac[[1]][fac[[1]]$SNP %in% q_act_2$SNP, ]
  
  #put all of them in the same order
  gwas_2 <- gwas[match(q_act_2$SNP, gwas$SNP),]
  
  #substract the Chisquare
  Q <- gwas_2$chisq - q_act_2$chisq
  df <- gwas_2$chisq_df - q_act_2$chisq_df
  #put all of them in the same order
  factor_loci[factor_loci$trait==tt,][match(q_act_2$SNP, factor_loci[factor_loci$trait==tt,]$SNP, nomatch = NA),]$Q_chisq_pval<- pchisq(Q,df,lower.tail=FALSE)
  print(table(pchisq(Q,df,lower.tail=FALSE)<5e-8))
  
}


factor_loci$het<- ifelse(factor_loci$Q_chisq_pval<5e-8, TRUE, FALSE)

fwrite(factor_loci, 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_qsnp/genomic_regions_table_heterogeneity_column.txt', sep = '\t', col.names = T)


pdf('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_qsnp/heterogeneity_index_genomic_regions.pdf', width = 10, height = 10)
ggplot(data=factor_loci, aes(x=trait, fill=het) ) +
  geom_bar(stat='count',position = position_stack(reverse=T),  color='black')+
  scale_fill_manual(values = c("grey80", "white")) +
  geom_text(aes(label = paste0("n=", ..count..)),position= position_stack(vjust = 0.5, reverse=T),stat='count')+
  labs(y = 'Number of genomic regions', x = '')+
  theme_classic() +
  theme(legend.position="bottom") + ggtitle('Q heterogeneoity index lead SNPs genomic regions')
dev.off()












