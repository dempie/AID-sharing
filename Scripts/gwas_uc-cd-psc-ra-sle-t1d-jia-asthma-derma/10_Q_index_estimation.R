library(GenomicSEM)
library(data.table)
library(ggplot2)

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


#-------- calculation of Q index for the lead SNP ------------------------------
q_ind <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/10_Q_index_estimation/output_model_estimation.RDS') 

fac <- list()
for(i in 1:3){
  tt<- c('f1', 'f2', 'f3')[i]
  fac[[tt]]<- fread(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/factor',i,'_gwas_final_withQindex.txt' ))
}

loci.table$Q_chisq_pval <- rep(NA, nrow(loci.table))



for(i in 1:3){
  tt <- c('f1', 'f2', 'f3')[i]
  q_act <- q_ind[[tt]][[1]]
  q_act_2 <- q_act[q_act$SNP %in% loci.table[loci.table$trait==tt,]$SNP,]
  gwas <- fac[[tt]][fac[[tt]]$SNP %in% q_act_2$SNP, ]
  
  #put all of them in the same order
  gwas_2 <- gwas[match(q_act_2$SNP, gwas$SNP),]
  
  #substract the Chisquare
  Q <- gwas_2$chisq - q_act_2$chisq
  df <- gwas_2$chisq_df - q_act_2$chisq_df
  #put all of them in the same order
  loci.table[loci.table$trait==tt,][match(q_act_2$SNP, loci.table[loci.table$trait==tt,]$SNP, nomatch = NA),]$Q_chisq_pval<- pchisq(Q,df,lower.tail=FALSE)
  print(table(pchisq(Q,df,lower.tail=FALSE)<5e-8))
  
}


loci.table$het<- ifelse(loci.table$Q_chisq_pval<5e-8, TRUE, FALSE)

fwrite(loci.table, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/10_Q_index_estimation/loci_table_nicola_heterogeneity_column.txt', sep = '\t', col.names = T)

#plot the heterogeneity stakced bar plot
a<- ggplot(data=loci.table, aes(x=trait, fill=het) ) + 
  geom_bar(stat='count',position = position_stack(reverse=T),  color='black')+ 
  scale_fill_manual(values = c("grey80", "white")) + 
  geom_text(aes(label = paste0("n=", ..count..)),position= position_stack(vjust = 0.5, reverse=T),stat='count')+
  labs(y = 'Number of loci', x = '')+
  theme_classic() +
  theme(legend.position="bottom") + ggtitle('Conditional loci')



factor_loci <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/factor_loci_moloc_closest_gene_info.txt', data.table = F)


b<- ggplot(data=loci.table[loci.table$SNP %in% factor_loci$SNP,], aes(x=trait, fill=het) ) + 
  geom_bar(stat='count',position = position_stack(reverse=T),  color='black')+ 
  scale_fill_manual(values = c("grey80", "white")) + 
  geom_text(aes(label = paste0("n=", ..count..)),position= position_stack(vjust = 0.5, reverse=T),stat='count')+
  labs(y = 'Number of loci', x = '')+
  theme_classic() +
  theme(legend.position="bottom") + ggtitle('GWAS loci')






######à a function to test if eveyrhting is ok with the Q het 
test <- function(){
  x <- sample(nrow(loci.table),1)
  snp <- loci.table[x,]$SNP
  tt <- loci.table[x, ]$trait
  Q <- fac[[tt]][fac[[tt]]$SNP==snp,]$chisq  - q_ind[[tt]][[1]][q_ind[[tt]][[1]]$SNP==snp,]$chisq
  df <- fac[[tt]][fac[[tt]]$SNP==snp,]$chisq_df  - q_ind[[tt]][[1]][q_ind[[tt]][[1]]$SNP==snp,]$chisq_df
  print(pchisq(Q,df,lower.tail=FALSE) == loci.table[loci.table$SNP==snp& loci.table$trait==tt ,]$Q_chisq_pval)
  
  
}

for(i in 1:500){test()}  













