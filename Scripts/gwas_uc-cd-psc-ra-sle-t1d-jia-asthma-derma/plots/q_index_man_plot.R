library(data.table)
library(qqman)
library(dplyr)
library(RColorBrewer)
library(miamiplot)


#manathan plots
gw <- list()
g <- list()
q <- list()
for(i in c(1:3)){
  tt <- c('f1', 'f2', 'f3')[i]
  gw[[tt]] <- fread(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/',tt,'_munged_build37.txt' ),
                    data.table = F)
  # g[[tt]] <-  gw[[tt]][gw[[tt]]$Pval_Estimate < 0.5 & (!is.na(gw[[tt]]$Pval_Estimate)),   ]
  # g[[tt]] <-  gw[[tt]][gw[[tt]]$Pval_Estimate < 0.05 & (!is.na(gw[[tt]]$Pval_Estimate)),   ]
}


#add the Q index 
q_index <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/factor1_gwas_final_withQindex.txt', data.table = F) 
q_index <- q_index[q_index$SNP %in% gw$f1$SNP,] #keep only the snps that passed the munging 
for(i in 1:3){
  tt <- c('f1', 'f2', 'f3')[i]
  gw[[tt]]$Q_chisq_pval <- rep(NA, nrow(gw[[tt]]))
  gw[[tt]][match(q_index$SNP ,gw[[tt]]$SNP),]$Q_chisq_pval <- q_index$Q_chisq_pval
}


#prepare the data for plotting, merge 
q_miami <- data.frame('rsid'= q_index$SNP,'chr'=q_index$CHR, 'pos'= q_index$BP, 'study'= rep('q', nrow(q_index)),'p'=q_index$Q_chisq_pval )
plot_miami <- list()
for(i in 1:3){
    tt <- c('f1', 'f2', 'f3')[i]
    plot_miami[[tt]] <- data.frame('rsid'= gw[[tt]]$SNP,'chr'=gw[[tt]]$CHR, 'pos'= gw[[tt]]$BP, 'study'=gw[[tt]]$LHS, 'p'= gw[[tt]]$P ) 
    plot_miami[[tt]] <- rbind(plot_miami[[tt]], q_miami)
}


ggmiami(data = gwas_results, 
        split_by = "study", split_at = 'A', p = "pval",  upper_ylab = "Study A", lower_ylab = "Study B", suggestive_line = F, chr_colors = NULL,
        upper_chr_colors = brewer.pal(12, 'Paired')[c(1,2)],lower_chr_colors = c('black', 'grey'),  genome_line_color = "black" )


ggmiami(data = plot_miami$f1, 
        split_by = "study", split_at = 'study', p = "p",  upper_ylab = "Study A", lower_ylab = "Study B", suggestive_line = F, chr_colors = NULL,
        upper_chr_colors = brewer.pal(12, 'Paired')[c(1,2)],lower_chr_colors = c('black', 'grey'),  genome_line_color = "black" )


