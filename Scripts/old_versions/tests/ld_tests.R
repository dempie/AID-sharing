library(LDlinkR)
library(data.table)

factor_loci <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/factor_loci_moloc_closest_gene_info.txt', data.table = F)


ldmatrix <- list()
for(i in 1:3){
  tt <- c('f1', 'f2','f3')[i]
  
     #cycle through the chromosomes
      for(k in 1:length(unique(factor_loci[factor_loci$trait==tt,]$chr))) {
            
            ch <- unique(factor_loci[factor_loci$trait==tt,]$chr)[k] 
            snps <- factor_loci[factor_loci$trait==tt & factor_loci$chr==ch, ]$SNP 
            
            #run only if there are at least 2 SNPS to compare
            if(length(snps)>=2){
                  ldmatrix[[tt]][[ch]]<- LDmatrix(snps = snps , token='d92873e28b8f')
            } else {
                  ldmatrix[[tt]][[ch]] <- paste0(snps, 'is the only one in this chromosome', collapse = '_')
            }
            
            
            }
  
}


saveRDS( ldmatrix ,'outputs/tests/ld_matrix_factor_loci.RDS')



#---- use LD to search for Heterogeniety loci ---------------------------------


#filter for only snps that are significant for P or Q


head(f1)
q_sign <- f1[f1$Q_CHISQ_PVAL< 5e-8, c('SNP', 'CHR', 'BP', 'Q_CHISQ_PVAL')]

dim(q_sign)

tr_sign <- list()
for(i in 1:3){
  tt <- c('f1', 'f2', 'f3')[i]
  tr_sign[[tt]] <- gw_q[[tt]][gw_q[[tt]]$P< 5e-8,]
}


matr <- list()
ld_matr <- list()
for(i in 1:3){
  tt <- c('f1', 'f2','f3')[i]
  
  #for each chromosome
  for(k in 1:length(unique(factor_loci[factor_loci$trait==tt,]$chr))) {
    ch <- unique(factor_loci[factor_loci$trait==tt,]$chr)[k] 
    
    #for each locus in the chromosome for that trait
    for(j in 1:nrow(factor_loci[factor_loci$chr==ch & factor_loci$trait==tt,])){
      
      lead <- factor_loci[factor_loci$trait==tt & factor_loci$chr==ch,]$SNP[j]
      lead_bp <- factor_loci[factor_loci$trait==tt & factor_loci$SNP==lead,]$BP
      
      SNPS <- q_sign[q_sign$CHR==ch & data.table::between(q_sign$BP,ifelse(lead_bp-1e6 < 0, 0,lead_bp-1e6 ) , lead_bp +1e6),]$SNP
      
      if(data.table::between(length(c(SNPS, lead)), 2, 10000)){
        #ld_matr[[tt]][[lead]]<- LDmatrix(c(SNPS, lead), token='d92873e28b8f')
        matr[[tt]][[lead]] <- c(c(SNPS, lead))
      } else {
        matr[[tt]][[lead]] <- ifelse(length(SNPS)<1, paste0('less than 1 SNP'), paste0('more than 1000') )
      }
      
    }
  }
}

saveRDS(ld_matr,'outputs/tests/ld_matrix_heterogeneity_factor.RDS')
saveRDS(matr,'outputs/tests/ld_matrix_heterogeneity_factor_snps.RDS')


































