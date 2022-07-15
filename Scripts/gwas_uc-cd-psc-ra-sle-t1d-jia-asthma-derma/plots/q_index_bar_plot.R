library(data.table)
library(RColorBrewer)
library(ggplot2)
library(GenomicRanges)
library(LDlinkR)


factor_loci <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/factor_loci_moloc_closest_gene_info.txt', data.table = F)
#------------------------------Q index locus breaker ---------------------------

gw_q <- list()
for(i in c(1:3)){
  tt <- c('f1', 'f2', 'f3')[i]
  gw_q[[tt]] <- fread(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/q_index_munged/',tt,'_munged_q_index_build37.txt' ),
                      data.table = F)
}


f1 <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/q_index_munged/f1_munged_q_index_build37.txt', data.table = F)
a <- locus.breaker(f1, p.label = 'Q_CHISQ_PVAL')
head(f1)
a$check <- rep(NA, nrow(a))
a <- locus.breaker(f1, p.label = 'Q_CHISQ_PVAL')
for(i in 1:nrow(a)){
      st <- a$start[i]
      en <- a$end[i]
      chrom <- a$chr[i]
      a$check[i] <- any(f1[f1$CHR==chrom & between(f1$BP, st, en),]$Q_CHISQ_PVAL< 5e-8)
}




a_range <- GRanges(seqnames = a$chr, IRanges(start = as.numeric(a$start), end = as.numeric(a$end)))
f_range <- GRanges(seqnames = factor_loci$chr, IRanges(start = as.numeric(factor_loci$start), end=as.numeric(factor_loci$end)))

ovl <- findOverlaps(f_range, a_range)


factor_loci$q_signif_loci_ovl <- rep(FALSE, nrow(factor_loci)) 
factor_loci$q_signif_loci_ovl[ovl@from] <- TRUE






factor_loci$trait <- toupper(factor_loci$trait) 

ggplot(data=factor_loci, aes(x=trait, fill=q_signif_loci_ovl) ) + 
      geom_bar(stat='count',position = position_stack(reverse=T),  color='black')+ 
      scale_fill_manual(values = c("grey80", "white")) + 
      geom_text(aes(label = paste0("n=", ..count..)),position= position_stack(vjust = 0.5, reverse=T),stat='count')+
      labs(y = 'Number of loci', x = '')+
      theme_classic() +
      theme(legend.position="bottom")


#instead of looking at the overlapp, look if there are SNPs inisde the locus are have significant Q snp
# factor_loci$het_method_2 <- rep(NA, nrow(factor_loci))
# for(i in 1:nrow(factor_loci)){
#       s <- factor_loci$start[i]
#       e <- factor_loci$end[i]
#       chr <- factor_loci$chr[i]
#       factor_loci$het_method_2[i] <- any(gw_q$f1[gw_q$f1$CHR==chr & between(gw_q$f1$BP, s,e ),]$Q_CHISQ_PVAL <= 5e-8)
# }
# 











































