library(data.table)
library(GenomicSEM)
library(tidyr)
library(dplyr)
library(corrplot)



prepare_munge <- function(sum_stats, rsID, the_effect_allele, the_non_effect_allele, pvalue, the_OR=NA, the_Beta=NA, the_SE=NA, the_chr=NA, the_bp=NA, to_remove=NA, path = NA, the_INFO=NA, the_MAF=NA){
  #an error if arguments are not provided 
  if (missing(sum_stats) | missing(rsID) | missing(the_effect_allele) | missing(the_non_effect_allele) |missing(pvalue) ) {
    
    stop( 'At least one argument is missing')
    
  } else {
    
    require(dplyr)
    require(data.table)
    sum_stats <- sum_stats  %>% rename(c(SNP = all_of(rsID), A1 = all_of(the_effect_allele), A2 = all_of(the_non_effect_allele), p = all_of(pvalue) ))
    sum_stats$SNP <- tolower(sum_stats$SNP)
    sum_stats$p <- as.numeric(sum_stats$p)
    
    
    #conditional options 
    #remove columns
    if(!is.na(to_remove[1])){sum_stats <- select(sum_stats,-(all_of(to_remove)))} 
    
    #rename SE column
    if(!is.na(the_SE)){ 
      sum_stats <- rename(sum_stats, SE=all_of(the_SE))
      sum_stats$SE <- as.numeric(sum_stats$SE)
    }
    
    #rename the effect column
    if(!is.na(the_OR)){
      sum_stats <- rename(sum_stats, OR=all_of(the_OR))
      sum_stats$OR <- as.numeric(sum_stats$OR)
    }
    
    #rename the effect column
    if(!is.na(the_MAF)){
      sum_stats <- rename(sum_stats, MAF=all_of(the_MAF))
      sum_stats$MAF <- as.numeric(sum_stats$MAF)
    }
    
    #rename the effect column
    if(!is.na(the_INFO)){
      sum_stats <- rename(sum_stats, INFO=all_of(the_INFO))
      sum_stats$INFO <- as.numeric(sum_stats$INFO)
    }
    
    if (!is.na(the_Beta)){
      sum_stats <- rename(sum_stats, Beta=all_of(the_Beta))
      sum_stats$Beta <- as.numeric(sum_stats$Beta)
    } 
    
    if(is.na(the_OR) & is.na(the_Beta) ) {stop('Effect column not specified ')}
    
    
    #rename the CHR column
    if(!is.na(the_chr)){ sum_stats <-  rename(sum_stats, CHR=all_of(the_chr))}
    
    #rename the BP column
    if(!is.na(the_bp)){ sum_stats <- rename(sum_stats, BP=all_of(the_bp))}
    
    #save the file if a path is provided
    if(is.na(path)){
      invisible(sum_stats)
      
      
    } else {
      
      fwrite(sum_stats, path, sep = '\t', col.names = T, row.names = F, quote = F)
      invisible(sum_stats)
    }
  }
}

#----------------------------------------------------------------------------

file_names <- c('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/03_gwas_estimation/summarystats_f1_.txt',
                'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/03_gwas_estimation/summarystats_f2_.txt', 
                'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/03_gwas_estimation/summarystats_f3_.txt')

N <- readRDS('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/03_gwas_estimation/effective_sample_size_factors.RDS')
munge(files = file_names, hm3 = 'SNP/w_hm3.snplist', 
      N = N, 
      trait.names = c('f1', 'f2', 'f3'),
      parallel = T,
      cores = 3)
system(' mv *.sumstats.gz outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_ldsc_factors/munged')
system(' mv *log  outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_ldsc_factors/log')

#-------------------------------------------------------------------------------

# a <- fread('Summary_Stats/vuckovic-2020_lymphocytecount_build37_GCST90002388_buildGRCh37.tsv')
# b <- fread('Summary_Stats/vuckovic-2020_eosinophilcounts_build37_GCST90002381_buildGRCh37.tsv')
# c <- fread('Summary_Stats/vuckovic-2020_monocytecount_build37_GCST90002393_buildGRCh37.tsv')
# 
# 
# for(i in 1:3) {
#   prepare_munge(list(a,b,c)[[i]],
#                 rsID = 'variant_id', 
#                 the_effect_allele = 'effect_allele',
#                 the_non_effect_allele = 'other_allele', 
#                 the_Beta = 'beta',
#                 pvalue = 'p_value', 
#                 the_SE = 'standard_error', 
#                 the_chr = 'CHR_num',
#                 the_bp = 'base_pair_location',
#                 to_remove = c('VARIANT', 'GENPOS', 'MLOG10P', 'effect_allele_frequency', 'ALT_MINOR', 'MA_FREQ', 'R2', 'GWSIG', 'INFO'), 
#                 path = paste0('outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/10_ldsc_factors/ready_for_munge/', c('lymphocyte','eosinophil', 'monocyte')[i], '_counts.build37.txt'))
# }
# 
# 
# 
# files = c('outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/10_ldsc_factors/ready_for_munge/lymphocyte_counts.build37.txt',
#           'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/10_ldsc_factors/ready_for_munge/monocyte_counts.build37.txt',
#           'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/10_ldsc_factors/ready_for_munge/eosinophil_counts.build37.txt'
# )
# 
# 
# munge(files = files,
#       hm3 = 'SNP/w_hm3.snplist',
#       N=c(408112,408112,408112),
#       trait.names = c('lymphocyte_count', 'monocyte_counts', 'eosinophil_counts'),
#       parallel = T,
#       cores = 2)
# 
# system(' mv *.sumstats.gz outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/10_ldsc_factors/munged/')
# system(' mv *log  outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/10_ldsc_factors/munged/log/')
# 
# 




#------------ run ldsc for allergies -------------------------------------------


munged_files <- c('outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/10_ldsc_factors/munged/lymphocyte_count.sumstats.gz',
                  'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/10_ldsc_factors/munged/monocyte_counts.sumstats.gz',
                  'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/10_ldsc_factors/munged/eosinophil_counts.sumstats.gz',
                  
                  'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_ldsc_factors/munged/f1.sumstats.gz',
                  'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_ldsc_factors/munged/f2.sumstats.gz', 
                  'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_ldsc_factors/munged/f3.sumstats.gz', 
                  'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/02_sumstats_function/munged_files/t1d.sumstats.gz',
                  'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/02_sumstats_function/munged_files/crohn.sumstats.gz',
                  'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/02_sumstats_function/munged_files/uc.sumstats.gz',
                  'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/02_sumstats_function/munged_files/psc.sumstats.gz',
                  'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/02_sumstats_function/munged_files/jia.sumstats.gz', 
                  'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/02_sumstats_function/munged_files/sle.sumstats.gz', 
                  'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/02_sumstats_function/munged_files/ra.sumstats.gz', 
                  'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/02_sumstats_function/munged_files/asthma.sumstats.gz', 
                  'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/02_sumstats_function/munged_files/derma.sumstats.gz',
                  
                  'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output/psoriasis_sliz-2021.sumstats.gz',
                  'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output/allergies_ferreira-2017.sumstats.gz')



names = c('lympho','monocyte','eosinophil','f1', 'f2', 'f3', 't1d','crohn', 'uc', 'psc', 'jia', 'sle', 'ra',  'asthma', 'derma', 'psoriasis', 'allergies')


ldsc_output <- ldsc(traits = munged_files, 
                    sample.prev = c(NA,NA,NA,.5, .5, .5,.5, .5, .5, 0.1928, 0.2643, .5, .5, ((64538)/(64538 + 329321)), .5, .5),    
                    population.prev = c(NA,NA,NA,0.00001, 0.00001,0.00001 ,0.03, 0.001, 0.001, 0.000050, 0.000447, 0.000500,  0.0046, 0.035700, 0.20, 0.20,0.20), 
                    trait.names = names,
                    ld = "ldscores/eur_w_ld_chr",
                    wld= "ldscores/eur_w_ld_chr", stand = T)

saveRDS(ldsc_output, 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_ldsc_factors/ldsc_output.RDS')
ldsc_output <- readRDS('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_ldsc_factors/ldsc_output.RDS')

### plot -----------------------------------------------------------------------

rownames(ldsc_output$S_Stand) <- colnames(ldsc_output$S_Stand)


corrplot(ldsc_output$S_Stand, order = 'hclust',addCoef.col = 'black', method = 'square',type = 'upper', is.corr = F,
         tl.col='black', outline=F,
         number.cex= 1.5,
         tl.cex=1.5,
         cl.cex=1.5,
         tl.srt= 45,
         addgrid.col='grey',
         title= 'LDSC genetic correlation matrix of autoimmune diseases', mar=c(0,0,1,0)
)


rownames( ldsc_output$S_Stand)

#-------------------------------------------------------------------------------
#plot only f1

pdf('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_ldsc_factors/correlation_matrix_F1.pdf', height = 16, width = 16)

tp <- ldsc_output$S_Stand[c('f1', 'crohn', 'uc', 'psc'), c('f1', 'crohn', 'uc', 'psc') ]

corrplot(tp, order = 'hclust',addCoef.col = 'black', method = 'square',type = 'upper', is.corr = T,
         tl.col='black', outline=F,
         number.cex= 1.5,
         tl.cex=1.5,
         cl.cex=1.5,
         tl.srt= 45,
         addgrid.col='grey',
         title= 'LDSC genetic correlation matrix of autoimmune diseases', mar=c(0,0,1,0)
         
)
dev.off()

#-------------------------------------------------------------------------------
#plot only f2

pdf('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_ldsc_factors/correlation_matrix_F2.pdf', height = 16, width = 16)

tp <- ldsc_output$S_Stand[c('f2', 't1d', 'sle', 'ra', 'jia'), c('f2', 't1d', 'sle', 'ra', 'jia')]

corrplot(tp, order = 'hclust',addCoef.col = 'black', method = 'square',type = 'upper', is.corr = T,
         tl.col='black', outline=F,
         number.cex= 1.5,
         tl.cex=1.5,
         cl.cex=1.5,
         tl.srt= 45,
         addgrid.col='grey',
         title= 'LDSC genetic correlation matrix of autoimmune diseases', mar=c(0,0,1,0)
         
)
dev.off()




#-------------------------------------------------------------------------------
#plot only F3
pdf('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_ldsc_factors/correlation_matrix_F3.pdf', height = 16, width = 16)
tp <- ldsc_output$S_Stand[ c('f1', 'f2' ,'f3', 'psoriasis', 'allergies'), c( 'f1', 'f2' ,'f3',  'psoriasis', 'allergies') ]

corrplot(tp, order = 'hclust',addCoef.col = 'black', method = 'square',type = 'upper', is.corr = T,
         tl.col='black', outline=F,
         number.cex= 1.5,
         tl.cex=1.5,
         cl.cex=1.5,
         tl.srt= 45,
         addgrid.col='grey',
         title= 'LDSC genetic correlation matrix of autoimmune diseases', mar=c(0,0,1,0)
         
)
dev.off()

#-------------------------------------------------------------------------------


#plot cell couts
pdf('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_ldsc_factors/correlation_matrix_cell_counts.pdf', height = 16, width = 16)
tp <- ldsc_output$S_Stand[c('f1', 'f2', 'f3',  'lympho','monocyte','eosinophil'), c('f1', 'f2', 'f3',   'lympho','monocyte','eosinophil') ]
corrplot(tp, order = 'original',addCoef.col = 'black', method = 'square',type = 'upper', is.corr = T,
         tl.col='black', outline=F,
         number.cex= 1.5,
         tl.cex=1.5,
         cl.cex=1.5,
         tl.srt= 45,
         addgrid.col='grey',
         title= 'LDSC genetic correlation matrix of autoimmune diseases', mar=c(0,0,1,0)
)
dev.off()

library(ggplot2)
library(ggpubr)
library(RColorBrewer)



a1 <- tp[c(1,2,3),c(4,5,6)]
a2 <- tp[c(1,2,3),c(4,5,6)]
a3 <- tp[c(1,2,3),c(4,5,6)]

a1 <- as.data.frame(a1)
a2 <- as.data.frame(a2)
a3 <- as.data.frame(a3)

a1[,4] <- rownames(a1)
a2[,4] <- rownames(a2)
a3[,4] <- rownames(a3)

brewer.pal(6,'Paired')[c(2,4,6)]



lympho <- ggplot(a1, aes(x=V4, y=lympho ,fill=V4))+ geom_bar(stat='identity', width=0.5) +ylim(-0.1,+0.4)+geom_hline(yintercept=0,  color = "grey")+scale_fill_manual(values=brewer.pal(6,'Paired')[c(2,4,6)])+theme_classic()
mono <- ggplot(a1, aes(x=V4, y=monocyte, fill=V4))+ geom_bar(stat='identity', width=0.5) +ylim(-0.1,+0.4)+geom_hline(yintercept=0, color = "grey")+scale_fill_manual(values=brewer.pal(6,'Paired')[c(2,4,6)])+theme_classic()
eos <- ggplot(a1, aes(x=V4, y=eosinophil, fill=V4))+ geom_bar(stat='identity', width=0.5)  +ylim(-0.1,+0.4)+geom_hline(yintercept=0, color = "grey")+theme_classic()+ scale_fill_manual(values=brewer.pal(6,'Paired')[c(2,4,6)])+theme_classic()

pdf('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_ldsc_factors/barplot_cell_counts.pdf', height = 6, width = 15)
ggarrange(lympho, mono, eos, ncol = 3, nrow = 1, labels = c('lymphocyte counts', 'monocyte counts', 'eosinohil counts'))
dev.off()






























