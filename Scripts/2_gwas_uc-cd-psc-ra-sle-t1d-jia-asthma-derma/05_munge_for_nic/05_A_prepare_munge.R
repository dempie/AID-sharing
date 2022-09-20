# in this script all the necessary steps for  preparation for the mungin for the colocalizzation will be performed 
library(data.table)
library(dplyr)
library(tidyr)




#all the sumstats in the gwas have Beta, except for psc and asthma that has OR, 
#to that ones I just add the beta column

#-------------------------------------------------------------------------------

#---- Create a function to prepare the summary stats ---------------------------
## create a function to prepare the GWAS for the munging depends on dplyr and data.table
##the function changes the name of the columns as required by the library(MungeSumstats) package
#and saves the files in the provided path (if the path is provided)
#the rename function is from dplyr
prepare_munge <- function(sum_stats, rsID, the_effect_allele, the_non_effect_allele, pvalue, the_OR=NA, the_Beta=NA, the_SE=NA, the_chr=NA, the_bp=NA, to_remove=NA, path = NA){
  #an error if arguments are not provided
  if (missing(sum_stats) | missing(rsID) | missing(the_effect_allele) | missing(the_non_effect_allele) |missing(pvalue) ) {
    
    stop( 'At least one argument is missing')
    
  } else {
    
    require(dplyr)
    require(data.table)
    sum_stats <- sum_stats  %>% rename(c(SNP = all_of(rsID), EFFECT_ALLELE  = all_of(the_effect_allele), NON_EFFECT_ALLELE = all_of(the_non_effect_allele), p = all_of(pvalue) ))
    sum_stats$SNP <- tolower(sum_stats$SNP)
    sum_stats$p <- as.numeric(sum_stats$p)
    sum_stats$EFFECT_ALLELE <- toupper(as.character( sum_stats$EFFECT_ALLELE))
    sum_stats$NON_EFFECT_ALLELE <- toupper(as.character( sum_stats$NON_EFFECT_ALLELE))
    
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

#-------f1 prepare GWAS -------------------------------------------------------
f1 <- fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/03_gwas_estimation/summarystats_f1_.txt', data.table = F)
head(f1)

f1 <- prepare_munge(f1,
                    rsID = 'SNP',
                    the_effect_allele = 'A1',
                    the_non_effect_allele = 'A2',
                    pvalue = 'P',
                    the_Beta = 'Beta',
                    the_SE = 'SE',
                    the_chr = 'CHR',
                    path = 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/ready_for_munge/f1_ready_for_munge.txt'
)

head(f1)
#-------f2 prepare GWAS -------------------------------------------------------
f2 <- fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/03_gwas_estimation/summarystats_f2_.txt', data.table = F)
head(f2)

f2 <- prepare_munge(f2,
                    rsID = 'SNP',
                    the_effect_allele = 'A1',
                    the_non_effect_allele = 'A2',
                    pvalue = 'P',
                    the_Beta = 'Beta',
                    the_SE = 'SE',
                    path = 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/ready_for_munge/f2_ready_for_munge.txt'
)
head(f2)

#-------f3 prepare GWAS -------------------------------------------------------
f3 <- fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/03_gwas_estimation/summarystats_f3_.txt', data.table = F)
head(f3)

f3 <- prepare_munge(f3,
                    rsID = 'SNP',
                    the_effect_allele = 'A1',
                    the_non_effect_allele = 'A2',
                    pvalue = 'P',
                    the_Beta = 'Beta',
                    the_SE = 'SE',
                    path = 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/ready_for_munge/f3_ready_for_munge.txt'
)

#----------uc prepare GWAS -----------------------------------------------------
uc <- fread('Summary_Stats/delange-2017_uc_build37_45975_20161107.txt', data.table = F)
uc <- uc %>% separate(col = 'MarkerName', sep = '_|_', remove = F, convert = T, extra = 'warn', into = c('SNP', 'Extra1', 'Extra2') ) %>% select(-c(Extra1, Extra2))  %>% separate(col='SNP', sep= ':', convert=F,  into = c('CHR', 'BP'), remove = F)
head(uc)    
uc <- select(uc,-c('SNP'))

prepare_munge(uc,
              rsID = 'MarkerName',
              the_effect_allele = 'Allele2',
              the_non_effect_allele = 'Allele1',
              pvalue = 'P.value',
              the_Beta = 'Effect',
              the_SE = 'StdErr',
              to_remove = c( 'Min_single_cohort_pval', 'Pval_GWAS3', 'Pval_IIBDGC', 'Pval_IBDseq', 'HetPVal', 'HetDf', 'HetChiSq', 'HetISq'  ),
              path = 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/ready_for_munge/uc_ready_for_munge_build37.txt'
)
head(uc)

#-----crohn prepare GWAS -------------------------------------------------------
cd <- fread('Summary_Stats/delange-2017_cd_build37_40266_20161107.txt', data.table = F)
cd <- cd %>% separate(col = 'MarkerName', sep = '_|_', remove = F, convert = T, extra = 'warn', into = c('SNP', 'Extra1', 'Extra2') ) %>% select(-c(Extra1, Extra2))  %>% separate(col='SNP', sep= ':', convert=F,  into = c('CHR', 'BP'), remove = F)
head(cd)    
cd <- select(cd,-c('SNP'))


prepare_munge(cd,
              rsID = 'MarkerName',
              the_effect_allele = 'Allele2',
              the_non_effect_allele = 'Allele1',
              pvalue = 'P.value',
              the_Beta = 'Effect',
              the_SE = 'StdErr',
              to_remove = c( 'Min_single_cohort_pval', 'Pval_GWAS3', 'Pval_IIBDGC', 'Pval_IBDseq', 'HetPVal', 'HetDf', 'HetChiSq', 'HetISq'  ),
              path = 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/ready_for_munge/cd_ready_for_munge_build37.txt'
)


#-----psc prepare GWAS----------------------------------------------------------
psc <- fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_summary_stats/psc_ji-2016.txt', data.table = F)
head(psc)


psc$beta <- log(psc$OR)
psc <- rename(psc, the_effect_in_OR=OR )

prepare_munge(psc,
              rsID = 'SNP',
              the_effect_allele = 'A1',
              the_non_effect_allele = 'A2',
              pvalue = 'P',
              the_Beta = 'beta',
              the_SE = 'SE',
              path = 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/ready_for_munge/psc_ready_for_munge_build37.txt'
)


#----jia prepare GWAS ----------------------------------------------------------
jia <- fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_summary_stats/jia_beta_lopezisac-2020.txt', data.table = F)
head(jia)

prepare_munge(jia,
              rsID = 'SNP',
              the_effect_allele = 'A1',
              the_non_effect_allele = 'A2',
              pvalue = 'p',
              the_Beta = 'Beta',
              the_SE = 'SE',
              path = 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/ready_for_munge/jia_ready_for_munge_build37.txt'
)

#-----sle prepare gwas ---------------------------------------------------------
sle <- fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_summary_stats/sle_beta_bentham-2015.txt', data.table = F)
head(sle)

prepare_munge(sle,
              rsID = 'SNP',
              the_effect_allele = 'A1',
              the_non_effect_allele = 'A2',
              pvalue = 'p',
              the_Beta = 'Beta',
              the_SE = 'SE',
              path = 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/ready_for_munge/sle_ready_for_munge_build37.txt'
)

#-----t1d prepare gwas ---------------------------------------------------------
t1d <- fread('Summary_Stats/chiou-2021_t1d_build38_GCST90014023_buildGRCh38.tsv', data.table=F)
head(t1d)
prepare_munge(t1d,
              rsID = 'variant_id',
              the_effect_allele = 'effect_allele',
              the_non_effect_allele = 'other_allele',
              pvalue = 'p_value',
              the_Beta = 'beta',
              the_SE = 'standard_error',
              the_chr = 'chromosome',
              the_bp = 'base_pair_location',
              path = 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/ready_for_munge/t1d_ready_for_mung_build38.txt'
)

#-----asthma prepare GWAS-------------------------------------------------------
asthma <- fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_summary_stats/asthma_han-2020.txt', data.table = F)
asthma$beta <- log(asthma$OR)
asthma <- rename(asthma, the_effect_in_OR=OR)

prepare_munge(asthma,
              rsID = 'SNP',
              the_effect_allele = 'A1',
              the_non_effect_allele = 'A2',
              pvalue = 'p',
              the_Beta = 'beta',
              the_SE = 'SE',
              path = 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/ready_for_munge/asthma_ready_for_munge_build37.txt'
)


#----------ra preare GWAS ------------------------------------------------------
ra <- fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_summary_stats/ra_eu_okada-2014_chr_bp.txt', data.table = F)

prepare_munge(ra,
              rsID = 'SNP',
              the_effect_allele = 'A1',
              the_non_effect_allele = 'A2',
              pvalue = 'p',
              the_Beta = 'Beta',
              the_SE = 'SE',
              path = 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/ready_for_munge/ra_ready_for_munge_build37.txt'
)


#-------------derma prepare GWAS
derma <- fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_summary_stats/derma_sliz-2021.txt', data.table = F)
head(derma)
prepare_munge(derma,
              rsID = 'SNP',
              the_effect_allele = 'A1',
              the_non_effect_allele = 'A2',
              pvalue = 'p',
              the_Beta = 'Beta',
              the_SE = 'SE',
              path = 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/ready_for_munge/derma_ready_for_munge_build38.txt'
)

#---------- MungeSumstats ------------------------------------------------------
#check the build
library(MungeSumstats)
a1 <- list.files('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/ready_for_munge/')
my_paths <- as.list(paste0('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/ready_for_munge/', a1))
my_paths <- my_paths [-10] #remove README
a1 <- a1[-10] #remove README
names(my_paths) <- a1 #give names

my_paths <- my_paths[-c(2,12)] #uc and cd do not have rsID and genome build cannot be inferred by the function, I checked in the paper and it is build37

builds <- get_genome_builds(my_paths) #check the genome build, only t1d and derma are build 38

#-------------------------------------------------------------------------------



