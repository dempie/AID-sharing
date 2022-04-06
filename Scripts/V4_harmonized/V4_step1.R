#  LD score regression on auotoimmunity GWAS
## Version4  will be a replication of the version 3 with harmonized GWAS summary stats


#The tutorial and info on the package and how to run the code are here:  
# https://github.com/GenomicSEM/GenomicSEM/wiki/3.-Models-without-Individual-SNP-effects

# i
#--- Load the libraries--------------------------------------
library(data.table)
library(GenomicSEM)
library(tidyr)
library(dplyr)
library(corrplot)
library(ggplot2)

#---- Create a function to prepare the summary stats ---------------------------
## create a function to prepare the GWAS for the munging depends on dplyr and data.table
##the function changes the name of the columns as required by the genomicSEM package
#and saves the files in the provided path (if the path is provided)

prepare_munge_2 <- function(sum_stats, rsID, effect_allele, non_effect_allele, pvalue, effect_size, to_remove=NA, path = NA){
  #an error if arguments are not provided 
  if (missing(sum_stats) | missing(rsID) | missing(effect_allele) | missing(non_effect_allele) |missing(pvalue) | missing(effect_size) ) {
    stop( 'At least one argument is missing')
    
  } else {
    
    require(dplyr)
    require(data.table)
    sum_stats <- sum_stats  %>% rename(c(SNP = rsID, A1 = effect_allele, A2 = non_effect_allele, p = pvalue, effect = effect_size))
    #conditional remove
    if(is.na(to_remove[1])){
      #do nothing 
    } else {
      sum_stats <- select(sum_stats,-(to_remove))
    }
    
    #save the file if a path is provided
    if(is.na(path)){
      return(sum_stats)
    } else {
      fwrite(sum_stats, path, sep = '\t', col.names = T, row.names = F, quote = F)
      return(sum_stats)
    }
  }
}
#------end of function ---------------------------------------------------------

#Divide into 2 groups 
#---- Group1 load the datasets -------------------------------------------------

cron_harmonized <- fread('Summary_Stats/harmonized/crohn_delange-2017_harmonized_28067908-GCST004132-EFO_0000384.h.tsv.gz', data.table = F)  
uc_harmonized <- fread('Summary_Stats/harmonized/uc_delange-2017_harmonized_28067908-GCST004133-EFO_0000729.h.tsv.gz', data.table = F)
allergies_harmonized <- fread('Summary_Stats/harmonized/allergies_ferreira-2017_harmonized_29083406-GCST005038-EFO_0003785.h.tsv.gz', data.table = F)
asthma_harmonized <- fread('Summary_Stats/harmonized/asthma_demeanis-2018_harmonized_29273806-GCST006862-EFO_0000270.h.tsv.gz', data.table = F)
celiac_harmonized <- fread('Summary_Stats/harmonized/celiac_dubois-2020_harmonized_20190752-GCST000612-EFO_0001060.h.tsv.gz' , data.table = F)

#---- crohn GWAS ---------------------------------------------------------------

head(cron_harmonized)

cron_harmonized <- select(cron_harmonized, c('variant_id', 'hm_effect_allele' , 'hm_other_allele', 'hm_beta', 'p_value' ))
prepare_munge_2(cron_harmonized,
                rsID = 'variant_id', 
                effect_allele = "hm_effect_allele", 
                non_effect_allele = 'hm_other_allele', 
                effect = 'hm_beta', 
                pvalue = 'p_value', 
                path = 'Outputs/V4_harmonized/Ready_for_munge/cron_harmonised.txt'
                              )


#---- uc GWAS ------------------------------------------------------------------
 
head(uc_harmonized)
uc_harmonized <- select(uc_harmonized, c('variant_id', 'hm_effect_allele' , 'hm_other_allele', 'hm_beta', 'p_value' ))
prepare_munge_2(uc_harmonized,
                rsID = 'variant_id', 
                effect_allele = "hm_effect_allele", 
                non_effect_allele = 'hm_other_allele', 
                effect = 'hm_beta', 
                pvalue = 'p_value', 
                path = 'Outputs/V4_harmonized/Ready_for_munge/uc_harmonised.txt')


#---- allergies GWAS -----------------------------------------------------------

head(allergies_harmonized)
allergies_harmonized <- select(allergies_harmonized, c('variant_id', 'hm_effect_allele' , 'hm_other_allele', 'hm_beta', 'p_value' ))
prepare_munge_2(allergies_harmonized,
                rsID = 'variant_id', 
                effect_allele = "hm_effect_allele", 
                non_effect_allele = 'hm_other_allele', 
                effect = 'hm_beta', 
                pvalue = 'p_value', 
                path = 'Outputs/V4_harmonized/Ready_for_munge/allergies_harmonised.txt')

#---- asthma GWAS --------------------------------------------------------------

head(asthma_harmonized)
asthma_harmonized <- select(asthma_harmonized, c('variant_id', 'hm_effect_allele' , 'hm_other_allele', 'hm_beta', 'p_value' ))

prepare_munge_2(asthma_harmonized,
                rsID = 'variant_id', 
                effect_allele = "hm_effect_allele", 
                non_effect_allele = 'hm_other_allele', 
                effect = 'hm_beta', 
                pvalue = 'p_value', 
                path = 'Outputs/V4_harmonized/Ready_for_munge/asthma_harmonised.txt')

#---- celiac GWAS --------------------------------------------------------------

head(celiac_harmonized)
celiac_harmonized <- select(celiac_harmonized, c('variant_id', 'hm_effect_allele' , 'hm_other_allele', 'hm_beta', 'p_value' ))

prepare_munge_2(celiac_harmonized,
                rsID = 'variant_id', 
                effect_allele = "hm_effect_allele", 
                non_effect_allele = 'hm_other_allele', 
                effect = 'hm_beta', 
                pvalue = 'p_value', 
                path = 'Outputs/V4_harmonized/Ready_for_munge/celiac_harmonised.txt')


#---- Group1 munge -------------------------------------------------------------

vector_files <- c('Outputs/V4_harmonized/Ready_for_munge/cron_harmonised.txt', 
                  'Outputs/V4_harmonized/Ready_for_munge/uc_harmonised.txt', 
                  'Outputs/V4_harmonized/Ready_for_munge/allergies_harmonised.txt', 
                  'Outputs/V4_harmonized/Ready_for_munge/asthma_harmonised.txt', 
                  'Outputs/V4_harmonized/Ready_for_munge/celiac_harmonised.txt')

munge(vector_files, 
      trait.names = c('crohn', 'uc', 'allergies', 'asthma', 'celiac'), 
      hm3 = 'SNP/w_hm3.snplist', 
      N = c(40266, 45975, 360838, 142486,  15283)
      )


#---- Group2 load the datasets--------------------------------------------------
jia_harmonized <- fread('Summary_Stats/harmonized/jia_lopezisac-2020_harmonized_33106285-GCST90010715-EFO_0002609.h.tsv.gz' ,data.table = F)
ms_harmonized <- fread('Summary_Stats/harmonized/ms_ims-2011_harmonized_21833088-GCST001198-EFO_0003885.h.tsv.gz', data.table = F)
pbc_harmonized <- fread('Summary_Stats/harmonized/pbc_cordel-2015_harmonized_26394269-GCST003129-EFO_1001486.h.tsv.gz' , data.table = F)
psc_harmonized <- fread( 'Summary_Stats/harmonized/psc_ji-2016_harmonized_27992413-GCST004030-EFO_0004268.h.tsv.gz', data.table = F)
ra_harmonized <-  fread('Summary_Stats/harmonized/ra_okada-2014_harmonized_24390342-GCST002318-EFO_0000685.h.tsv.gz' , data.table = F)
sle_harmoninez <- fread('Summary_Stats/harmonized/sle_bentham-2015_harmonized_26502338-GCST003156-EFO_0002690.h.tsv.gz', data.table = F)


#---- jia GWAS -----------------------------------------------------------------

head(jia_harmonized)

#no effect reported, everythin is na

#---- ms GWAS ------------------------------------------------------------------

head(ms_harmonized)
ms_harmonized <- select(ms_harmonized, c('variant_id', 'hm_effect_allele' , 'hm_other_allele', 'hm_beta', 'p_value' , 'standard_error'))

prepare_munge_2(ms_harmonized,
                 rsID = 'variant_id', 
                 effect_allele = "hm_effect_allele", 
                 non_effect_allele = 'hm_other_allele', 
                 effect = 'hm_beta', 
                 pvalue = 'p_value', 
                 path = 'Outputs/V4_harmonized/Ready_for_munge/ms_harmonised.txt')
#---- pbc GWAS -----------------------------------------------------------------

head(pbc_harmonized)
pbc_harmonized <- select(pbc_harmonized, c('variant_id', 'hm_effect_allele' , 'hm_other_allele', 'hm_beta', 'p_value' , 'standard_error'))

prepare_munge_2(pbc_harmonized,
                 rsID = 'variant_id', 
                 effect_allele = "hm_effect_allele", 
                 non_effect_allele = 'hm_other_allele', 
                 effect = 'hm_beta', 
                 pvalue = 'p_value', 
                 path = 'Outputs/V4_harmonized/Ready_for_munge/pbc_harmonised.txt')

#--- psc GWAS ------------------------------------------------------------------
head(psc_harmonized)

psc_harmonized <- select(psc_harmonized, c('variant_id', 'hm_effect_allele' , 'hm_other_allele', 'hm_odds_ratio', 'p_value', 'standard_error'))

prepare_munge_2(psc_harmonized,
                 rsID = 'variant_id', 
                 effect_allele = "hm_effect_allele", 
                 non_effect_allele = 'hm_other_allele', 
                 effect = 'hm_odds_ratio', 
                 pvalue = 'p_value', 
                 path = 'Outputs/V4_harmonized/Ready_for_munge/psc_harmonised.txt')
#---- ra GWAS ------------------------------------------------------------------

head(ra_harmonized)
ra_harmonized <- select(ra_harmonized, c('variant_id', 'hm_effect_allele' , 'hm_other_allele', 'hm_odds_ratio', 'p_value' ))

prepare_munge_2(ra_harmonized,
                 rsID = 'variant_id', 
                 effect_allele = "hm_effect_allele", 
                 non_effect_allele = 'hm_other_allele', 
                 effect = 'hm_odds_ratio', 
                 pvalue = 'p_value', 
                 path = 'Outputs/V4_harmonized/Ready_for_munge/ra_harmonised.txt')
#---- sle GWAS -----------------------------------------------------------------

head(sle_harmoninez)
sle_harmonized <- select( sle_harmoninez, c('variant_id', 'hm_effect_allele' , 'hm_other_allele', 'hm_beta', 'p_value' ))

prepare_munge_2(sle_harmonized ,
                 rsID = 'variant_id', 
                 effect_allele = "hm_effect_allele", 
                 non_effect_allele = 'hm_other_allele', 
                 effect = 'hm_beta', 
                 pvalue = 'p_value', 
                 path = 'Outputs/V4_harmonized/Ready_for_munge/sle_harmonized.txt')

#---- Group 2 munge-------------------------------------------------------------

vector_files <- c('Outputs/V4_harmonized/Ready_for_munge/ms_harmonised.txt', 
                  'Outputs/V4_harmonized/Ready_for_munge/pbc_harmonised.txt', 
                  'Outputs/V4_harmonized/Ready_for_munge/psc_harmonised.txt', 
                  'Outputs/V4_harmonized/Ready_for_munge/ra_harmonised.txt', 
                  'Outputs/V4_harmonized/Ready_for_munge/sle_harmonized.txt')

munge(vector_files, 
      trait.names = c('ms', 'pbc', 'psc', 'ra', 'sle'), 
      hm3 = 'SNP/w_hm3.snplist', 
      N = c(26621, 13239, 14890, 103638, 14257)
)


#------ run LDSR---------------------------------------------------------------

traits <- c( 'Outputs/V4_harmonized/Munged/crohn.sumstats.gz', 
             'Outputs/V4_harmonized/Munged/uc.sumstats.gz', 
              
             'Outputs/V4_harmonized/Munged/asthma.sumstats.gz', 
             'Outputs/V4_harmonized/Munged/celiac.sumstats.gz',
             
             'Outputs/V4_harmonized/Munged/ms.sumstats.gz', 
             'Outputs/V4_harmonized/Munged/pbc.sumstats.gz', 
             'Outputs/V4_harmonized/Munged/psc.sumstats.gz', 
             'Outputs/V4_harmonized/Munged/ra.sumstats.gz', 
             'Outputs/V4_harmonized/Munged/sle.sumstats.gz'
             ) 

trait.names <- c('crohn', 'uc', 'asthma', 'celiac',
                 'ms', 'pbc', 'psc', 'ra', 'sle' )

sample.prev <- round(c( (12194/28072),(12366/33609), (19954/107715), (4533/10750),
                       (9772/16849), (2764/10475), (2871/12019), (19234/61565), (5201/9066) )
                       ,2)
#


population.prev <-  round(c( (100/100000), (30/100000), (0.0357), (0.014), 
                             (35.9/100000), (10/100000), (5/100000), (460/100000), (50/100000) )
                          ,5)

ld <- "ldscores/eur_w_ld_chr"
wld <- "ldscores/eur_w_ld_chr"

# run the ld function
LDS_output <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, stand = T)

#save the output 
saveRDS(LDS_output, file = 'Outputs/V4_harmonized/LDS_output')
output <- readRDS('Outputs/V4_harmonized/LDS_output')

diag(output$S)
rownames(output$S_Stand) <- colnames(output$S_Stand)
pdf('test.pdf', height = 14, width = 14)
corrplot(output$S_Stand, order = 'hclust', addCoef.col = 'black', is.corr = F) 
mtext('Harmonized Summary Statistics', at=5, line=3, cex=1)
dev.off()
#----- plotting ----------------------------------------------------------------






