#  LD score regression on auotoimmunity GWAS
## this script is for running ldsc on all the traits prepared in script 01_prepare-sumstats.R
## 

#The tutorial and info on the package and how to run the code are here:  
# https://github.com/GenomicSEM/GenomicSEM/wiki/3.-Models-without-Individual-SNP-effects

# in this script I do the LD score regreeion of the summary stats munged in V3_step1
#--- Load the libraries--------------------------------------
library(data.table)
library(GenomicSEM)
library(stats)
library(tidyr)
library(dplyr)
library(corrplot)
library(qgraph)

#---- LD score regression function----------------------------------------------

traits <- c('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output/allergies_ferreira-2017.sumstats.gz',  
            'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output/ms_andlauer.sumstats.gz',
            'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output/alzheimer_kunkle-2019.sumstats.gz',
            'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output/armfat.sumstats.gz',
            
            'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output/asthma_han-2020.sumstats.gz',
            'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output/celiac_dubois-2010.sumstats.gz',
            'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output/crohn_delange-2017.sumstats.gz',
            'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output/jia_lopezisac-2020.sumstats.gz',
            'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output/pbc_cordell-2015.sumstats.gz',
            
            'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output/psc_ji-2016.sumstats.gz', 
            'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output/ra_okada-2014_only-eu.sumstats.gz',
            'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output/sle_bentham-2015.sumstats.gz', 
            'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output/thyro_saevarsdottir-2020.sumstats.gz',
            'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output/uc_delange-2017.sumstats.gz', 
            
            'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output/t1d_chiou-2021.sumstats.gz' ,
            'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output/derma_sliz-2021.sumstats.gz',
            'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output/psoriasis_sliz-2021.sumstats.gz', 
            'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output/vitiligo_jin-2016.sumstats.gz'
            
            ) 

trait.names <- c( 'allergies','ms', 'alzheimer', 'armfat',  
                  'asthma', 'celiac', 'crohn', 'jia', 'pbc',  
                 'psc', 'ra_eu', 'sle', 'thyro', 'uc', 
                 't1d', 'derma', 'psoriasis', 'vitiligo')
sample.prev <- round( c(.5, .5, .5, NA,
                        (64538/(64538 + 239321)),.5, .5, (3305/(3305 + 9196)), .5, 
                        ( 2871 /(2871 + 12019)), .5, .5, .5, .5, 
                        .5, .5, .5, .5)
                     ,2)
population.prev <-  round(c((0.20 ),(35.9/100000), (0.058), (NA),
                            (0.0357),  (0.014), (100/100000), (44.7/100000), (10/100000),
                            (5/100000), (460/100000), (50/100000), (0.05), (30/100000), 
                            (0.095), (0.15), (0.02), (0.002)
                            ),5)

ld <- "ldscores/eur_w_ld_chr"
wld <- "ldscores/eur_w_ld_chr"

# run the ld function
LDS_output <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, stand = T)
system('mv *.log outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/02_ldsc_alltraits/')
#save the output 
saveRDS(LDS_output, file = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/02_ldsc_alltraits/ldsc_output_all-traits.RDS')
output2 <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/02_ldsc_alltraits/ldsc_output_all-traits.RDS')

#---- heritability -------------------------------------------------------------
cbind(colnames(output2$S_Stand), (diag(output2$S)) )

#-----plot the final matrix ----------------------------------------------------

rownames(output2$S_Stand) <- colnames(output2$S_Stand)


pdf(file = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/02_ldsc_alltraits/correlation-matrix_complete.pdf', height = 14, width = 14 )
corrplot( output2$S_Stand, order = 'hclust', addCoef.col = 'black', is.corr = T, type = 'upper')
dev.off()

pdf(file = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/02_ldsc_alltraits/network_complete.pdf', height = 14, width = 14 )
qgraph(output2$S_Stand,threshold=0.4,layout="spring", diag=F)
dev.off()




