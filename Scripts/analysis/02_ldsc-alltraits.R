#  LD score regression on auotoimmunity GWAS
## version3  will be a replication of the version 2 that might be wrong in the 
## allele orientation 

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

traits <- c('outputs/version3/01_output_prepare-sumstats/Munged-Sumstats/allergies.sumstats.gz',  
            'outputs/version3/01_output_prepare-sumstats/Munged-Sumstats/ms_2.sumstats.gz' ,
            'outputs/version3/01_output_prepare-sumstats/Munged-Sumstats/alzheimer_1.sumstats.gz',
            'outputs/version3/01_output_prepare-sumstats/Munged-Sumstats/alzheimer_2.sumstats.gz',
            'outputs/version3/01_output_prepare-sumstats/Munged-Sumstats/armfat.sumstats.gz',
            'outputs/version3/01_output_prepare-sumstats/Munged-Sumstats/asthma_1.sumstats.gz',
            'outputs/version3/01_output_prepare-sumstats/Munged-Sumstats/asthma_2.sumstats.gz',
            
            'outputs/version3/01_output_prepare-sumstats/Munged-Sumstats/celiac.sumstats.gz',
            'outputs/version3/01_output_prepare-sumstats/Munged-Sumstats/crohn.sumstats.gz',
            'outputs/version3/01_output_prepare-sumstats/Munged-Sumstats/jia.sumstats.gz',
            'outputs/version3/01_output_prepare-sumstats/Munged-Sumstats/ms_1.sumstats.gz',
            'outputs/version3/01_output_prepare-sumstats/Munged-Sumstats/pbc.sumstats.gz',
            
            'outputs/version3/01_output_prepare-sumstats/Munged-Sumstats/psc.sumstats.gz', 
            'outputs/version3/01_output_prepare-sumstats/Munged-Sumstats/ra.sumstats.gz',
            'outputs/version3/01_output_prepare-sumstats/Munged-Sumstats/sle.sumstats.gz', 
            'outputs/version3/01_output_prepare-sumstats/Munged-Sumstats/thyro.sumstats.gz',
            'outputs/version3/01_output_prepare-sumstats/Munged-Sumstats/uc.sumstats.gz', 
            
            'outputs/version3/01_output_prepare-sumstats/Munged-Sumstats/t1d_chiou-2021.sumstats.gz' ,
            'outputs/version3/01_output_prepare-sumstats/Munged-Sumstats/derma_sliz-2021.sumstats.gz',
            'outputs/version3/01_output_prepare-sumstats/Munged-Sumstats/psoriasis_sliz-2021.sumstats.gz'
            ) 

trait.names <- c( 'allergies','ms_2','alzheimer_1', 'alzheimer_2', 'armfat', 'asthma_1', 'asthma_2', 
                 'celiac', 'crohn', 'jia', 'ms_1', 'pbc',  
                 'psc', 'ra', 'sle', 'thyro', 'uc', 
                 't1d', 'derma', 'psoriasis')
sample.prev <- round( c(.5, .5, .5, .5, NA, .5, (64538/(64538 + 239321)),
                        .5, .5, (3305/(3305 + 9196)), (9772 /(9772 + 16849)), .5, 
                        ( 2871 /(2871 + 12019)), .5, .5, .5, .5, 
                        .5, .5, .5)
                     ,2)
population.prev <-  round(c((0.20 ),(35.9/100000), (0.058), (0.058), (NA), (0.0357), (0.0357), 
                            (0.014), (100/100000), (44.7/100000), (35.9/100000), (10/100000),
                            (5/100000), (460/100000), (50/100000), (0.05), (30/100000), 
                            (0.095), (0.15), (0.02)
                            ),5)

ld <- "ldscores/eur_w_ld_chr"
wld <- "ldscores/eur_w_ld_chr"

# run the ld function
LDS_output <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, stand = T)

#save the output 
saveRDS(LDS_output, file = 'outputs/version3/02_output_ldsc-all-traits/LDS_output_complete')
output2 <- readRDS('outputs/version3/02_output_ldsc-all-traits/LDS_output_complete')

#---- heritability -------------------------------------------------------------
cbind(colnames(output2$S_Stand), (diag(output2$S)) )

#-----plot the final matrix ----------------------------------------------------

rownames(output2$S_Stand) <- colnames(output2$S_Stand)


pdf(file = 'outputs/version3/02_output_ldsc-all-traits/correlation-matrix_complete.pdf', height = 14, width = 14 )
corrplot( output2$S_Stand, order = 'hclust', addCoef.col = 'black', is.corr = F)
dev.off()

qgraph(output2$S_Stand,threshold=0.4,layout="spring")



#---- remnove m1 that is interacting with everything because of the low h2
t1 <- output2$S_Stand[-11, -11]


pdf(file = 'Graphs/V3/Correlation-matrix_noMS_1.pdf', height = 14, width = 14 )
corrplot(t1, order = 'hclust', addCoef.col = 'black')
dev.off()


#remove autocrrealtion for qraph
?qgraph
t_no_auto <- t1
diag(t_no_auto) <- rep(0, 16)

pdf(file = 'Graphs/V3/Network_noMS_1.pdf', height = 14, width = 14 )
qgraph(t_no_auto,threshold=0.4,layout="spring")
dev.off()









