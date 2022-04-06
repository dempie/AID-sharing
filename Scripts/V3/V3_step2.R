#  LD score regression on auotoimmunity GWAS
## Version3  will be a replication of the version 2 that might be wrong in the 
## allele orientation 

#The tutorial and info on the package and how to run the code are here:  
# https://github.com/GenomicSEM/GenomicSEM/wiki/3.-Models-without-Individual-SNP-effects

# in this script I do the LD score regreeion of the summary stats munged in V£_step1
#--- Load the libraries--------------------------------------
library(data.table)
library(GenomicSEM)
library(tidyr)
library(dplyr)
library(corrplot)

#---- LD score regression function----------------------------------------------

traits <- c(  
            'Outputs/Version3/Munged-Sumstats/alzheimer_1.sumstats.gz',
            'Outputs/Version3/Munged-Sumstats/alzheimer_2.sumstats.gz',
            'Outputs/Version3/Munged-Sumstats/armfat.sumstats.gz',
            'Outputs/Version3/Munged-Sumstats/asthma_1.sumstats.gz',
            'Outputs/Version3/Munged-Sumstats/asthma_2.sumstats.gz',
            
            'Outputs/Version3/Munged-Sumstats/celiac.sumstats.gz',
            'Outputs/Version3/Munged-Sumstats/crohn.sumstats.gz',
            'Outputs/Version3/Munged-Sumstats/jia.sumstats.gz',
            'Outputs/Version3/Munged-Sumstats/ms.sumstats.gz',
            'Outputs/Version3/Munged-Sumstats/pbc.sumstats.gz',
            
            'Outputs/Version3/Munged-Sumstats/psc.sumstats.gz', 
            'Outputs/Version3/Munged-Sumstats/ra.sumstats.gz',
            'Outputs/Version3/Munged-Sumstats/sle.sumstats.gz', 
            'Outputs/Version3/Munged-Sumstats/thyro.sumstats.gz',
            'Outputs/Version3/Munged-Sumstats/uc.sumstats.gz'
            ) 

trait.names <- c('alzheimer_1', 'alzheimer_2', 'armfat', 'asthma_1', 'asthma_2', 
                 'celiac', 'crohn', 'jia', 'ms', 'pbc', 
                 'psc', 'ra', 'sle', 'thyro', 'uc')
sample.prev <- round(c((21982/41944), (86531/676486), (NA), (19954/107715), (64538/239321), 
                       (4533/10750), (12194/28072), (3305/9196), (9772/16849), (2764/10475),
                       (2871/12019), (19234/61565), (5201/9066), (30324/725172), (12366/33609)
                       ),2)
population.prev <-  round(c((0.058), (0.058), (NA), (0.0357), (0.0357), 
                            (0.014), (100/100000), (44.7/100000), (35.9/100000), (10/100000),
                            (5/100000), (460/100000), (50/100000), (0.05), (30/100000)
                            ),5)

ld <- "ldscores/eur_w_ld_chr"
wld <- "ldscores/eur_w_ld_chr"

# run the ld function
LDS_output <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, stand = T)

#save the output 
saveRDS(LDS_output, file = 'Outputs/Version3/LDS_output')
output2 <- readRDS('Outputs/Version3/LDS_output')

diag(output2$S)
rownames(output2$S_Stand) <- colnames(output2$S_Stand)

pdf('test.2.pdf', height = 14, width = 14)
corrplot(output2$S_Stand, order = 'hclust', addCoef.col = 'black', is.corr = F)
mtext('Non harmonized Summary Statistics', at=5, line=3, cex=1)

dev.off()

install.packages('qgraph')
library(qgraph)


qgraph(output2$S_Stand,threshold=0.5,layout="spring")

summary(lm(iris$Sepal.Length~iris$Petal.Length))






