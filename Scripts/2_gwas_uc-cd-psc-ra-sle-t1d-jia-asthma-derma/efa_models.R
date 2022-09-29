library(GenomicSEM)
model <- readRDS('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/02_sumstats_function/ldsc_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS')
model


require(Matrix)
Ssmooth<-as.matrix((nearPD(model$S, corr = FALSE))$mat)

#run EFA with promax rotation and 2 factors using the factanal function in the stats package
require(stats)
EFA<-factanal(covmat = Ssmooth, factors = 3, rotation = "promax")

#print the loadings
EFA$loadings