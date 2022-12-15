library(GenomicSEM)
model <- readRDS('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/02_sumstats_function/ldsc_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS')



require(Matrix)
Ssmooth<-as.matrix((nearPD(model$S, corr = FALSE))$mat)

#run EFA with promax rotation and 2 factors using the factanal function in the stats package
require(stats)
EFA_1 <-factanal(covmat = Ssmooth, factors = 1, rotation = "promax")
EFA_2 <- factanal(covmat = Ssmooth, factors = 2, rotation = "promax")
EFA_3 <- factanal(covmat = Ssmooth, factors = 3, rotation = "promax")
EFA_4 <- factanal(covmat = Ssmooth, factors = 4, rotation = "promax")
EFA_5 <- factanal(covmat = Ssmooth, factors = 5, rotation = "promax")


#print the loadings
EFA_1$loadings
EFA_2$loadings
EFA_3$loadings
EFA_4$loadings
EFA_5$loadings


#plot the loadings
plotting <- data.frame('factors'= 1:5, 'cumulative var'=c(0.254, 0.426, 0.589, 0.670, 0.735 ))
plot(plotting$factors, plotting$cumulative.var, type = 'b', ylim = c(0,1), xlab = 'number of factors', main = 'Cumalative variance')
barplot(plotting$cumulative.var, names.arg = 1:5, ylim = c(0,1), main = 'Cumulative variance', xlab = 'number of factors')



# the more factors

ldsc_ok <- readRDS('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/02_sumstats_function/ldsc_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS')

aid_model_3 <-'F1 =~ NA*crohn + uc  + psc  
F2 =~ NA*jia + sle + ra + t1d 
F3 =~ NA*asthma + derma 

F1~~F2
F1~~F3
F2~~F3



F1 ~~ 1*F1
F2 ~~ 1*F2
F3 ~~ 1*F3
derma~~a*derma
a>0.001
'

aid_factor <-usermodel(ldsc_ok, estimation = "DWLS", model = aid_model_3, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
aid_factor

aid_model_4 <-'F1 =~ NA*crohn + uc  + psc  
F2 =~ NA*jia + sle + ra + t1d 
F3 =~ NA*asthma + derma 
F4 =~ NA*crohn + uc + jia

F1~~F2
F1~~F3
F2~~F3
F1~~F4
F2~~F4
F3~~F4

F4 ~~ 1*F4
F1 ~~ 1*F1
F2 ~~ 1*F2
F3 ~~ 1*F3
derma~~a*derma
a>0.001
'


#run the model
aid_factor <-usermodel(ldsc_ok, estimation = "DWLS", model = aid_model_4, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
aid_factor


