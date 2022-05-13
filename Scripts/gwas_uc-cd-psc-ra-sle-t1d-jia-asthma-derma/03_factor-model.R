#  LD score regression on auotoimmunity GWAS
## version3  will be a replication of the version 2 that might be wrong in the 
## allele orientation 

#The tutorial and info on the package and how to run the code are here:  
# https://github.com/GenomicSEM/GenomicSEM/wiki/3.-Models-without-Individual-SNP-effects

# in this script I do the LD score regreeion of the summary stats munged in V3_step1 
# with only the studies we want to put into the factor model

#--- Load the libraries--------------------------------------
library(data.table)
library(GenomicSEM)
library(Matrix)
library(tidyr)
library(dplyr)
library(corrplot)
library(qgraph)
library(semPlot)
library(lavaan)


output2 <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/02_ldsc_alltraits/ldsc_output_all-traits.RDS')
rownames(output2$S_Stand) <- colnames(output2$S_Stand)


#----The model------------------------------------------------------------------
#remove allergies 

#Specify the Genomic confirmatory factor model




#ok very bello
aid_model <-'F1 =~ NA*crohn + uc  + psc  
F2 =~ NA*jia + sle + ra_eu + t1d 
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
#NO
# aid_model <-'F1 =~ NA*crohn + uc  + psc   
#              F2 =~ NA*jia + pbc + sle + ra_eu + t1d 
#              F3 =~ NA*asthma_2 + psoriasis + derma
#              
# F1~~F2
# F1 ~~ 1*F1
# F2 ~~ 1*F2
# F3 ~~ 1*F3
# derma~~a*derma
# a>0.001
# '

#run the model
aid_factor <-usermodel(output2, estimation = "DWLS", model = aid_model, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

#print the two factor model results
aid_factor

saveRDS(aid_factor, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/03_factor-model/model_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS')


#----Nicola's function for plotting SEM-----------------------------------------

semPlotModel_GSEM=function(gsem.object=GWISoutput , est.label="STD_All"){ 
  require(semPlot)
  object=gsem.object$results
  object$free=0
  numb=1:length(which(object$op!="~~"))
  object$free[which(object$op!="~~")]=numb
  varNames <- lavaanNames(object, type = "ov")
  factNames <- lavaanNames(object, type = "lv")
  factNames <- factNames[!factNames %in% varNames]
  n <- length(varNames)
  k <- length(factNames)
  if (is.null(object$label)) 
    object$label <- rep("", nrow(object))
  semModel <- new("semPlotModel")
  object$est <- object[,est.label]
  if (is.null(object$group)) 
    object$group <- ""
  semModel@Pars <- data.frame(label = object$label, lhs = ifelse(object$op == 
                                                                   "~" | object$op == "~1", object$rhs, object$lhs), edge = "--", 
                              rhs = ifelse(object$op == "~" | object$op == "~1", object$lhs, 
                                           object$rhs), est = object$est, std = NA, group = object$group, 
                              fixed = object$free==0, par = object$free, stringsAsFactors = FALSE)
  semModel@Pars$edge[object$op == "~~"] <- "<->"
  semModel@Pars$edge[object$op == "~*~"] <- "<->"
  semModel@Pars$edge[object$op == "~"] <- "~>"
  semModel@Pars$edge[object$op == "=~"] <- "->"
  semModel@Pars$edge[object$op == "~1"] <- "int"
  semModel@Pars$edge[grepl("\\|", object$op)] <- "|"
  semModel@Thresholds <- semModel@Pars[grepl("\\|", semModel@Pars$edge), 
                                       -(3:4)]
  semModel@Pars <- semModel@Pars[!object$op %in% c(":=", "<", 
                                                   ">", "==", "|", "<", ">"), ]
  semModel@Vars <- data.frame(name = c(varNames, factNames), 
                              manifest = c(varNames, factNames) %in% varNames, exogenous = NA, 
                              stringsAsFactors = FALSE)
  semModel@ObsCovs <- list()
  semModel@ImpCovs <- list()
  semModel@Computed <- FALSE
  semModel@Original <- list(object)
  return(semModel)
  
}

pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/03_factor-model/path_uc-cd-psc-ra-sle-t1d-jia-asthma-derm_factor.pdf', height = 7, width = 14)
semPaths(semPlotModel_GSEM(aid_factor) , what = 'est' , residuals = T, edge.color = 'black',
         sizeMan = 7, sizeLat = 7, label.cex=1, edge.label.cex = 0.8)
dev.off()




