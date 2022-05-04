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

#---- Crete a table that contains all the info on the GWAS----------------------

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
            'outputs/version3/01_output_prepare-sumstats/Munged-Sumstats/uc.sumstats.gz'
) 

trait.names <- c( 'allergies','ms_2','alzheimer_1', 'alzheimer_2', 'armfat', 'asthma_1', 'asthma_2', 
                  'celiac', 'crohn', 'jia', 'ms_1', 'pbc', 
                  'psc', 'ra', 'sle', 'thyro', 'uc')
sample.prev <-  c(.5, .5, .5, .5, NA, .5, (64538/(64538 + 239321)),
                        .5, .5, (3305/(3305 + 9196)), (9772 /(9772 + 16849)), .5, 
                        ( 2871 /(2871 + 12019)), .5, .5, .5, .5 )
                      
population.prev <-  c((0.20 ),(35.9/100000), (0.058), (0.058), (NA), (0.0357), (0.0357), 
                            (0.014), (100/100000), (44.7/100000), (35.9/100000), (10/100000),
                            (5/100000), (460/100000), (50/100000), (0.05), (30/100000) )

GWAS_info <- data.frame(trait.names, traits, sample.prev, population.prev)
row.names(GWAS_info) <- GWAS_info$trait.names

saveRDS(GWAS_info, file = 'outputs/version3/GWAS_info_table')

#----remove the studies I do not want to include in the next step---------------
#remove asthma_1 (demeanis is the smallest one compared to asthma_2 that is han e al)
#remove Alzh_1 that is the smallest (kunkle), leave Alzh_2 that is the bigger one (Wightman)
#remove m1 that is interacting with everything because of the low h2
#remove armfat that is not needed, it's only my internal control fo heritability (around 0.20)

GWAS_info_step3 <- GWAS_info[!(row.names(GWAS_info) %in% c('asthma_1', 'alzheimer_1', 'ms_1', 'armfat')),]


#----run ldsc function with only the studies we want----------------------------
ld <- "ldscores/eur_w_ld_chr"
wld <- "ldscores/eur_w_ld_chr"


ldsc_output_step3 <- ldsc(traits = GWAS_info_step3$traits, 
                          trait.names = GWAS_info_step3$trait.names, 
                          sample.prev = GWAS_info_step3$sample.prev, 
                          population.prev = GWAS_info_step3$population.prev,
                          ld = ld, 
                          wld =  wld ,
                          stand = T)

saveRDS(ldsc_output_step3, file='outputs/version3/03_output_prepare-sumstats/ldsc_output_step3')
ldsc_output_step3 <- readRDS('outputs/version3/03_output_prepare-sumstats/ldsc_output_step3')

rownames(ldsc_output_step3$S_Stand) <-  colnames(ldsc_output_step3$S_Stand) 

#visaul inspection to select the right model
corrplot(ldsc_output_step3$S_Stand, order = 'hclust', addCoef.col = 'black', type = 'upper')

qgraph(ldsc_output_step3$S_Stand, layout= 'spring', vsize=1.5, threshold=0.5, 
       node.width=2, diag=F, label.cex=1, color='black')

qgraph(ldsc_output_step3$S_Stand, layout= 'spring', vsize=1.5, threshold=0.6, 
       node.width=2, diag=F, label.cex=1, color='black')
#from qgraph a two factor model might be adequate, 
#factor 1 ms, psc, uc, and chrohn 
#factor 2 pbc, sle, jia 

#----The model------------------------------------------------------------------

#Specify the Genomic confirmatory factor model

aid_model <-'F1 =~ NA*crohn + uc  + psc 
             F2 =~ NA*jia + pbc + sle + ra
F1~~F2'

#run the model
aid_factor <-usermodel(ldsc_output_step3, estimation = "DWLS", model = aid_model, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

#print the two factor model results
aid_factor
#the statistic related to the fitting seem good sot the two factor model will be used to run the GWAS

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

semPaths(plot1, what = 'est' , residuals = T, groups = lista, 
         color = c('darksalmon', 'darkslategray3'), edge.color = 'black',
         sizeMan = 7, sizeLat = 7, label.cex=2, edge.label.cex = 1)





