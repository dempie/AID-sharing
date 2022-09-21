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
library(RColorBrewer)
library(gridExtra)

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

names_plot = c('Type 1 Diabetes','Crohn\'s disease', 'Ulcerative colitis', 
               'Primary sclerosing cholangitis ', 'Juvenile idiopathic arthritis', 'Systemic lupus erythematosus', 
               'Rheumatoid arthritis',  'Asthma', 'Atopic dermatitis')

#----- load the correlation matrix ---------------------------------------------

ldsc_ok <- readRDS('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/02_sumstats_function/ldsc_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS')


#---- Three factor model ---------------------------------------------------------
aid_model <-'F1 =~ NA*crohn + uc  + psc  
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

aid_factor <-usermodel(ldsc_ok, estimation = "DWLS", model = aid_model, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
signif(aid_factor$modelfit,4)



pdf('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure1/three_factor_model.pdf', height = 20, width = 18)
semPaths(semPlotModel_GSEM(aid_factor), 
         what = 'path' , 
         whatLabels= 'est',
         residuals = T, 
         sizeMan = 8, 
         sizeInt = 10,
         label.cex=1, 
         theme="colorblind", 
         rotation = 4, 
         layout = "tree2", 
         sizeLat = 12,
         edge.color = "black",
         edge.label.cex = 1,
         curve = 2,
         nodeLabels = c('CD', 'UC', 'PSC', 'JIA', 
                        'SLE', 'RA','T1D' ,'AST', 'ECZ' , 'F1', 'F2', 'F3' ), 
         label.cex = 1.5, 
         color = list( lat=brewer.pal(12, 'Paired')[c(1,3,5)]),
         groups=list(c("crohn","uc","psc"),c("jia","sle","ra_eu","t1d"), c( "asthma" ,"derma")),
         height = 18, width = 18, 
         edge.label.position=0.5,
         asize=3,
         esize=1,
         
)
dev.off()


pdf('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure1/three_factor_model_fit.pdf', height = 16, width = 16)
grid.table(signif(aid_factor$modelfit,4))
dev.off()

#------ two factor model -----------------------------------------------------

two_f_model <-'F1 =~ NA*crohn + uc  + psc  +jia + sle + ra + t1d 
F2 =~ NA*asthma + derma 

F1~~F2

F1 ~~ 1*F1
F2 ~~ 1*F2

derma~~a*derma
a>0.001
'

two_factor <-usermodel(ldsc_ok, estimation = "DWLS", model = two_f_model, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
two_factor$modelfit


pdf('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure1/two_factor_model.pdf', height = 16, width = 16)
semPaths(semPlotModel_GSEM(two_factor), what = 'path' , 
         whatLabels= 'est',
         residuals = T, 
         sizeMan = 7, 
         sizeInt = 2,
         label.cex=1, 
         theme="colorblind", 
         rotation = 4, 
         layout = "tree2", 
         sizeLat = 7,
         edge.color = "black",
         edge.label.cex = 1,
         curve = 2,
         nodeLabels = c('CD', 'UC', 'PSC', 'JIA', 
                        'SLE', 'RA','T1D' ,'AST', 'ECZ' , 'F1', 'F2'), 
         label.cex = 1.5, 
         #color = list( lat=brewer.pal(12, 'Paired')[c(1,3,5)]),
         #groups=list(c("crohn","uc","psc"),c("jia","sle","ra_eu","t1d"), c( "asthma" ,"derma")),
         height = 10, width = 16, 
         edge.label.position=0.5,
         asize=3,
         esize=1)
dev.off()

pdf('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure1/two_factor_model_fit.pdf', height = 16, width = 16)
grid.table(signif(two_factor$modelfit),4)
dev.off()

#------ one factro model -------------------------------------------------------

one_f_model <-'F1 =~ NA*crohn + uc  + psc  +jia + sle + ra + t1d +asthma + derma 


F1 ~~ 1*F1

'

one_factor <-usermodel(ldsc_ok, estimation = "DWLS", model = one_f_model, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
one_factor$modelfit


pdf('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure1/one_factor_model.pdf', height = 16, width = 16)
semPaths(semPlotModel_GSEM(one_factor), what = 'path' , 
         whatLabels= 'est',
         residuals = T, 
         sizeMan = 7, 
         sizeInt = 2,
         label.cex=1, 
         theme="colorblind", 
         rotation = 4, 
         layout = "tree2", 
         sizeLat = 7,
         edge.color = "black",
         edge.label.cex = 1,
         curve = 2,
         nodeLabels = c('CD', 'UC', 'PSC', 'JIA', 
                        'SLE', 'RA','T1D' ,'AST', 'ECZ' , 'F'), 
         label.cex = 1.5, 
         #color = list( lat=brewer.pal(12, 'Paired')[c(1,3,5)]),
         #groups=list(c("crohn","uc","psc"),c("jia","sle","ra_eu","t1d"), c( "asthma" ,"derma")),
         height = 10, width = 16, 
         edge.label.position=0.5,
         asize=3,
         esize=1)
dev.off()

pdf('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure1/one_factor_model_fit.pdf', height = 16, width = 16)
grid.table(signif(one_factor$modelfit),4)
dev.off()















