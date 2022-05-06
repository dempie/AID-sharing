#  LD score regression on auotoimmunity GWAS
## graphing of V3 result 

#---- Libraries ----------------------------------------------------------------

library(corrplot)
library(qgraph)
library(dplyr)

#---- Load the dataset ---------------------------------------------------------

output2 <- readRDS('outputs/version3/04_output_sumstats-function/ldsc_output_04_noRA')

#look at the matrix

rownames(output2$S_Stand) <- colnames(output2$S_Stand) 
corrplot(output2$S_Stand, is.corr = F, addCoef.col = 'black', type = 'upper')


#- heritability --
cbind(colnames(output2$S_Stand), (diag(output2$S)) )

#---- remove columns we do not want --------------------------------------------
#remove asthma_1 (demeanis is the smallest one compared to asthma_2 that is han e al)
#remove Alzh_1 that is the smallest (kunkle), leave Alzh_2 that is the bigger one (Wightman)
#remove m1 that is interacting with everything because of the low h2
#remove armfat that is not needed, it's only my internal control fo heritability (around 0.20)
t1<- output2$S_Stand[ !rownames(output2$S_Stand) %in% c('ms_1', 'asthma_2', 'armfat', 'alzheimer_1') , !colnames(output2$S_Stand) %in% c('ms_1', 'asthma_2', 'armfat', 'alzheimer_1') ]




pre <- colnames(t1)

colnames(t1) <-c('Allergies',  "Multiple sclerosis ", 'Alzheimer\'s disease', 'Asthma', 'Celiac disease',
                 'Crohn\'s disease', 'Juvenile idiopathic arthritis', 'Primary biliary cholangitis', 'Primary sclerosing cholangitis',
                 'Rheumatoid arthritis' , 'Systemic lupus erythematosus', 'Autoimmune thyroid disease', 'Ulcerative colitis')

rownames(t1) <- colnames(t1)

after <- colnames(t1)
rbind(pre, after)

blue_red <- colorRampPalette(c('navy', 'white', 'red3')) 
corrplot(t1, order = 'hclust', addCoef.col = 'black', 
         type = 'upper', tl.col = 'black', col = blue_red(200),
         number.cex=0.75, 
         tl.cex= 0.8)


#export as PDF
pdf('Graphs/Correlation_matrix_V3.pdf', height = 14, width = 14)
corrplot(t1, order = 'hclust', addCoef.col = 'black', 
         type = 'upper', tl.col = 'black', col = blue_red(200),
         number.cex=1, 
         tl.cex= 1.2,
         tl.srt = 45)
dev.off()

#without correlation coefficients
pdf('Graphs/Correlation_matrix_V3_nocoeff.pdf', height = 14, width = 14)
corrplot(t1, order = 'hclust', 
         type = 'upper', tl.col = 'black', col = blue_red(200),
         number.cex=1, 
         tl.cex= 1.2,
         tl.srt = 45)
dev.off()





#---- remove autocorrealtion for qraph

rownames(net$S_Stand) <- colnames(net$S_Stand) 
corrplot(net$S_Stand, is.corr = F, addCoef.col = 'black', type = 'upper')


#- heritability --
cbind(colnames(net2$S_Stand), (diag(ne2$S)) )

#---- remove columns we do not want --------------------------------------------
#remove asthma_1 (demeanis is the smallest one compared to asthma_2 that is han e al)
#remove Alzh_1 that is the smallest (kunkle), leave Alzh_2 that is the bigger one (Wightman)
#remove m1 that is interacting with everything because of the low h2
#remove armfat that is not needed, it's only my internal control fo heritability (around 0.20)
t1<- net$S_Stand[ !rownames(net$S_Stand) %in% c('ms_1', 'asthma_2', 'armfat', 'alzheimer_1') , !colnames(net$S_Stand) %in% c('ms_1', 'asthma_2', 'armfat', 'alzheimer_1') ]




pre <- colnames(t1)

colnames(t1) <-c('Allergies',  "Multiple sclerosis ", 'Alzheimer\'s disease', 'Asthma', 'Celiac disease',
                 'Crohn\'s disease', 'Juvenile idiopathic arthritis', 'Primary biliary cholangitis', 'Primary sclerosing cholangitis',
                 'Rheumatoid arthritis' , 'Systemic lupus erythematosus', 'Autoimmune thyroid disease', 'Ulcerative colitis')

rownames(t1) <- colnames(t1)

after <- colnames(t1)
rbind(pre, after)

pdf('/home/pietro.demela/graph.pdf', height = 14, width = 14)
qgraph(output2$S_Stand, layout= 'spring', vsize=1.5, threshold=0.45, 
       node.width=2, diag=F, label.cex=1, color='black')
dev.off()



#-------- Nicola's function for SEM plots --------------------------------------

semPlotModel_GSEM=function(gsem.object=GWISoutput , est.label="STD_All"){ 
  require(semPlot)
  require(lavaan)
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


#the fucntion creates the object that is then accepted by semPath function from
#semPlot package. Here an example of how to plot it
plot1 <- semPlotModel_GSEM(gsem.object = aid_factor)

#for having different colors for the different groups create a list of groups
f1 <- c('crohn', 'uc',  'psc' )
f2 <- c('jia', 't1d' ,'sle', )
lista <- list(f1, f2)


pdf(file = 'Graphs/V3/Two_factor_model_1', height = 14, width = 14)
semPaths(plot1, what = 'par', residuals = T, groups = lista, 
         color = c('darksalmon', 'darkslategray3'), edge.color = 'black',
         sizeMan = 10, sizeLat = 10, label.cex=2, edge.label.cex = 2)
dev.off()


