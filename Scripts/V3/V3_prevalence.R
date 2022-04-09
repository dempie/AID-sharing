#  LD score regression on auotoimmunity GWAS
## Version3  will be a replication of the version 2 that might be wrong in the 
## calculation of the prevaelnces as reported in the tutorial: 
# https://github.com/GenomicSEM/GenomicSEM/wiki/2.1-Calculating-Sum-of-Effective-Sample-Size-and-Preparing-GWAS-Summary-Statistics

#---- load the csv that i generated taking the info from the papers ------------

#------create a function for prevalence -------------------

calculate_prevalence <- function(path_of_csv){
  #read the file
  preva_csv <- read.csv(file = path_of_csv , header = T, sep = ',')
  preva_csv <-  preva_csv[, c('cases', 'controls')]
  preva_csv_noNA <- preva_csv[complete.cases(preva_csv), ]
  
  #calculated effective prevalence 
    #calculate sample prevalence for each cohort
    preva_csv_noNA$v <-preva_csv_noNA$cases/(preva_csv_noNA$cases+preva_csv_noNA$controls)
    #calculate cohort specific effective sample size
    preva_csv_noNA$EffN<-4*preva_csv_noNA$v*(1-preva_csv_noNA$v)*(preva_csv_noNA$cases+preva_csv_noNA$controls)
    #calculate sum of effective sample size: 
    eff_sample_size <- round(sum(preva_csv_noNA$EffN), 3)
    
  #print number of cases and controls, and effective sample size
  cat(paste0( 'Number of cases  ' , sum(preva_csv_noNA$cases), '\n', 
              'Number of contros  ', sum(preva_csv_noNA$controls), '\n', 
              'Effective sample size  ', round(sum(preva_csv_noNA$EffN), 3)
  ))
  # output the prevalenc when assigned to a variable 
  invisible(round(eff_sample_size, 3))
}





#---- ms imsc ----------------------------------------------------------------
## data generated in the paper 
#calculate sample prevalence for each cohort
# prevalence = cases / (cases + controls)
ms_imsc_sample_size <- ( 9772+ 16849)
ms_imsc_prevalence <-  9772 / ( 9772+ 16849) 

#---- ms andelauer-------------------
#calculate sample prevalence for each cohort





#--------

