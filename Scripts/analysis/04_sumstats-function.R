#  LD score regression on auotoimmunity GWAS
## version3  will be a replication of the version 2 that might be wrong in the 
## allele orientation 

#The tutorial and info on the package and how to run the code are here:  
#https://github.com/GenomicSEM/GenomicSEM/wiki/5.-User-Specified-Models-with-SNP-Effects

#This script, step4, is for preparing the sumstats to run the GWAS by using the two factor model
#specified in step3. 
# 'Please also note that the sumstats function requires that the variables be listed -
#in the same order as they are listed for the ldsc function.'

#This means that I have to run again the analysis from the munge step for only 
#the GWAS we are keeping in the model. 
#The GWAS included in the two factor model are croh, uc, psc, jia, pbc, sle, ra.
#as in V3_step3


#--- Load the libraries--------------------------------------
library(data.table)
library(GenomicSEM)
library(tidyr)
library(dplyr)
library(corrplot)

#----CHeck the standard errors--------------------------------------------------

#----- Function to check if the SE are STD error of log(OD)-----
is_se_logB <- function(BETA,SE, PVALUE) {
    p_calculated <- 2*pnorm((abs(BETA) / SE),lower.tail = F)
    p_reported <- PVALUE
    data.frame(p_calculated, p_reported)
}

#------- Function to compute SE from CI or from Betas-----

add_SE <- function(OR, PVALUE, upper_CI, use_chi = T, add_se = T, sum_stat) {
  
  #calculate SE from Confidence intervals 
  if(use_chi == F){
      #add SE column
      if(add_se == T) {
        sum_stat2 <- mutate(sum_stat, SE = (log(upper_CI) - log(OR) )/qnorm(0.975) )
        return(sum_stat2) }
    
    
    #calculate SE from X2 distribution
  } else {
      chi <-  qchisq(1-PVALUE, df=1)
      se_column <- sqrt((log(OR)^2/chi))

      #add SE and log(OD) column
      if(add_se == T) {
        sum_stat2 <- sum_stat %>% mutate(SE = sqrt(log(OR)^2/chi)) %>% mutate(Beta= log(OR) )
        return(sum_stat2) }
  }
}
   

#---- crohn SE -----------------------------------------------------------------
cd <- fread('outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/crohn_delange-2017.txt', data.table = F)
head(cd) 
is_se_logB(cd$effect, cd$StdErr, cd$p) #SE logistic beta (effect column is beta as there are negative values)

#----- uc SE --------------------------------------------------------------------
uc <- fread('outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/uc_delange-2017.txt', data.table = F)
head(uc)
is_se_logB(uc$effect, uc$StdErr, uc$p) #SE logistic beta (effect column is beta as there are negative values)

#----- psc SE -------------------------------------------------------------------
psc <- fread('outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/psc_ji-2016.txt', data.table = F)
head(psc)

is_se_logB(log(psc$Effect), psc$SE, psc$P) #SE of logistic beta

#----- jia SE ------------------------------------------------------------------

jia  <- fread('outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/jia_lopezisac-2020.txt', data.table = F) 

#there are a bunch of columns we do not care and the SE have to be calculated from the CI
jia_ok <- jia %>% select(-c(effect, all_OR_lower, all_OR_upper, alternate_ids ))
head(jia_ok)

colnames(jia_ok)[8] <- 'beta'
colnames(jia_ok)[9] <- 'se'

fwrite(jia_ok, file='outputs/version3/04_output_sumstats-function/Sumstats_ready_for_munge/jia_lopez-2016_beta_se.txt', 
       sep = '\t', col.names = T, row.names = F, quote = F)

#----- sle SE ------------------------------------------------------------------

sle <- fread('outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/sle_beta_bentham-2015.txt', data.table = F)
head(sle)
is_se_logB(BETA = sle$effect, SE = sle$se, PVALUE = sle$p) #SE of logistic beta (effect column is beta as there are negative values)

#----- ra SE -------------------------------------------------------------------

ra <- fread('outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/ra_okada-2014.txt', data.table = F)
head(ra)
length(which(ra$effect == 1)) #1324526 this means that they have rounded to 1 and when computing z this causes inf and -inf 
dim(ra) #9739303  8
#remove the effects == 1
ra <- ra[!(ra$effect == 1),]
dim(ra) #8414777       8
head(ra)


ra <- add_SE( OR = ra$effect, PVALUE = ra$p, use_chi = T, sum_stat = ra,
                 upper_CI = ra$`OR_95%CIup`, add_se = T)

summary(se_column)
OR_SE_issue <- ra[which(ra$SE ==0 ), ]
ra$CHI <-  qchisq(1- (ra$p) , df=1)

head(ra)
length(which(ra$Beta== 0))
is_se_logB(BETA = ra_ok$Beta, SE = ra_ok$SE, PVALUE = ra_ok$p) #SE of logistic beta
dim(ra)
#remove OD column and CI columns
ra <- ra %>% select(-c(effect))
fwrite(ra, file= 'outputs/version3/04_output_sumstats-function/Sumstats_ready_for_munge/ra_okada-2014_SE_RoundIssues.txt',
       sep = '\t', col.names = T, row.names = F, quote = F)

length(which(ra$SE == 0))

#----ra_ha----------------------------------------------------------------------

ra_ha <- fread('outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/ra_ha_rsID.txt', data.table = F)
dim(ra_ha)
is_se_logB(BETA= ra_ha$effect, SE= ra_ha$standard_error, PVALUE = ra_ha$p) #se of logistic beta


#---- Munge the summary stats --------------------------------------------------

#Function for sample prevalence calculation-------------------------------------
#a function for doing it rapidly from csv files 
#csv files generated by me, from info in the papers
calculate_prevalence <- function(path_of_csv){
  #read the file and remove the NA columns and rows
  preva_csv <- read.csv(file = path_of_csv , header = T, sep = ';')
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
  
  #output the prevalenc when assigned to a variable 
  invisible(round(eff_sample_size, 3))
}

#---- munge function -----------------------------------------------------------

file_names <- c('outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/crohn_delange-2017.txt',
                'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/uc_delange-2017.txt', 
                'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/psc_ji-2016.txt',       
                'outputs/version3/04_output_sumstats-function/Sumstats_ready_for_munge/jia_lopez-2016_beta_se.txt',
                'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/sle_beta_bentham-2015.txt',
                'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/ra_ha_rsID.txt' 
                ) 

crohn_p <- calculate_prevalence('Prevalences/CSV_prevalences/crohn_delange-2017.csv')
uc_p <- calculate_prevalence('Prevalences/CSV_prevalences/uc_delange-2017.csv')
psc_p <- 2871 + 12019
jia_p <- 3305 + 9196 
sle_p <- calculate_prevalence('Prevalences/CSV_prevalences/sle_benthman-2015.csv')
ra_p <- 84687

prevalences <- c(crohn_p, uc_p, psc_p, jia_p, sle_p, ra_p)
trait.names = c('croh', 'uc', 'psc', 'jia', 'sle', 'ra') 
hm3 = 'SNP/w_hm3.snplist'

munge(files = file_names, hm3 = 'SNP/w_hm3.snplist', N = prevalences, trait.names = trait.names )

#----LDSC-----------------------------------------------------------------------

# Load the table that contains all the info on the GWAS
GWAS_info <- readRDS('outputs/version3/GWAS_info_table')
#keep only the one I want ot use here
GWAS_info_2 <- GWAS_info[ c('croh', 'uc', 'psc', 'jia', 'sle', 'ra') , ]

#-----ldsc function ------------------------------------------------------------
munged_files <- c('outputs/version3/04_output_sumstats-function/munged/croh.sumstats.gz', 
                  'outputs/version3/04_output_sumstats-function/munged/uc.sumstats.gz', 
                  'outputs/version3/04_output_sumstats-function/munged/psc.sumstats.gz',
                  'outputs/version3/04_output_sumstats-function/munged/jia.sumstats.gz', 
                  'outputs/version3/04_output_sumstats-function/munged/sle.sumstats.gz', 
                  'outputs/version3/04_output_sumstats-function/munged/ra.sumstats.gz')

names = c('croh', 'uc', 'psc', 'jia', 'sle', 'ra') 
prevalences <- c(crohn_p, uc_p, psc_p, jia_p, sle_p, ra_p)

ldsc_step4 <- ldsc(traits = munged_files, sample.prev = c(.5, .5, 0.1928140, 0.2643788, 0.5, 0.5),  
                   population.prev = GWAS_info_2$population.prev, trait.names = names,
                   ld = "ldscores/eur_w_ld_chr",
                   wld= "ldscores/eur_w_ld_chr", stand = T)

saveRDS(ldsc_step4, 'outputs/version3/04_output_sumstats-function/ldsc_V3_step4')
ldsc4_step4 <- readRDS('outputs/version3/04_output_sumstats-function/ldsc_V3_step4')
rownames(ldsc4_step4$S_Stand) <- colnames(ldsc_step4$S_Stand)
corrplot(ldsc4_step4$S_Stand, order = 'hclust', addCoef.col = 'black', type = 'upper')


#---- Two factor model ---------------------------------------------------------
aid_model <-'F1 =~ NA*croh + uc  + psc 
            F2 =~ NA*jia  + sle + ra
F1~~F2
F1 ~~ 1*F1
F2 ~~ 1*F2'


#run the model
aid_factor <-usermodel(ldsc_step4, estimation = "DWLS", model = aid_model, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)


#save the two factor model results
saveRDS(aid_factor, file = 'outputs/version3/04_output_sumstats-function/aid_twofactor') 
aid_factor <- readRDS('outputs/version3/04_output_sumstats-function/aid_twofactor')

#----Sumstat function-----------------------------------------------------------

file_names <- c('outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/crohn_delange-2017.txt',
                'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/uc_delange-2017.txt', 
                'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/psc_ji-2016.txt',       
                'outputs/version3/04_output_sumstats-function/Sumstats_ready_for_munge/jia_lopez-2016_beta_se.txt',
               # 'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/pbc_cordell-2015_beta_SE.txt',
                'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/sle_beta_bentham-2015.txt',
                'outputs/version3/04_output_sumstats-function/Sumstats_ready_for_munge/ra_okada-2014_SE_RoundIssues.txt') 

prevalences <- c(crohn_p, uc_p, psc_p, jia_p, sle_p, ra_p) #calculated above
se_logit <-c(T,T,T,T,T,T,T)
OLS <- c(F,F,F,F,F,F)
lin_pr <- c(F,F,F,F,F,F)



aid_sumstats_noPBC <-sumstats(files = file_names, ref = 'SNP/reference.1000G.maf.0.005.txt.gz',
                        trait.names = c('croh', 'uc', 'psc', 'jia', 'sle', 'ra'),
                        se.logit =  c(T,T,T,T,T,T),
                        OLS = c(F,F,F,F,F,F), 
                        linprob = c(F,F,F,F,F,F), 
                        N = prevalences, 
                        parallel = T,
                        keep.indel= T 
                        )

saveRDS(aid_sumstats_noPBC, file= 'outputs/version3/04_output_sumstats-function/aid_sumstats_noPBC')
saveRDS(aid_sumstats, file='outputs/version3/04_output_sumstats-function/aid_sumstats')
aid_sumstats <- readRDS('outputs/version3/04_output_sumstats-function/aid_sumstats') 


#----------------------------------------------------------------------------

#----QC plots --------------------------------------------------------------



















