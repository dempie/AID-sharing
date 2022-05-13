#  LD score regression on auotoimmunity GWAS
## in this script summary statistics will be prepare to be run in ldsc function


#The tutorial and info on the package and how to run the code are here:  
# https://github.com/GenomicSEM/GenomicSEM/wiki/3.-Models-without-Individual-SNP-effects

#--- Load the libraries--------------------------------------
library(data.table)
library(GenomicSEM)
library(tidyr)
library(dplyr)
        


# Since the dataset are quite heavy in terms of memory I perform the preparation
# of them in groups of 5 GWAS

#---- Group 1 Load the datasets -----------------------------------

ms_1 <- fread('Summary_Stats/imsgc-2011_ms_build37_imsgc_2011_21833088_ms_efo0003885_1_gwas.sumstats.tsv.gz', data.table = F)
ms_2 <- fread('Summary_Stats/andlauer-2016_ms_GCST003566_buildGRCh37.tsv.gz', data.table = F)
sle <- fread('Summary_Stats/bentham-2015_sle_build37_26502338_sle_efo0002690_1_gwas.sumstats.tsv.gz', data.table=F)
pbc <- fread('Summary_Stats/cordell-2015_pbc_build37_26394269_pbc_efo1001486_1_gwas.sumstats.tsv.gz', data.table = F)
crohn <- fread('Summary_Stats/delange-2017_cd_build37_40266_20161107.txt', data.table = F)
uc <- fread('Summary_Stats/delange-2017_uc_build37_45975_20161107.txt', data.table = F)
armfat <- fread('Summary_Stats/continuous-23123-both_sexes-irnt.tsv.gz', data.table = F)



#---- Create a function to prepare the summary stats ---------------------------
## create a function to prepare the GWAS for the munging depends on dplyr and data.table
##the function changes the name of the columns as required by the genomicSEM package
#and saves the files in the provided path (if the path is provided)

prepare_munge <- function(sum_stats, rsID, the_effect_allele, the_non_effect_allele, pvalue, the_OR=NA, the_Beta=NA, the_SE=NA, the_chr=NA, the_bp=NA, to_remove=NA, path = NA){
  #an error if arguments are not provided 
  if (missing(sum_stats) | missing(rsID) | missing(the_effect_allele) | missing(the_non_effect_allele) |missing(pvalue) ) {
   
     stop( 'At least one argument is missing')
    
      } else {
        
          require(dplyr)
          require(data.table)
          sum_stats <- sum_stats  %>% rename(c(SNP = all_of(rsID), A1 = all_of(the_effect_allele), A2 = all_of(the_non_effect_allele), p = all_of(pvalue) ))
          sum_stats$SNP <- tolower(sum_stats$SNP)
          sum_stats$p <- as.numeric(sum_stats$p)
         
          
          #conditional options 
                #remove columns
                if(!is.na(to_remove[1])){sum_stats <- select(sum_stats,-(all_of(to_remove)))} 
                      
                #rename SE column
                if(!is.na(the_SE)){ 
                      sum_stats <- rename(sum_stats, SE=all_of(the_SE))
                      sum_stats$SE <- as.numeric(sum_stats$SE)
                      }
              
                #rename the effect column
                if(!is.na(the_OR)){
                      sum_stats <- rename(sum_stats, OR=all_of(the_OR))
                      sum_stats$OR <- as.numeric(sum_stats$OR)
                      }
                      
                if (!is.na(the_Beta)){
                      sum_stats <- rename(sum_stats, Beta=all_of(the_Beta))
                      sum_stats$Beta <- as.numeric(sum_stats$Beta)
                } 
                
                if(is.na(the_OR) & is.na(the_Beta) ) {stop('Effect column not specified ')}
                
          
                #rename the CHR column
                if(!is.na(the_chr)){ sum_stats <-  rename(sum_stats, CHR=all_of(the_chr))}
              
                #rename the BP column
                if(!is.na(the_bp)){ sum_stats <- rename(sum_stats, BP=all_of(the_bp))}
          
                #save the file if a path is provided
                if(is.na(path)){
                    invisible(sum_stats)
                  
                  
                } else {
                  
                  fwrite(sum_stats, path, sep = '\t', col.names = T, row.names = F, quote = F)
                  invisible(sum_stats)
                }
              }
    }


#----------create a function for QC of the GWAS of the summary stats------------
#this function expects the sum_stats formatted by prepare_munge function

qc_summary_stats <- function(sum_stats, plots=F){
          #compute number of SNPs
          n_SNPs <- nrow(sum_stats)
          
          #compute the number of unique rsIDs 
          sum_stats_rsIDs<- sum_stats[grep('rs',sum_stats$SNP),]
          n_rsIDs <- length(unique(sum_stats_rsIDs$SNP))
          
          SNP_per_chr <- vector(mode='integer', length = 22)
          rsID_per_chr <-  vector(mode='integer', length = 22)
          
          #compute the number of unique rsIDs and SNPs per chromosome 
          for(i in c(1:22)){
            
                SNP_per_chr[i] <- nrow(sum_stats[sum_stats$CHR== i ,]) 
                rsID_per_chr[i] <- length(unique(sum_stats_rsIDs[sum_stats_rsIDs$CHR == i, ]$SNP))
          }
          
          qc_metrics_SNPs <- cbind(Chromosome = c(1:22) ,n_SNPs = SNP_per_chr )
          qc_metrics_rsIDs <- cbind(Chromosome = c(1:22) ,n_rsIDs = rsID_per_chr )
          
          cat(paste0( 'Total number of SNP  ' , n_SNPs , '\n', 
                      'Total number of SNP with rsID  ', n_rsIDs, '\n')
          )
          
          output <- list(qc_SNPs = qc_metrics_SNPs, qc_rsIDs=qc_metrics_rsIDs )
          
          
          
          if(plots==T){

            barplot(t(output$qc_SNPs) , main = 'Number of unique rsID per chromosome',  ylab = 'Number of SNP', 
                    names.arg = output$qc_SNPs[,1], cex.names = 0.8, 
                    legend.text =  paste0( 'Total number of SNP with rsID  ', n_rsIDs))
          }
          
          return(output)   
          
}

#---- calculate prevalence function --------------------------------------------
# https://github.com/GenomicSEM/GenomicSEM/wiki/2.1-Calculating-Sum-of-Effective-Sample-Size-and-Preparing-GWAS-Summary-Statistics

#a function for doing it rapidly from csv files (generated by me, from info in the papers)
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
              'Effective sample size  ', round(sum(preva_csv_noNA$EffN)), '\n')
  )
  
  #output the prevalenc when assigned to a variable 
  invisible(round(eff_sample_size, 3))
}


#------end of function ---------------------------------------------------------

#----- ms GWAS -----------------------------------------------------------------

head(ms_1)
ms <- prepare_munge(ms_1,  
                    rsID = 'rsid' , 
                    the_effect_allele = 'effect_allele' , ##manually checked from the paper seemed to correspond (when it did not correspond the effect was in the opposite direction)   
                    the_non_effect_allele = 'other_allele', 
                    pvalue = 'p',
                    the_Beta = 'beta' , 
                    the_SE = 'se',
                    the_chr = 'chrom', 
                    the_bp = 'pos',
                    to_remove = c('OR', 'OR_lower', 'OR_upper'),
                    path = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/ms_imsgc-2011.txt'
                    )
                    
head(ms)
qc_summary_stats(ms, T) #very few SNPs, not

#---- ms_2 GWAS ----------------------------------------------------------------

head(ms_2)
dim(ms_2) # 7968107       9
ms_x<- prepare_munge(ms_2,  
                    rsID = 'variant_id' , 
                    the_effect_allele = 'other_allele' , #I believe they misreported the alleles, they seem flipped comparing the effect direction with GWAS catalog and the paper 
                    the_non_effect_allele = 'effect_allele', 
                    pvalue = 'p_value',
                    the_Beta = 'beta' , 
                    the_SE = 'standard_error',
                    the_chr = 'chromosome', 
                    the_bp = 'base_pair_location',
                    path = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/ms_andlauer-2016.txt'
                    )
head(ms_x)
qc_summary_stats(ms_x, T)
#---- sle GWAS -----------------------------------------------------------------

# the same sumstats but with betas and SE
head(sle)
sle <- prepare_munge(sle,
              rsID = 'rsid', 
              the_effect_allele = 'effect_allele', #manually confirmed on the paper 
              the_non_effect_allele = 'other_allele', 
              pvalue = 'p',
              the_Beta= 'beta', 
              the_SE = 'se',
              the_chr = 'chrom', 
              the_bp = 'pos',
              to_remove = c('OR', 'OR_lower', 'OR_upper'),
              path =  'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/sle_beta_bentham-2015.txt') 

qc_summary_stats(sle, T)

#---- pbc GWAS -----------------------------------------------------------------
head(pbc)
dim(pbc) #1134141      11

pbc <- prepare_munge(pbc,
                     rsID = 'rsid',
                     the_effect_allele = 'effect_allele', #manually checked from the paper seemed to correspond (when it did not correspond the effect was in the opposite direction) 
                     the_non_effect_allele = 'other_allele', 
                     pvalue = 'p', 
                     the_Beta =  'beta',
                     the_SE = 'se',
                     the_chr = 'chrom', 
                     the_bp = 'pos',
                     to_remove = c('OR', 'OR_lower', 'OR_upper'), 
                     path =  'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/pbc_cordell-2015.txt')

head(pbc)

qc_summary_stats(pbc, T) #very small study, not ok 


#---- crohn GWAS----------------------------------------------------------------

head(crohn)
dim(crohn) # 9570787      15

#crohn requires to add the rsIDs
#load the file with the reference SNP and prepare it for merging
referenceSNP <- fread('SNP/reference.1000G.maf.0.005.txt.gz')
referenceSNP <- referenceSNP %>% unite(CHR, BP, sep= ':', na.rm = F, remove = T, col = 'chrPosition' ) %>% select(-c(MAF, A1,A2))

#prepare crohn for merging 
crohn <- crohn %>% separate(col = 'MarkerName', sep = '_|_', remove = T, convert = F, extra = 'warn', into = c('SNP', 'Extra1', 'Extra2') ) %>% select(-c(Extra1, Extra2))  %>% separate(col='SNP', sep= ':', convert=F,  into = c('CHR', 'BP'), remove = F)
crohn <- merge.data.table(crohn, referenceSNP, 
                          by.x = 'SNP', by.y = 'chrPosition', 
                          all.x = T, all.y = F, sort = T)

dim(crohn) #9570787      18


#prepare function

crohn <- prepare_munge(crohn, 
                       rsID = 'SNP',
                       the_effect_allele = 'Allele2',    #manually checked on the paper for Risk allele
                       the_non_effect_allele = 'Allele1', 
                       pvalue = 'P.value',
                       the_Beta = 'Effect',
                       the_SE = 'StdErr',
                       the_chr = 'CHR', 
                       the_bp = 'BP',
                       to_remove = c( 'Min_single_cohort_pval', 'Pval_GWAS3', 'Pval_IIBDGC', 'Pval_IBDseq', 'HetPVal', 'HetDf', 'HetChiSq', 'HetISq'  ),
                       path =  'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/crohn_delange-2017.txt')

head(crohn)
dim(crohn) #9570787       7

qc_summary_stats(crohn, T)

#---- uc GWAS-------------------------------------------------------------------

head(uc)


#prepare uc for merging 
uc <- uc %>% separate(col = 'MarkerName', sep = '_|_', remove = F, convert = F, extra = 'warn', into = c('SNP', 'Extra1', 'Extra2') ) %>% select(-c(Extra1, Extra2)) %>% separate(col='SNP', sep= ':', convert=F,  into = c('CHR', 'BP'), remove = F)  
uc <- merge.data.table(uc, referenceSNP, 
                          by.x = 'SNP', by.y = 'chrPosition', 
                          all.x = T, all.y = F, sort = T)

uc <- uc %>% select(-c(SNP, MarkerName)) %>%rename( 'SNP'=SNP.y)
#prepare function
uc <- prepare_munge(uc, 
                      rsID = 'SNP',
                      the_effect_allele = 'Allele2', #  manual check seemed to indicate Allele 2 as risk, but not 100% sure
                      the_non_effect_allele = 'Allele1', 
                      pvalue = 'P.value',
                      the_Beta  = 'Effect',
                      the_SE = 'StdErr',
                      the_chr = 'CHR', 
                      the_bp = 'BP',
                      to_remove = c('Min_single_cohort_pval', 'Pval_GWAS3', 'Pval_IIBDGC', 'Pval_IBDseq', 'HetPVal', 'HetDf', 'HetChiSq', 'HetISq'  ),
                      path =  'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/uc_delange-2017.txt')


head(uc)
qc_summary_stats(uc, T)


#---- armfat GWAS --------------------------------------------------------------
head(armfat)

armfat <- armfat %>% unite(chr, pos, sep = ':', remove = T, col = 'chrPosition') 
armfat <- merge.data.table(armfat, referenceSNP, 
                           by.x = 'chrPosition', by.y = 'chrPosition', 
                           all.x = T, all.y = F, sort = F )

#eliminate the snp without rsID otherwise the file is quite big
armfat <- armfat[ (!is.na(armfat$SNP)),]
#prepare_function

armfat <- prepare_munge(armfat, 
              rsID = 'SNP',
              the_effect_allele = 'alt', #readme file said alt is the effect one
              the_non_effect_allele = 'ref',
              pvalue = 'pval_meta', 
              the_Beta =  'beta_meta',
              path =  'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/armfat.txt'
              )

armfat_p <- 492874
armfat <-  'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/armfat.txt'
munge(armfat, c('armfat'),   hm3 = 'SNP/w_hm3.snplist', N= c(armfat_p))

#move output to its folder
system('mv armfat.sumstats.gz armfat_munge.log /outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output ')

#----calculate sample prevalence of group 1-------------------------------------


# calculate the prevalence
ms_1_p <- 9772 + 16849
ms_2_p <- calculate_prevalence('Prevalences/CSV_prevalences/ms_andlauer-2016.csv') #13297.26
sle_p <- calculate_prevalence('Prevalences/CSV_prevalences/sle_benthman-2015.csv') #13218.73
pbc_p <- calculate_prevalence('Prevalences/CSV_prevalences/pbc_cordell-2015.csv') #8380
crohn_p <- calculate_prevalence('Prevalences/CSV_prevalences/crohn_delange-2017.csv') #3923.98
uc_p <- calculate_prevalence('Prevalences/CSV_prevalences/uc_delange-2017.csv') #36082.15



#---- Group 1 munge function  --------------------------------------------------

Prevalences <- c(ms_1_p, ms_2_p, sle_p, pbc_p, crohn_p, uc_p)
vector_files <- c( 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/ms_imsgc-2011.txt',
                   'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/ms_andlauer-2016.txt',
                   'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/sle_beta_bentham-2015.txt',
                   'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/pbc_cordell-2015.txt',
                   'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/crohn_delange-2017.txt',
                   'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/uc_delange-2017.txt'
                  )

munge(vector_files,
      trait.names = c('ms_imgsc-2011', 'ms_andlauer' , 'sle_bentham-2015', 'pbc_cordell-2015', 'crohn_delange-2017', 'uc_delange-2017'),
      hm3 = 'SNP/w_hm3.snplist',
      N = Prevalences)

#move output into its folder
system('mv ms_imgsc-2011.sumstats.gz ms_andlauer.sumstats.gz  sle_bentham-2015.sumstats.gz pbc_cordell-2015.sumstats.gz crohn_delange-2017.sumstats.gz uc_delange-2017.sumstats.gz ms_imgsc-2011_ms* outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output')

#---- END of Group1 ------------------------------------------------------------

#---- Group 2 Load the datasets ------------------------------------------------

asthma_demeanis <- fread('Summary_Stats/demeanis-2017_asthma_build37_TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv', data.table = F)
asthma_han <- fread('Summary_Stats/han-2020_asthma_build7_HanY_prePMID_asthma_UKBB.txt.gz', data.table = F)
celiac <- fread('Summary_Stats/dubois-2010_celiac_build37_20190752_cel_efo0001060_1_gwas.sumstats.tsv.gz', data.table = F)
allergies <- fread('Summary_Stats/ferreira-2017_allergies_build37_SHARE-without23andMe.LDSCORE-GC.SE-META.v0', data.table = F)
psc <- fread('Summary_Stats/ji-2016_psc_build37_ipscsg2016.result.combined.full.with_header.txt', data.table = F)

#---- asthma_demeanis GWAS------------------------------------------------------

head(asthma_demeanis)

dim(asthma_demeanis) #2 001 280      23

asthma_demeanis <- prepare_munge(asthma_demeanis,
              rsID = 'rsid',
              the_effect_allele = 'alternate_allele', #manually checked on the paper
              the_non_effect_allele = 'reference_allele',
              pvalue = 'European_ancestry_pval_fix',
              the_Beta=  'European_ancestry_beta_fix',
              the_chr = 'chr',
              the_bp = 'position', 
              the_SE = 'European_ancestry_se_fix',
              to_remove = c("Multiancestry_beta_fix","Multiancestry_se_fix","Multiancestry_pval_fix" ,"Multiancestry_beta_rand" ,
                             "Multiancestry_se_rand","Multiancestry_pval_rand","Multiancestry_HetQtest","Multiancestry_df_HetQtest",
                             "Multiancestry_pval_HetQtest","European_ancestry_beta_rand","European_ancestry_se_rand",
                             "European_ancestry_pval_rand","European_ancestry_HetQtest","European_ancestry_df_HetQtest",
                             "European_ancestry_pval_HetQtest"),
              path =  'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/asthma_demeanis-2018.txt')

qc_summary_stats(asthma_demeanis, T)
head(asthma_demeanis)
dim(asthma_demeanis) #2 001 280       8

#---- asthma_han GWAS -----------------------------------------------------------
head(asthma_han)
dim(asthma_han) #9 572 556      12

asthma_han$SE <- (log(asthma_han$OR_95U) - log(asthma_han$OR) )/qnorm(0.975)

 asthma_han <- prepare_munge(asthma_han, 
                            rsID = 'SNP',
                            the_effect_allele = 'EA', #manually checked from the paper seemed to correspond (when it did not correspond the effect was in the opposite direction
                            the_non_effect_allele = 'NEA',
                            pvalue = 'P',
                            the_OR =  'OR',
                            the_chr = 'CHR',
                            the_SE = 'SE',
                            the_bp = 'BP',
                            path =  'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/asthma_ban-2020.txt'
                            )
head(asthma_han)
dim(asthma_han) #9572556      12


#add SE of logistic beta to Han-2020 for GWAS estimation


sum(is.na(asthma_han$SE))
sum(asthma_han$SE==0)
summary(asthma_han$SE) #the calculated standard errors stay in the range between 0.005788 and 0.050769 , and no there are no NA or 0 values

summary(log(asthma_han$OR)/asthma_han$SE) #also the ranges of the z scores are acceptable 



#function to caluclate p values and ocmpare with the reported ones
is_se_logB <- function(BETA,SE, PVALUE) {
  p_calculated <- 2*pnorm((abs(BETA) / SE),lower.tail = F)
  p_reported <- PVALUE
  data.frame(p_calculated, p_reported)
}

#check if the pvalues correspond with the new calcualted SE
a <- runif(1, min=1, max=nrow(asthma_han)) #take a randow selection of the dataset and compare the p-values

is_se_logB(log(asthma_han$OR), asthma_han$SE, asthma_han$p) [a:(a+20), ]

#save the asthma_han dataset with Standar errors

#-----celiac GWAS---------------------------------------------------------------

head(celiac)
dim(celiac) #523398     11

celiac <- prepare_munge(celiac,
                         rsID = 'rsid', 
                          the_effect_allele = 'effect_allele', #manually checked from the paper seemed to correspond
                          the_non_effect_allele = 'other_allele',
                          pvalue = 'p',
                          the_SE = 'se', 
                          the_chr = 'chrom', 
                          the_bp = 'pos',
                          the_Beta =  'beta', 
                          to_remove = c('OR','OR_lower', 'OR_upper'),
                          path =  'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/celiac_dubois-2020.txt'
                          )

qc_summary_stats(celiac, T) #very small study, not
head(celiac)
dim(celiac) #523398      8

#---- allergies GWAS -----------------------------------------------------------

head(allergies)
dim(allergies) #8307659      13

#renmae the SNP column in order not to have the same name with the SNP. 
colnames(allergies)[1] <- 'position'


allergies <- prepare_munge(allergies,
          rsID = 'RS_ID', 
          the_effect_allele = 'EFFECT_ALLELE', #manually checked from the paper seemed to correspond
          the_non_effect_allele = 'OTHER_ALLELE',
          pvalue = 'PVALUE',
          the_SE = 'SE', 
          the_chr = 'CHR', 
          the_bp = 'BP',
          the_Beta  = 'BETA', 
          path = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/allergies_ferreira-2017.txt'
) 

qc_summary_stats(allergies, T)


#----- psc GWAS ----------------------------------------------------------------


head(psc)
dim(psc) #7891602      13

#psc columns are not read by munge function, so I create a new dataframe 
psc_ok <- data.frame( 'SNP' = psc$SNP , 
                      'A2' = psc$allele_0  , 
                      'A1' = psc$allele_1 , #manually checked on the paper
                      'Pos' = psc$pos, 
                      'OR' = psc$or , 
                      'SE' = psc$se ,
                      'P' = psc$p , 
                      'CHR'= psc$`#chr`, 
                        'BP'=psc$pos)

fwrite(psc_ok, file= 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/psc_ji-2016.txt',
       col.names = T, row.names = F, sep = '\t', quote = F)

qc_summary_stats(psc_ok, T)

#----calculate sample prevalence of group 2-------------------------------------
# 
asthma_deme_1_p <- calculate_prevalence('Prevalences/CSV_prevalences/asthma_demeanis-2018.csv') #only the Europeans 54627.63
asthma_han_2_p <- 64538 + 239321
celiac_p <- calculate_prevalence('Prevalences/CSV_prevalences/celiac_dubois-2010.csv') #12656.83
allergies_p <- calculate_prevalence('Prevalences/CSV_prevalences/allergies_ferreira-2017.csv') # 206237.7
psc_p <- 2871 + 12019

#---- Group2 munge function ----------------------------------------------------

Prevalences_group2 <- c(asthma_deme_1_p, asthma_han_2_p, celiac_p, allergies_p, psc_p )
vector_files <- c( 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/asthma_demeanis-2018.txt',
                   'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/asthma_ban-2020.txt',
                   'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/celiac_dubois-2020.txt',
                   'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/allergies_ferreira-2017.txt',
                   'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/psc_ji-2016.txt')

munge(vector_files,
      trait.names = c('asthma_demeanis-2018', 'asthma_han-2020', 'celiac_dubois-2010', 'allergies_ferreira-2017', 'psc_ji-2016'),
      hm3 = 'SNP/w_hm3.snplist',
      N = Prevalences_group2)

system('mv asthma_demeanis-2018.sumstats.gz asthma_han-2020.sumstats.gz celiac_dubois-2010.sumstats.gz allergies_ferreira-2017.sumstats.gz psc_ji-2016.sumstats.gz asthma_deme* outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output')

#---- END of Group2 ------------------------------------------------------------


#---- Group 3 Load the datasets ------------------------------------------------

alzh_kunkle <- fread('Summary_Stats/kunkle-2019_alzheimer_build37_Kunkle_etal_Stage1_results.txt', data.table = F)
jia <- fread('Summary_Stats/lopezisac-2020_jia_build37_GCST90010715_buildGRCh37.tsv', data.table = F)
thyro <- fread('Summary_Stats/saevarsdottir_auto-thyroid_build37_AITD2020', data.table = F)

#---- alzheimer_kunkle GWAS-----------------------------------------------------

head(alzh_kunkle)
dim(alzh_kunkle) #11480632        8

alzh_kunkle <- prepare_munge(alzh_kunkle,
                               rsID = 'MarkerName',
                               the_effect_allele = 'Effect_allele',  #checked in the paper and also in the readme file
                               the_non_effect_allele = 'Non_Effect_allele', 
                               the_Beta = 'Beta', 
                                the_SE = 'SE', 
                                the_chr = 'Chromosome', 
                                the_bp = 'Position',
                               pvalue = 'Pvalue',
                               path =  'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/alzheimer_kunkle-2019.txt')
head(alzh_kunkle)
dim(alzh_kunkle) #11480632        8
qc_summary_stats(alzh_kunkle, T)
#---- jia GWAS -----------------------------------------------------------------

head(jia)
dim(jia) #7461261      13
#no information about the effect or non-effect allele
#so I checked some rsID as they were reported risk/nonrisk in the paper and on GWAS catalog

jia[ grep('rs6679677', jia$variant_id),] # A listed as risk allele in GWAS catalog and in the paper with OR 1.36, so alleleB is the effect
jia[ grep('rs7731626', jia$variant_id),] # G listed as risk allele in GWAS catalog and in the paper but with opposite effect (in the paper 1.22)
#allela B interpreted as the effect allele



jia <- prepare_munge(jia, 
                    rsID = 'variant_id',
                    the_effect_allele = 'alleleB', #checked on the paper (when it did not correspond the effect was in the opposite direction)
                    the_non_effect_allele = 'alleleA',
                    pvalue = 'p_value', 
                    the_Beta = 'frequentist_add_beta_1' ,
                    the_SE = 'frequentist_add_se_1', 
                    the_chr = 'chromosome', 
                    the_bp = 'position', 
                    to_remove = c('all_maf'  , 'all_OR', 'all_OR_lower', 'all_OR_upper', 'alternate_ids'),
                    path=  'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/jia_beta_lopezisac-2020.txt')


dim(jia) #7461261      13

qc_summary_stats(jia, T)

#---- thyro GWAS ---------------------------------------------------------------

head(thyro)
#thyro columns are not read by munge function
thyro_ok <- data.frame('SNP'= thyro$rsID, 
                            'A1' = thyro$A1,  #checked on the paper 
                            'A2'= thyro$A0, 
                            'P' = thyro$P, 
                            'OR' = thyro$`OR-A1`, 
                            'BP' = thyro$Pos, 
                              'CHR'=thyro$Chr )


#the file is big, remove the rows without rsID
thyro_ok <- thyro_ok[(!is.na(thyro_ok$SNP)),]

fwrite(thyro_ok, file= 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/thyro_saevarsdottir-2020.txt',
       col.names = T, row.names = F, sep = '\t', quote = F)



#----calculate sample prevalence of group 3-------------------------------------

alzh_kunkle_1_p <- calculate_prevalence('Prevalences/CSV_prevalences/ad_kunkle-2019.csv') #51911
jia_p <- 3305 + 9196
thyro_p <- calculate_prevalence('Prevalences/CSV_prevalences/atd_saevrasdpttir-2017.csv') #114296.2

Prevalences_group3 <- c(alzh_kunkle_1_p, jia_p, thyro_p)
#---- Group 3 munge ------------------------------------------------------------

vector_files <- c( 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/alzheimer_kunkle-2019.txt',
                   'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/jia_beta_lopezisac-2020.txt',
                   'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/thyro_saevarsdottir-2020.txt')

munge(vector_files,
      trait.names = c('alzheimer_kunkle-2019', 'jia_lopezisac-2020', 'thyro_saevarsdottir-2020'),
      hm3 = 'SNP/w_hm3.snplist',
      N = Prevalences_group3)

system('mv alzheimer_kunkle-2019* jia_lopezisac-2020* thyro_saevarsdottir-2020*  outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output')

#---END of group 3--------------------------------------------------------------

#---Group 4 load the dataset----------------------------------------------------

t1d <- fread('Summary_Stats/chiou-2021_t1d_build38_GCST90014023_buildGRCh38.tsv', data.table = F, nThread = 32 )


#----t1d GWAS-------------------------------------------------------------------

head(t1d)
median(t1d$p_value)
median(t1d$beta)

#the p_value column is not properly read by munge function as it is a character
#So transform it into as.numeric, with data.table is the only way that works!!!!!


t1d_ok <-  data.frame(rsID = t1d$variant_id,  
                      A1 = t1d$effect_allele, 
                      A2 = t1d$other_allele, 
                      p = as.double(t1d$p_value),
                      CHR = t1d$chromosome,
                      BP = t1d$base_pair_location, 
                      effect_allele_frequency= t1d$effect_allele_frequency,
                      SE = t1d$standard_error, 
                      sample_size = t1d$sample_size, 
                      effect = as.double(t1d$beta)
)

t1_ok_no_X <- t1d_ok[!(t1d_ok$CHR=='X'),] #remove X chromosome as it will cause the munging to fail!!!!!!!!!!



fwrite(t1_ok_no_X, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/t1d_chiou-2021.txt',
       sep = '\t', col.names = T, row.names = F)

trait.names <- c('t1d_chiou-2021')
traits <- c( 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/t1d_chiou-2021.txt')
prevalence_t1d <- calculate_prevalence('Prevalences/CSV_prevalences/t1d_chiou-2021.csv')  #46453
munge(files = traits, hm3 = 'SNP/w_hm3.snplist', N = prevalence_t1d, trait.names = trait.names)

#move the file in the output folder
system('mv t1d* outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output')


#----derma GWAS-----------------------------------------------------------------
derma <- fread('Summary_Stats/sliz-2021_atopic-dermatitis_build38_GCST90027161_buildGRCh38.tsv.gz', data.table = F)
head(derma) 

derma <- prepare_munge(derma, rsID = 'variant_id',
              the_Beta = 'beta',
              the_effect_allele = 'effect_allele', 
              the_non_effect_allele = 'other_allele',
              the_SE = 'standard_error', 
              the_chr = 'chromosome', 
              the_bp = 'base_pair_location', 
              pvalue = 'p_value', 
              path =  'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/derma_sliz-2021.txt')

trait.names <- c('derma_sliz-2021')
traits <- c( 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/derma_sliz-2021.txt')
prevalence_derma <- calculate_prevalence('Prevalences/CSV_prevalences/derma_sliz-2021.csv')
munge(files = traits, hm3 = 'SNP/w_hm3.snplist', N = prevalence_derma, trait.names = trait.names)

#move the file in the output folder
system('mv derma* outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output')

qc_summary_stats(derma, T)
#------- psoriasis -------------------------------------------------------------
psoriasis <- fread('Summary_Stats/sliz-2021_atopic-dermatitis_build38_GCST90027161_buildGRCh38.tsv.gz')
head(psoriasis)

psoriasis <- prepare_munge(psoriasis, 
              the_effect_allele = 'effect_allele',
              the_non_effect_allele = 'other_allele',
              the_Beta =  'beta',
              pvalue = 'p_value',
              the_SE = 'standard_error' ,
              the_chr = 'chromosome', 
              the_bp = 'base_pair_location',
              rsID = 'variant_id',
              path =  'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/psoriasis_stuart-2022.txt'
)

trait.names <- c('psoriasis_sliz-2021')
traits <- c( 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/psoriasis_stuart-2022.txt')
prevalence_psoriasis <- 43401
munge(files = traits, hm3 = 'SNP/w_hm3.snplist', N = prevalence_psoriasis, trait.names = trait.names)

##move the file in the output folder
system(' mv pso* outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output')

qc_summary_stats(psoriasis, T)

#----vitiligo GWAS ------------------------------------------------------------------

list_vitiligo <- lapply(c(1:22), function(x)fread(paste0('Summary_Stats/jin-2016_vitiligo_build37/GWAS123chr',x,'cmh.txt.gz')))
str(list_vitiligo)

#sum chunks contain 2 columns that we do not need and makes the join impossible (columns V13, V14 and, s)
#remove them

merge_vitiligo <- function(list_chunks){
  
  for(i in c(1:length(list_chunks))){
        
        #remove s column if exist
        if('s' %in% colnames(list_chunks[[i]]) ){
              list_chunks[[i]] <- select(list_chunks[[i]], -c('s'))
        }
        #remove V13 column if exists
        if('V13' %in% colnames(list_chunks[[i]]) ){
              list_chunks[[i]] <- select(list_chunks[[i]], -c('V13'))
        }
        #rmeove V14 column if exists
        if('V14' %in% colnames(list_chunks[[i]]) ){
              list_chunks[[i]] <- select(list_chunks[[i]], -c('V14'))
        }
    
  }    
    
  return(list_chunks)
}

#prepare for merging
a <- merge_vitiligo(list_vitiligo)   

#chr1 has a P value name that it's different from the others, change it
a[[1]]<- rename(a[[1]], 'P' = 'CMH P') 

#merge them
vitiligo <- do.call(rbind, a) 
vitiligo$SNP <- tolower(vitiligo$SNP) #convert to lower case

dim(vitiligo) # 8721242      12
length(unique(vitiligo$CHR)) #22 chromosomes


vitiligo <- prepare_munge(vitiligo, 
              rsID = 'SNP', 
              the_effect_allele = 'A1', #checked on the paper and confirmed
              the_non_effect_allele = 'A2', 
              pvalue = 'P', 
              the_Beta = 'ORX',
              the_SE = 'SE', 
              the_chr = 'CHR',
              the_bp = 'BP',
              path =  'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/vitiligo_jin-2016.txt' )

vitiligo_prev<- calculate_prevalence('Prevalences/CSV_prevalences/vitiligo_jin-2016.csv')

munge( 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/vitiligo_jin-2016.txt',
      hm3 ='SNP/w_hm3.snplist',
      trait.names = 'vitiligo_jin-2016',
      N =  vitiligo_prev)

#move to output folder
system('mv vitiligo* outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output ')

qc_summary_stats(vitiligo, T)


#-----okada european -----------------------------------------------------------


okada_euro <- fread('Summary_Stats/okada-2014_ra-european_build37_MegaGWAS_summary_European.txt.gz', data.table = F)
head(okada_euro)
okada_euro<- rename(okada_euro, 
      SNPID = V1,
       BP =V2,
       Neighboring_gene =V3,
       A1 = V4,
       A2 = V5,
       No.studies =V6,
       No.RAcases =V7,
       No.controls =V8,
       A1_freq_cases =V9,
       A2_freq_controls =V10,
       Beta_allele_1 =V11,
       SE = V12,
       P_of_allele_1 =V13
       )
head(okada_euro)
dim(okada_euro) #8 514 610      13
prepare_munge(okada_euro, 
              rsID = 'SNPID',
              the_effect_allele = 'A1',
              the_non_effect_allele = 'A2', 
              the_Beta = 'Beta_allele_1', 
              pvalue = 'P_of_allele_1',
              path =  'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/ra_eu_okada-2014.txt')


prevalence_ra_okada_eu <- calculate_prevalence('Prevalences/CSV_prevalences/ra_okada-2014_only-eu.csv')

munge( 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/ra_eu_okada-2014.txt',
      hm3 ='SNP/w_hm3.snplist',
      trait.names = 'ra_okada-2014_only-eu',
      N =  prevalence_ra_okada_eu)

system('mv ra_oka* outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/munge_output')

#----- systemic sclerosis GWAS--------------------------------------------------
























                        