#  LD score regression on auotoimmunity GWAS
## in this script summary statistics will be prepare to be run in ldsc function


#The tutorial and info on the package and how to run the code are here:  
# https://github.com/GenomicSEM/GenomicSEM/wiki/3.-Models-without-Individual-SNP-effects

#--- Load the libraries--------------------------------------
library(data.table)
library(tidyr)
library(dplyr)
      
#---- Create a function to prepare the summary stats ---------------------------
## create a function to prepare the GWAS for the munging depends on dplyr and data.table
##the function changes the name of the columns as required by the genomicSEM package
#and saves the files in the provided path (if the path is provided)

prepare_munge <- function(sum_stats, rsID, the_effect_allele, the_non_effect_allele, pvalue, the_OR=NA, the_Beta=NA, the_SE=NA, the_chr=NA, the_bp=NA, to_remove=NA, path = NA, the_INFO=NA, the_MAF=NA){
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
          
                  #rename the effect column
                  if(!is.na(the_MAF)){
                    sum_stats <- rename(sum_stats, MAF=all_of(the_MAF))
                    sum_stats$MAF <- as.numeric(sum_stats$MAF)
                  }
          
                  #rename the effect column
                  if(!is.na(the_INFO)){
                    sum_stats <- rename(sum_stats, INFO=all_of(the_INFO))
                    sum_stats$INFO <- as.numeric(sum_stats$INFO)
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







#---- sle GWAS -----------------------------------------------------------------
sle <- fread('Summary_Stats/bentham-2015_sle_build37_26502338_sle_efo0002690_1_gwas.sumstats.tsv.gz', data.table=F)


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
              path =  'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/01_qc_summary_stats/sle_beta_bentham-2015.txt') 

qc_summary_stats(sle, T)


#---- crohn GWAS----------------------------------------------------------------
crohn <- fread('Summary_Stats/delange-2017_cd_build37_40266_20161107.txt', data.table = F)

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

crohn <- crohn %>% select(-c(SNP)) %>%rename( 'SNP'=SNP.y)

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
                       path =  'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/01_qc_summary_stats/crohn_delange-2017.txt')

head(crohn)
dim(crohn) #9570787       7

qc_summary_stats(crohn, T)

#---- uc GWAS-------------------------------------------------------------------
uc <- fread('Summary_Stats/delange-2017_uc_build37_45975_20161107.txt', data.table = F)

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
                      the_effect_allele = 'Allele2', #  manual checked for risk allele
                      the_non_effect_allele = 'Allele1', 
                      pvalue = 'P.value',
                      the_Beta  = 'Effect',
                      the_SE = 'StdErr',
                      the_chr = 'CHR', 
                      the_bp = 'BP',
                      to_remove = c('Min_single_cohort_pval', 'Pval_GWAS3', 'Pval_IIBDGC', 'Pval_IBDseq', 'HetPVal', 'HetDf', 'HetChiSq', 'HetISq'  ),
                      path =  'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/01_qc_summary_stats/uc_delange-2017.txt')


head(uc)
qc_summary_stats(uc, T)


#---- asthma_han GWAS ----------------------------------------------------------

asthma_han <- fread('Summary_Stats/han-2020_asthma_build7_HanY_prePMID_asthma_UKBB.txt.gz', data.table = F)

head(asthma_han)
dim(asthma_han) #9 572 556      12

#add SE of logistic beta to Han-2020 for GWAS estimation
asthma_han$SE <- (log(asthma_han$OR_95U) - log(asthma_han$OR) )/qnorm(0.975)

asthma_han <- prepare_munge(asthma_han, 
                            rsID = 'SNP',
                            the_effect_allele = 'EA', #manually checked from the paper 
                            the_non_effect_allele = 'NEA',
                            pvalue = 'P',
                            the_OR =  'OR',
                            the_chr = 'CHR',
                            the_SE = 'SE',
                            the_bp = 'BP',
                            to_remove = c('EAF', 'INFO', 'N'),
                            path =  'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/01_qc_summary_stats/asthma_han-2020.txt'
                            )
head(asthma_han)
dim(asthma_han) #9572556      13




sum(is.na(asthma_han$SE))
sum(asthma_han$SE==0) #0
summary(asthma_han$SE) #the calculated standard errors stay in the range between 0.005788 and 0.050769 , and no there are no NA or 0 values
summary(log(asthma_han$OR)/asthma_han$SE) #also the ranges of the z scores are acceptable 



#function to caluclate p values and ocmpare with the reported ones
is_se_logB <- function(BETA,SE, PVALUE) {
  p_calculated <- 2*pnorm((abs(BETA) / SE),lower.tail = F)
  p_reported <- PVALUE
  #data.frame(p_calculated, p_reported), 
  a <- lm(p_calculated~p_reported) #fit a linear model to see if they are perfectly correlated 
  print(summary(a))
  
  rand <- sample(1:length(BETA), 10000)
  plot(p_reported[rand], p_calculated[rand]) 
  abline(a, col="red") # regression line (y~x)
  abline(0, 1, col='blue', lty= 4 )
}


is_se_logB(log(asthma_han$OR), asthma_han$SE, asthma_han$p) 




#----- psc GWAS ----------------------------------------------------------------
psc <- fread('Summary_Stats/ji-2016_psc_build37_ipscsg2016.result.combined.full.with_header.txt', data.table = F)


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

fwrite(psc_ok, file= 'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/01_qc_summary_stats/psc_ji-2016.txt',
       col.names = T, row.names = F, sep = '\t', quote = F)

qc_summary_stats(psc_ok, T)



#---jia gwas--------------------------------------------------------------------

jia <- fread('Summary_Stats/lopezisac-2020_jia_build37_GCST90010715_buildGRCh37.tsv', data.table = F)



head(jia)
dim(jia) #7461261      13


jia <- prepare_munge(jia, 
                    rsID = 'variant_id',
                    the_effect_allele = 'alleleB', #checked on the paper
                    the_non_effect_allele = 'alleleA',
                    pvalue = 'p_value', 
                    the_Beta = 'frequentist_add_beta_1' ,
                    the_SE = 'frequentist_add_se_1', 
                    the_chr = 'chromosome', 
                    the_bp = 'position', 
                    to_remove = c( 'all_OR','all_maf', 'all_OR_lower', 'all_OR_upper', 'alternate_ids'),
                    path=  'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/01_qc_summary_stats/jia_beta_lopezisac-2020.txt')


dim(jia) #7461261      13

qc_summary_stats(jia, T)



#----t1d GWAS-------------------------------------------------------------------
t1d <- fread('Summary_Stats/chiou-2021_t1d_build38_GCST90014023_buildGRCh38.tsv', data.table = F, nThread = 4 )
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
                      SE = t1d$standard_error, 
                      effect = as.double(t1d$beta)
)

t1_ok_no_X <- t1d_ok[!(t1d_ok$CHR=='X'),] #remove X chromosome as it will cause the munging to fail!!!!!!!!!!



fwrite(t1_ok_no_X, 'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/01_qc_summary_stats/t1d_chiou-2021.txt',
       sep = '\t', col.names = T, row.names = F)



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
              to_remove =c('effect_allele_frequency'),
              path =  'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/01_qc_summary_stats/derma_sliz-2021.txt')



#-----okada european -----------------------------------------------------------


okada_normal <- fread('Summary_Stats/okada-2014_RA_build37_RA_GWASmeta_TransEthnic_v2.txt.gz.gz', data.table = F)
okada_euro <- fread('Summary_Stats/okada-2014_ra-european_build37_MegaGWAS_summary_European.txt.gz', data.table = F)

head(okada_euro)
dim(okada_euro) #8514610      13

okada_euro<- rename(okada_euro, 
                    SNPID = 'V1',
                    BP ='V2',
                    Neighboring_gene ='V3',
                    A1 = 'V4',
                    A2 = 'V5',
                    No.studies ='V6',
                    No.RAcases ='V7',
                    No.controls ='V8',
                    A1_freq_cases ='V9',
                    A2_freq_controls ='V10',
                    Beta_allele_1 ='V11',
                    SE = 'V12',
                    P_of_allele_1 ='V13'
)

#add CHR column, the format of the Base pair column makes things difficult. Just merge the Bp and CHR position from okada published
okada_normal <- select(okada_normal, c('Chr', 'Position(hg19)', 'SNPID')) 
okada_chr <- merge.data.table(okada_euro, okada_normal, 
                              by.x = 'SNPID', by.y = 'SNPID', 
                              all.x = F , all.y = F, sort = F) #put false so we are sure that each has a positon and chr

head(okada_chr)
dim(okada_chr) #8514610 
okada_chr <- select(okada_chr, -c(BP))
prepare_munge(okada_chr, 
              rsID = 'SNPID',
              the_effect_allele = 'A1',
              the_non_effect_allele = 'A2', 
              the_Beta = 'Beta_allele_1', 
              the_chr = 'Chr', 
              the_bp = 'Position(hg19)',
              the_SE = 'SE', 
              pvalue = 'P_of_allele_1',
              path =  'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/01_qc_summary_stats/ra_eu_okada-2014_chr_bp.txt')

















                        