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

prepare_munge<- function(sum_stats, rsID, the_effect_allele, the_non_effect_allele, pvalue, effect_size, to_remove=NA, path = NA){
  #an error if arguments are not provided 
  if (missing(sum_stats) | missing(rsID) | missing(the_effect_allele) | missing(the_non_effect_allele) |missing(pvalue) | missing(effect_size) ) {
    stop( 'At least one argument is missing')
    
  } else {
    
    require(dplyr)
    require(data.table)
    sum_stats <- sum_stats  %>% rename(c(SNP = rsID, A1 = all_of(the_effect_allele), A2 = all_of(the_non_effect_allele), p = all_of(pvalue), effect = all_of(effect_size)))
    #conditional remove
    if(is.na(to_remove[1])){
      #do nothing 
    } else {
      sum_stats <- select(sum_stats,-(all_of(to_remove)))
    }
    
    #save the file if a path is provided
    if(is.na(path)){
      return(sum_stats)
    } else {
      fwrite(sum_stats, path, sep = '\t', col.names = T, row.names = F, quote = F)
      return(sum_stats)
    }
  }
}


#------end of function ---------------------------------------------------------

#----- ms GWAS -----------------------------------------------------------------

head(ms_1)
dim(ms_1) # 472086     11
ms <- prepare_munge(ms_1,  
                    rsID = 'rsid' , 
                    the_effect_allele = 'effect_allele' , ##manually checked from the paper seemed to correspond (when it did not correspond the effect was in the opposite direction)   
                    the_non_effect_allele = 'other_allele', 
                    pvalue = 'p',
                    effect_size = 'OR' , 
                    to_remove = c('beta', 'se' ),
                    path = 'outputs/version3/01_output_prepare-sumstats/01_output_prepare-sumstats/Sumstats_ready_for_munge/ms_imsgc-2011.txt'
                    )
                    
head(ms_1)
dim(ms_1) #472086      8

#---- ms_2 GWAS ----------------------------------------------------------------

head(ms_2)
dim(ms_2) # 7968107       9
ms_x<- prepare_munge(ms_2,  
                    rsID = 'variant_id' , 
                    the_effect_allele = 'other_allele' , #I believe they misreported the alleles, they seem flipped comparing the effect direction with GWAS catalog and the paper 
                    the_non_effect_allele = 'effect_allele', 
                    pvalue = 'p_value',
                    effect_size = 'beta' , 
                    path = 'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/ms_andlauer-2016.txt'
                    )
head(ms_x)
#---- sle GWAS -----------------------------------------------------------------

head(sle)
dim(sle) #7915251      11
sle <- prepare_munge(sle,
                     rsID = 'rsid', 
                     the_effect_allele = 'effect_allele', #manually confirmed on the paper 
                     the_non_effect_allele = 'other_allele', 
                     pvalue = 'p',
                     effect_size = 'OR', 
                     to_remove = c('beta', 'se'),
                     path = 'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/sle_bentham-2015.txt') 

head(sle)
dim(sle) #7915251       9

# the same sumstats but with betas and SE
prepare_munge(sle,
              rsID = 'rsid', 
              the_effect_allele = 'effect_allele', #manually confirmed on the paper 
              the_non_effect_allele = 'other_allele', 
              pvalue = 'p',
              effect_size = 'beta', 
              to_remove = c('OR', 'OR_lower', 'OR_upper'),
              path = 'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/sle_beta_bentham-2015.txt') 


#---- pbc GWAS -----------------------------------------------------------------
head(pbc)
dim(pbc) #1134141      11

pbc_1 <- prepare_munge(pbc,
                     rsID = 'rsid',
                     the_effect_allele = 'effect_allele', #manually checked from the paper seemed to correspond (when it did not correspond the effect was in the opposite direction) 
                     the_non_effect_allele = 'other_allele', 
                     pvalue = 'p', 
                     effect_size = 'OR',
                     to_remove = c('se', 'beta'), 
                     path = 'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/pbc_cordell-2015.txt')

pbc <- fread('Summary_Stats/cordell-2015_pbc_build37_26394269_pbc_efo1001486_1_gwas.sumstats.tsv.gz', data.table = F)
pbc_2 <- prepare_munge(pbc,
                     rsID = 'rsid',
                     the_effect_allele = 'effect_allele', #manually checked from the paper seemed to correspond (when it did not correspond the effect was in the opposite direction) 
                     the_non_effect_allele = 'other_allele', 
                     pvalue = 'p', 
                     effect_size = 'beta',
                     to_remove = c('OR', 'OR_lower', 'OR_upper'), 
                     path = 'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/pbc_cordell-2015_beta_SE.txt')

head(pbc)
dim(pbc) #1134141       8


#---- crohn GWAS----------------------------------------------------------------

head(crohn)
dim(crohn) # 9570787      15

#crohn requires to add the rsIDs
#load the file with the reference SNP and prepare it for merging
referenceSNP <- fread('SNP/reference.1000G.maf.0.005.txt.gz')
referenceSNP <- referenceSNP %>% unite(CHR, BP, sep= ':', na.rm = F, remove = T, col = 'chrPosition' ) %>% select(-c(MAF, A1,A2))

#prepare crohn for merging 
crohn <- crohn %>% separate(col = 'MarkerName', sep = '_|_', remove = T, convert = F, extra = 'warn', into = c('SNP', 'Extra1', 'Extra2') ) %>% select(-c(Extra1, Extra2))  
crohn <- merge.data.table(crohn, referenceSNP, 
                          by.x = 'SNP', by.y = 'chrPosition', 
                          all.x = T, all.y = F, sort = T)
dim(crohn) #9570787      16
colnames(crohn)[1] <- 'Variant_rsID'

#prepare function

crohn <- prepare_munge(crohn, 
                       rsID = 'SNP.y',
                       the_effect_allele = 'Allele2',    #manually checked on the paper for Risk allele
                       the_non_effect_allele = 'Allele1', 
                       pvalue = 'P.value',
                       effect_size = 'Effect',
                       to_remove = c('Variant_rsID', 'Min_single_cohort_pval', 'Pval_GWAS3', 'Pval_IIBDGC', 'Pval_IBDseq', 'HetPVal', 'HetDf', 'HetChiSq', 'HetISq'  ),
                       path = 'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/crohn_delange-2017.txt')
head(crohn)
dim(crohn) #9570787       7


#---- uc GWAS-------------------------------------------------------------------

head(uc)
dim(uc) # 9570787      15

#prepare uc for merging 
uc <- uc %>% separate(col = 'MarkerName', sep = '_|_', remove = T, convert = F, extra = 'warn', into = c('Variant_rsID', 'Extra1', 'Extra2') ) %>% select(-c(Extra1, Extra2))  
uc <- merge.data.table(uc, referenceSNP, 
                          by.x = 'Variant_rsID', by.y = 'chrPosition', 
                          all.x = T, all.y = F, sort = F )
dim(uc) # 9570787      16

#prepare function
uc <- prepare_munge(uc, 
                       rsID = 'SNP',
                       the_effect_allele = 'Allele2', #  manual check seemed to indicate Allele 2 as risk, but not 100% sure
                       the_non_effect_allele = 'Allele1', 
                       pvalue = 'P.value',
                       effect_size = 'Effect',
                       to_remove = c('Variant_rsID', 'Min_single_cohort_pval', 'Pval_GWAS3', 'Pval_IIBDGC', 'Pval_IBDseq', 'HetPVal', 'HetDf', 'HetChiSq', 'HetISq'  ),
                       path = 'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/uc_delange-2017.txt')
head(uc)


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
              effect_size = 'beta_meta',
              path = 'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/armfat.txt'
              )

#----calculate sample prevalence of group 1-------------------------------------
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
              'Effective sample size  ', round(sum(preva_csv_noNA$EffN), 3)
  ))
  
  #output the prevalenc when assigned to a variable 
  invisible(round(eff_sample_size, 3))
}

# calculate the prevalence 
ms_1_p <- 9772 + 16849
ms_2_p <- calculate_prevalence('Prevalences/CSV_prevalences/ms_andlauer-2016.csv')
sle_p <- calculate_prevalence('Prevalences/CSV_prevalences/sle_benthman-2015.csv')
pbc_p <- calculate_prevalence('Prevalences/CSV_prevalences/pbc_cordell-2015.csv')
crohn_p <- calculate_prevalence('Prevalences/CSV_prevalences/crohn_delange-2017.csv')
uc_p <- calculate_prevalence('Prevalences/CSV_prevalences/uc_delange-2017.csv')


#---- Group 1 munge function  --------------------------------------------------

Prevalences <- c(ms_1_p, ms_2_p, sle_p, pbc_p, crohn_p, uc_p)
vector_files <- c('outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/ms_imsgc-2011.txt',
                  'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/ms_andlauer-2016.txt',
                  'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/sle_bentham-2015.txt',
                  'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/pbc_cordell-2015.txt',
                  'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/crohn_delange-2017.txt', 
                  'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/uc_delange-2017.txt'
                  )

munge(vector_files, 
      trait.names = c('ms_1', 'ms_2' , 'sle', 'pbc', 'crohn', 'uc', 'armfat'), 
      hm3 = 'SNP/w_hm3.snplist', 
      N = Prevalences)

armfat_p <- 492874
armfat <- 'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/armfat.txt'
munge(armfat, c('armfat'),   hm3 = 'SNP/w_hm3.snplist', N= c(armfat_p))

#---- END of Group1 ------------------------------------------------------------

#---- Group 2 Load the datasets ------------------------------------------------

asthma_demeanis <- fread('Summary_Stats/demeanis-2017_asthma_build37_TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv', data.table = F)
asthma_han <- fread('Summary_Stats/han-2020_asthma_build7_HanY_prePMID_asthma_UKBB.txt.gz', data.table = F)
celiac <- fread('Summary_Stats/dubois-2010_celiac_build37_20190752_cel_efo0001060_1_gwas.sumstats.tsv.gz', data.table = F)
allergies <- fread('Summary_Stats/ferreira-2017_allergies_build37_SHARE-without23andMe.LDSCORE-GC.SE-META.v0', data.table = F)
psc <- fread('Summary_Stats/ji-2016_psc_build37_ipscsg2016.result.combined.full.with_header.txt', data.table = F)

#---- asthma_demeanis GWAS------------------------------------------------------

head(asthma_demeanis)
dim(asthma_demeanis) #2001280      23

asthma_demeanis <- prepare_munge(asthma_demeanis,
              rsID = 'rsid',
              the_effect_allele = 'alternate_allele', #manually checked on the paper
              the_non_effect_allele = 'reference_allele',
              pvalue = 'European_ancestry_pval_fix',
              effect_size = 'European_ancestry_beta_fix',
              to_remove = c("Multiancestry_beta_fix","Multiancestry_se_fix","Multiancestry_pval_fix" ,"Multiancestry_beta_rand" ,
                             "Multiancestry_se_rand","Multiancestry_pval_rand","Multiancestry_HetQtest","Multiancestry_df_HetQtest",
                             "Multiancestry_pval_HetQtest","European_ancestry_beta_rand","European_ancestry_se_rand",
                             "European_ancestry_pval_rand","European_ancestry_HetQtest","European_ancestry_df_HetQtest",
                             "European_ancestry_pval_HetQtest"),
              path = 'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/asthma_demeanis-2018.txt')

head(asthma_demeanis)
dim(asthma_demeanis) #2001280       8

#---- asthma_han GWAS -----------------------------------------------------------
head(asthma_han)
dim(asthma_han) #9572556      12

asthma_han <- prepare_munge(asthma_han, 
                            rsID = 'SNP',
                            the_effect_allele = 'EA', #manually checked from the paper seemed to correspond (when it did not correspond the effect was in the opposite direction
                            the_non_effect_allele = 'NEA',
                            pvalue = 'P',
                            effect_size = 'OR',
                            path = 'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/asthma_ban-2020.txt'
                            )
head(asthma_han)
dim(asthma_han) #9572556      12

#-----celiac GWAS---------------------------------------------------------------

head(celiac)
dim(celiac) #523398     11

celiac <- prepare_munge(celiac,
                         rsID = 'rsid', 
                          the_effect_allele = 'effect_allele', #manually checked from the paper seemed to correspond
                          the_non_effect_allele = 'other_allele',
                          pvalue = 'p',
                          effect_size = 'beta', 
                          to_remove = c('OR','OR_lower', 'OR_upper'),
                          path = 'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/celiac_dubois-2020.txt'
                          )
head(celiac)
dim(celiac) #523398      8

#---- allergies GWAS -----------------------------------------------------------

head(allergies)
dim(allergies) #8307659      13

#renmae the SNP column in order not to have the same name with the SNP. 
colnames(allergies)[1] <- 'position'
colnames(allergies)[10] <- 'rsID'

fwrite(allergies, 'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/allergies_ferreira-2017.txt',
       col.names = T, row.names = F, sep = '\t', quote = F)


#----- psc GWAS ----------------------------------------------------------------


head(psc)
dim(psc) #7891602      13

#psc columns are not read by munge function, so I create a new dataframe 
psc_ok <- data.frame( 'SNP' = psc$SNP , 
                      'A2' = psc$allele_0  , 
                      'A1' = psc$allele_1 , #manually checked on the paper
                      'Pos' = psc$pos, 
                      'Effect' = psc$or , 
                      'SE' = psc$se ,
                      'P' = psc$p )

fwrite(psc_ok, file='outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/psc_ji-2016.txt',
       col.names = T, row.names = F, sep = '\t', quote = F)


#----calculate sample prevalence of group 2-------------------------------------

asthma_deme_1_p <- calculate_prevalence('Prevalences/CSV_prevalences/asthma_demeanis-2018.csv') #only the Europeans
asthma_han_2_p <- 64538 + 239321 
celiac_p <- calculate_prevalence('Prevalences/CSV_prevalences/celiac_dubois-2010.csv')
allergies_p <- calculate_prevalence('Prevalences/CSV_prevalences/allergies_ferreira-2017.csv')
psc_p <- 2871 + 12019

#---- Group2 munge function ----------------------------------------------------

Prevalences_group2 <- c(asthma_deme_1_p, asthma_han_2_p, celiac_p, allergies_p, psc_p )
vector_files <- c('outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/asthma_demeanis-2018.txt',
                  'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/asthma_ban-2020.txt',
                  'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/celiac_dubois-2020.txt',
                  'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/allergies_ferreira-2017.txt', 
                  'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/psc_ji-2016.txt')

munge(vector_files, 
      trait.names = c('asthma_1', 'asthma_2', 'celiac', 'allergies', 'psc'), 
      hm3 = 'SNP/w_hm3.snplist', 
      N = Prevalences_group2)

#---- END of Group2 ------------------------------------------------------------


#---- Group 3 Load the datasets ------------------------------------------------

alzh_kunkle <- fread('Summary_Stats/kunkle-2019_alzheimer_build37_Kunkle_etal_Stage1_results.txt', data.table = F)
jia <- fread('Summary_Stats/lopezisac-2020_jia_build37_GCST90010715_buildGRCh37.tsv', data.table = F)
ra <- fread('Summary_Stats/okada-2014_RA_build37_RA_GWASmeta_TransEthnic_v2.txt.gz.gz', data.table = F)
thyro <- fread('Summary_Stats/saevarsdottir_auto-thyroid_build37_AITD2020', data.table = F)
alzh_wightman <- fread('Summary_Stats/wightman-2021_alzheimer_build37_PGCALZ2sumstatsExcluding23andMe.txt.gz', data.table = F)
ssc_loper <- fread('Summary_Stats/lopez-2019_ssc_build37_Lopez-Isac_prePMID_META_GWAS_SSc.meta.txt')

#---- alzheimer_kunkle GWAS-----------------------------------------------------

head(alzh_kunkle)
dim(alzh_kunkle) #11480632        8

alzh_kunkle <- prepare_munge(alzh_kunkle,
                               rsID = 'MarkerName',
                               the_effect_allele = 'Effect_allele',  #checked in the paper and also in the readme file
                               the_non_effect_allele = 'Non_Effect_allele', 
                               effect_size = 'Beta', 
                               pvalue = 'Pvalue',
                               path = 'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/alzheimer_kunkle-2019.txt')
head(alzh_kunkle)
dim(alzh_kunkle) #11480632        8

#---- jia GWAS -----------------------------------------------------------------

head(jia)
dim(jia) #7461261      13
#no information about the effect or non-effect allele
#so I checked some rsID as they were reported risk/nonrisk in the paper and on GWAS catalog

jia[ grep('rs6679677', jia$variant_id),] # A listed as risk allele in GWAS catalog and in the paper
jia[ grep('rs7731626', jia$variant_id),] # G listed as risk allele in GWAS catalog and in the paper but with opposite effect
#allela A interpreted as the effect allele I think

jia <- prepare_munge(jia, 
                       rsID = 'variant_id',
                       the_effect_allele = 'alleleB', #checked on the paper (when it did not correspond the effect was in the opposite direction)
                       the_non_effect_allele = 'alleleA',
                       pvalue = 'p_value', 
                       effect_size = 'all_OR' , 
                       path= 'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/jia_lopezisac-2020.txt')

head(jia)
dim(jia) #7461261      13

#---- ra GWAS ------------------------------------------------------------------

head(ra)
dim(ra) #9739303       8

ra <- prepare_munge(ra,
                      rsID = 'SNPID', 
                      the_effect_allele = 'A1', #checked on the paper (when it did not correspond the effect was in the opposite direction)
                      the_non_effect_allele = 'A2' ,
                      effect_size = 'OR(A1)',
                      pvalue = 'P-val', 
                      path= 'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/ra_okada-2014.txt',
                      )
head(ra)
dim(ra) #9739303       8

#---- thyro GWAS ---------------------------------------------------------------

head(thyro)
#thyro columns are not read by munge function
thyro_ok <- data.frame('SNP'= thyro$rsID, 
                            'A1' = thyro$A1,  #checked on the paper 
                            'A2'= thyro$A0, 
                            'P' = thyro$P, 
                            'effect' = thyro$`OR-A1`, 
                            'pos' = thyro$Pos)

#the file is big, remove the rows without rsID
thyro_ok <- thyro_ok[(!is.na(thyro_ok$SNP)),]

fwrite(thyro_ok, file='outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/thyro_saevarsdottir-2020.txt',
       col.names = T, row.names = F, sep = '\t', quote = F)

#---- alzheimer wightman GWAS --------------------------------------------------

head(alzh_wightman)
dim(alzh_wightman)

#add rsID

##load the file with the reference SNP and prepare it for merging
#referenceSNP <- fread('SNP/reference.1000G.maf.0.005.txt.gz')
#referenceSNP <- referenceSNP %>% unite(CHR, BP, sep= ':', na.rm = F, remove = T, col = 'chrPosition' ) %>% select(-c(MAF, A1,A2))

alzh_wightman <- unite(alzh_wightman, chr, PosGRCh37, sep = ':', col = 'chrPosition', remove = T, na.rm = F)
alzh_wightman <- merge.data.table(alzh_wightman, referenceSNP,
                                  by.x = 'chrPosition', by.y = 'chrPosition', 
                                  all.x = T, all.y = F, sort = F)

length(grep('rs',alzh_wightman$SNP)) #9103904 SNPs with an rsID now

alzh_wightman <- prepare_munge(alzh_wightman, 
                                 rsID = 'SNP', 
                                 the_effect_allele = 'testedAllele', #still I am not convinced that this is the effect
                                 the_non_effect_allele = 'otherAllele', 
                                 pvalue = 'p',
                                 effect_size = 'z',
                                 path = 'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/alzh_wightman-2021.txt')
head(alzh_wightman)

#----calculate sample prevalence of group 3-------------------------------------

alzh_kunkle_1_p <- calculate_prevalence('Prevalences/CSV_prevalences/ad_kunkle-2019.csv')
jia_p <- 3305 + 9196 
ra_p <- calculate_prevalence('Prevalences/CSV_prevalences/ra_okada-2014.csv')
thyro_p <- calculate_prevalence('Prevalences/CSV_prevalences/atd_saevrasdpttir-2017.csv')
alzh_wightman_2_p <- calculate_prevalence('Prevalences/CSV_prevalences/ad_wightmn-2019.csv')

Prevalences_group3 <- c(alzh_kunkle_1_p, jia_p, ra_p, thyro_p, alzh_wightman_2_p)
#---- Group 3 munge ------------------------------------------------------------

vector_files <- c('outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/alzheimer_kunkle-2019.txt',
                  'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/jia_lopezisac-2020.txt',
                  'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/ra_okada-2014.txt',
                  'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/thyro_saevarsdottir-2020.txt', 
                  'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/alzh_wightman-2021.txt')

munge(vector_files, 
      trait.names = c('alzheimer_1', 'jia', 'ra', 'thyro', 'alzheimer_2'), 
      hm3 = 'SNP/w_hm3.snplist', 
      N = Prevalences_group3)

#---END of group 3--------------------------------------------------------------

#---Group 4 load the dataset----------------------------------------------------

t1d <- fread('Summary_Stats/chiou-2021_t1d_build38_GCST90014023_buildGRCh38.tsv', data.table = F)


#----t1d GWAS-------------------------------------------------------------------

head(t1d)

#the p_value column is not properly read by munge function as it is a character
#So transform it into as.numeric


t1d_ok <-  data.frame('rsID' = t1d$variant_id,
                  'Beta' = t1d$beta, 
                  'A1' = t1d$effect_allele, 
                  'A2' = t1d$other_allele, 
                  'p' = as.numeric(t1d$p_value),
                  'chr' = t1d$chromosome,
                  'position' = t1d$base_pair_location, 
                  'effect_allele_frequency'= t1d$effect_allele_frequency,
                  'SE' = t1d$standard_error, 
                  'sample_size' = t1d$sample_size
)
 
head(t1d_ok)
sum(t1d_ok$pvalue <0)


fwrite(t1d_ok, 'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/t1d_chiou-2021.txt',
       sep = '\t', col.names = T, row.names = F, quote = F)

trait.names <- c('t1d_chiou-2021')
traits <- c('outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/t1d_chiou-2021.txt') 
prevalence_t1d <- calculate_prevalence('Prevalences/CSV_prevalences/t1d_chiou-2021.csv')
munge(files = traits, hm3 = 'SNP/w_hm3.snplist', N = prevalence_t1d, trait.names = trait.names)


#----derma GWAS-----------------------------------------------------------------
derma <- fread('Summary_Stats/sliz-2021_atopic-dermatitis_build38_GCST90027161_buildGRCh38.tsv.gz', data.table = F)
head(derma) 
prepare_munge(derma, rsID = 'variant_id',
              effect_size = 'beta',
              the_effect_allele = 'effect_allele', 
              the_non_effect_allele = 'other_allele', 
              pvalue = 'p_value', 
              path = 'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/derma_sliz-2021.txt')

trait.names <- c('derma_sliz-2021')
traits <- c('outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/derma_sliz-2021.txt') 
prevalence_derma <- calculate_prevalence('Prevalences/CSV_prevalences/derma_sliz-2021.csv')
munge(files = traits, hm3 = 'SNP/w_hm3.snplist', N = prevalence_derma, trait.names = trait.names)

#-----ra_ha-------------------------------------------------------------

ra_ha <- fread('Summary_Stats/ha-2020_ra_build37_GCST90013534_buildGRCh37.tsv', data.table = T)
head(ra_ha)

#add rsID
referenceSNP <- fread('SNP/reference.1000G.maf.0.005.txt.gz')
referenceSNP <- referenceSNP %>% unite(CHR, BP, sep= ':', na.rm = F, remove = T, col = 'chrPosition' ) %>% select(-c(MAF, A1,A2))
ra_ha <- unite(ra_ha, chromosome, base_pair_location, sep=':', na.rm=F, remove=T, col = 'chrPosition' )

ra_ha_rsID <- merge.data.table(ra_ha, referenceSNP, 
                               by.x = 'chrPosition', by.y = 'chrPosition', 
                               all.x = T, all.y = F, sort = F )
#SE column rename

ra_ha_rsID<- rename(ra_ha_rsID, SE= standard_error )



ra_ha_rsID[grep('rs6705628', ra_ha_rsID$SNP), ]
exp(-0.1246) #0.88285, in the paper they reported 0.88 as OR, so it should be a logistic Beta. 

#effect allele checked on the paper ok!

prepare_munge(ra_ha_rsID, 
              the_effect_allele = 'effect_allele',
              the_non_effect_allele = 'other_allele',
              effect_size = 'beta',
              pvalue = 'p_value',
              rsID = 'SNP',
              path= 'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/ra_ha_rsID.txt')


trait.names <- c('ra_ha-2021')
traits <- c('outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/ra_ha_rsID.txt') 
prevalence_ra_ha <- 84687
munge(files = traits, hm3 = 'SNP/w_hm3.snplist', N = prevalence_ra_ha, trait.names = trait.names)

#------- psoriasis -------------------------------------------------------------
psoriasis <- fread('Summary_Stats/sliz-2021_atopic-dermatitis_build38_GCST90027161_buildGRCh38.tsv.gz')
head(psoriasis)

prepare_munge(psoriasis, 
              the_effect_allele = 'effect_allele',
              the_non_effect_allele = 'other_allele',
              effect_size = 'beta',
              pvalue = 'p_value',
              rsID = 'variant_id',
              path = 'outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/psoriasis_stuart-2022.txt'
)

trait.names <- c('psoriasis_sliz-2021')
traits <- c('outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/psoriasis_stuart-2022.txt') 
prevalence_psoriasis <- 43401
munge(files = traits, hm3 = 'SNP/w_hm3.snplist', N = prevalence_psoriasis, trait.names = trait.names)




