# in this script all the necessary steps for  preparation for the mungin for the colocalizzation will be performed 
library(data.table)
library(dplyr)
library(MungeSumstats)
library(GenomicRanges)
library(tidyr)




#all the sumstats in the gwas have Beta, except for psc and asthma that has OR, 
#to that ones I just add the beta column

#-------------------------------------------------------------------------------

#---- Create a function to prepare the summary stats ---------------------------
## create a function to prepare the GWAS for the munging depends on dplyr and data.table
##the function changes the name of the columns as required by the library(MungeSumstats) package
#and saves the files in the provided path (if the path is provided)

prepare_munge <- function(sum_stats, rsID, the_effect_allele, the_non_effect_allele, pvalue, the_OR=NA, the_Beta=NA, the_SE=NA, the_chr=NA, the_bp=NA, to_remove=NA, path = NA){
  #an error if arguments are not provided
  if (missing(sum_stats) | missing(rsID) | missing(the_effect_allele) | missing(the_non_effect_allele) |missing(pvalue) ) {

    stop( 'At least one argument is missing')

  } else {

    require(dplyr)
    require(data.table)
    sum_stats <- sum_stats  %>% rename(c(SNP = all_of(rsID), EFFECT_ALLELE  = all_of(the_effect_allele), NON_EFFECT_ALLELE = all_of(the_non_effect_allele), p = all_of(pvalue) ))
    sum_stats$SNP <- tolower(sum_stats$SNP)
    sum_stats$p <- as.numeric(sum_stats$p)
    sum_stats$EFFECT_ALLELE <- toupper(as.character( sum_stats$EFFECT_ALLELE))
    sum_stats$NON_EFFECT_ALLELE <- toupper(as.character( sum_stats$NON_EFFECT_ALLELE))

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

#-------f1 prepare GWAS -------------------------------------------------------
f1 <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/factor1_gwas_final_withQindex.txt', data.table = F)
head(f1)

f1 <- prepare_munge(f1,
                    rsID = 'SNP',
                    the_effect_allele = 'A1',
                    the_non_effect_allele = 'A2',
                    pvalue = 'Pval_Estimate',
                    the_Beta = 'est',
                    the_SE = 'SE',
                    to_remove = c('op', 'free', 'label', 'Z_Estimate', 'chisq', 'chisq_df' , 'chisq_pval', 'AIC', 'error', 'warning' , 'Q_chisq', 'Q_df', 'Q_chisq_pval'),
                    path = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/prepare_for_munge/f1_ready_for_munge.txt'
                    )

head(f1)
#-------f2 prepare GWAS -------------------------------------------------------
f2 <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/factor2_gwas_final_withQindex.txt', data.table = F)
head(f2)

f2 <- prepare_munge(f2,
                    rsID = 'SNP',
                    the_effect_allele = 'A1',
                    the_non_effect_allele = 'A2',
                    pvalue = 'Pval_Estimate',
                    the_Beta = 'est',
                    the_SE = 'SE',
                    to_remove = c('op', 'free', 'label', 'Z_Estimate', 'chisq', 'chisq_df' , 'chisq_pval', 'AIC', 'error', 'warning' , 'Q_chisq', 'Q_df', 'Q_chisq_pval'),
                    path = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/prepare_for_munge/f2_ready_for_munge.txt'
)
head(f2)

#-------f3 prepare GWAS -------------------------------------------------------
f3 <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/factor3_gwas_final_withQindex.txt', data.table = F)
head(f3)

f3 <- prepare_munge(f3,
                    rsID = 'SNP',
                    the_effect_allele = 'A1',
                    the_non_effect_allele = 'A2',
                    pvalue = 'Pval_Estimate',
                    the_Beta = 'est',
                    the_SE = 'SE',
                    to_remove = c('op', 'free', 'label', 'Z_Estimate', 'chisq', 'chisq_df' , 'chisq_pval', 'AIC', 'error', 'warning' , 'Q_chisq', 'Q_df', 'Q_chisq_pval'),
                    path = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/prepare_for_munge/f3_ready_for_munge.txt'
)

#----------uc prepare GWAS -----------------------------------------------------
uc <- fread('Summary_Stats/delange-2017_uc_build37_45975_20161107.txt', data.table = F)
uc <- uc %>% separate(col = 'MarkerName', sep = '_|_', remove = F, convert = T, extra = 'warn', into = c('SNP', 'Extra1', 'Extra2') ) %>% select(-c(Extra1, Extra2))  %>% separate(col='SNP', sep= ':', convert=F,  into = c('CHR', 'BP'), remove = F)
head(uc)    
uc <- select(uc,-c('SNP'))

    prepare_munge(uc,
                    rsID = 'MarkerName',
                    the_effect_allele = 'Allele2',
                    the_non_effect_allele = 'Allele1',
                    pvalue = 'P.value',
                    the_Beta = 'Effect',
                    the_SE = 'StdErr',
                    to_remove = c( 'Min_single_cohort_pval', 'Pval_GWAS3', 'Pval_IIBDGC', 'Pval_IBDseq', 'HetPVal', 'HetDf', 'HetChiSq', 'HetISq'  ),
                    path = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/prepare_for_munge/uc_ready_for_munge_build37.txt'
                    )
head(uc)

#-----crohn prepare GWAS -------------------------------------------------------
cd <- fread('Summary_Stats/delange-2017_cd_build37_40266_20161107.txt', data.table = F)
cd <- cd %>% separate(col = 'MarkerName', sep = '_|_', remove = F, convert = T, extra = 'warn', into = c('SNP', 'Extra1', 'Extra2') ) %>% select(-c(Extra1, Extra2))  %>% separate(col='SNP', sep= ':', convert=F,  into = c('CHR', 'BP'), remove = F)
head(cd)    
cd <- select(cd,-c('SNP'))


prepare_munge(cd,
              rsID = 'MarkerName',
              the_effect_allele = 'Allele2',
              the_non_effect_allele = 'Allele1',
              pvalue = 'P.value',
              the_Beta = 'Effect',
              the_SE = 'StdErr',
              to_remove = c( 'Min_single_cohort_pval', 'Pval_GWAS3', 'Pval_IIBDGC', 'Pval_IBDseq', 'HetPVal', 'HetDf', 'HetChiSq', 'HetISq'  ),
              path = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/prepare_for_munge/cd_ready_for_munge_build37.txt'
)


#-----psc prepare GWAS----------------------------------------------------------
psc <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/psc_ji-2016.txt', data.table = F)
head(psc)


psc$beta <- log(psc$OR)
psc <- rename(psc, the_effect_in_OR=OR )

prepare_munge(psc,
              rsID = 'SNP',
              the_effect_allele = 'A1',
              the_non_effect_allele = 'A2',
              pvalue = 'P',
              the_Beta = 'beta',
              the_SE = 'SE',
              path = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/prepare_for_munge/psc_ready_for_munge_build37.txt'
)


#----jia prepare GWAS ----------------------------------------------------------
jia <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/jia_beta_lopezisac-2020.txt', data.table = F)
head(jia)

prepare_munge(jia,
              rsID = 'SNP',
              the_effect_allele = 'A1',
              the_non_effect_allele = 'A2',
              pvalue = 'p',
              the_Beta = 'Beta',
              the_SE = 'SE',
              path = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/prepare_for_munge/jia_ready_for_munge_build37.txt'
)

#-----sle prepare gwas ---------------------------------------------------------
sle <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/sle_beta_bentham-2015.txt', data.table = F)
head(sle)

prepare_munge(sle,
              rsID = 'SNP',
              the_effect_allele = 'A1',
              the_non_effect_allele = 'A2',
              pvalue = 'p',
              the_Beta = 'Beta',
              the_SE = 'SE',
              path = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/prepare_for_munge/sle_ready_for_munge_build37.txt'
)

#-----t1d prepare gwas ---------------------------------------------------------
t1d <- fread('Summary_Stats/chiou-2021_t1d_build38_GCST90014023_buildGRCh38.tsv', data.table=F)
head(t1d)
prepare_munge(t1d,
              rsID = 'variant_id',
              the_effect_allele = 'effect_allele',
              the_non_effect_allele = 'other_allele',
              pvalue = 'p_value',
              the_Beta = 'beta',
              the_SE = 'standard_error',
              the_chr = 'chromosome',
              the_bp = 'base_pair_location',
              path = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/prepare_for_munge/t1d_ready_for_mung_build38.txt'
)

#-----asthma prepare GWAS-------------------------------------------------------
asthma <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/asthma_ban-2020.txt', data.table = F)
asthma$beta <- log(asthma$OR)
asthma <- rename(asthma, the_effect_in_OR=OR)

prepare_munge(asthma,
              rsID = 'SNP',
              the_effect_allele = 'A1',
              the_non_effect_allele = 'A2',
              pvalue = 'p',
              the_Beta = 'beta',
              the_SE = 'SE',
              path = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/prepare_for_munge/asthma_ready_for_munge_build37.txt'
)


#----------ra preare GWAS ------------------------------------------------------
ra <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/ra_eu_okada-2014_chr_bp.txt', data.table = F)

prepare_munge(ra,
              rsID = 'SNP',
              the_effect_allele = 'A1',
              the_non_effect_allele = 'A2',
              pvalue = 'p',
              the_Beta = 'Beta',
              the_SE = 'SE',
              path = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/prepare_for_munge/ra_ready_for_munge_build37.txt'
)


#-------------derma prepare GWAS
derma <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/derma_sliz-2021.txt', data.table = F)
head(derma)
prepare_munge(derma,
              rsID = 'SNP',
              the_effect_allele = 'A1',
              the_non_effect_allele = 'A2',
              pvalue = 'p',
              the_Beta = 'Beta',
              the_SE = 'SE',
              path = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/prepare_for_munge/derma_ready_for_munge_build38.txt'
)

#---------- MungeSumstats ------------------------------------------------------
#check the build
a1 <- list.files('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/prepare_for_munge/')
my_paths <- as.list(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/prepare_for_munge/', a1))
my_paths <- my_paths [-10] #remove README
a1 <- a1[-10] #remove README
names(my_paths) <- a1 #give names

my_paths <- my_paths[-c(2,12)] #uc and cd do not have rsID and genome build cannot be inferred by the function, I checked in the paper and it is build37

builds <- get_genome_builds(my_paths) #check the genome build, only t1d and derma are build 38

#-------------------------------------------------------------------------------

create_range <- function(res, chr='CHR', start='start', end='end', w=T, tag=NA ) {
  obj <-GRanges( seqnames =  res$chr ,IRanges(names = res$chr ,start = as.numeric(res$start), end = as.numeric(res$end))) 
  if( sum(obj@ranges@width> 2000000) >0 ) { warning(paste0( sum(obj@ranges@width> 2000000)) ,' ranges are wider than 2 MB') 
    print(obj[which(obj@ranges@width> 2000000),])
  }
  ifelse(w, return(obj), return(data.frame(names= obj@ranges@NAMES ,start=obj@ranges@start, width=obj@ranges@width))) 
  
}
#-------------------------------------------------------------------------------


sstats_names <- list.files('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/')
sstats_names  <- sstats_names [c(-8, -9)]
mypaths <- as.list(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/', sstats_names))

trait_names <- list()
for( i in c(1:length(sstats_names ))){
  trait_names[i] <- strsplit(sstats_names , '_')[[i]][[1]]
  
}

#------------------------------------

#this function identifies significant loci with the locus.breaker fuction
#then creates a table of all the loci for all traits and puts them togehter.
#then looks at the overlapps between the loci and creates a macro loci that do not overlapp among them and assigns them a unique number (column pan_loci)
#the function requires the sumstatsto have the columns: SNP, CHR,BP,EFFECT_ALLELE, NON_EFFECT_ALLELE, BETA,SE, P

locus_lister <- function(my_paths, gwas_names) {
    require(data.table)
    require(dplyr)
    require(GenomicRanges)
      loci <- list()
      for(i in c(1:length(my_paths))){
                sstats <- fread(my_paths[[i]], data.table = F)
                colnames(sstats) <- toupper(colnames(sstats))
                
                sstats <- select(sstats, c(SNP, CHR,BP, A2, A1, BETA,SE, P))
                loci[[i]] <-  locus.breaker(sstats)  #locus breaker function is in my R profile
                names(loci[i]) <- gwas_names[[i]]
                loci[[i]]$trait <- rep(gwas_names[i], nrow(loci[[i]]))
      }
      
      all_loci <- do.call(rbind, loci) #create the listof all the loci
      
      pan_loci <- reduce(create_range(all_loci)) #create the genomic ranges object with my function create_range
      pan_loci_non_reduced <- create_range(all_loci)
      
      overlapping <- findOverlaps(pan_loci_non_reduced, pan_loci) #find overlaps between all the loci and use it as an index of the unique non overlapping loci
      overlapping
      all_loci$pan_locus <- rep(0, nrow(all_loci)) #allocate the column 
      all_loci[overlapping@from,]$pan_locus <- overlapping@to  #assing the number as index of which macro loci is overlapping 
      
      return(all_loci)
  
}


#---------generate the list of all loci---------------------------------------------------

loci <- locus_lister(mypaths, trait_names)
dim(loci)

loci_all <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/06_gwas_analysis/loci_all_gwas.RDS')

loci_all$psc$start==loci[loci$trait=='psc','start']

head(loci_all$t1d, 20)
head(loci[loci$trait=='t1d',], 20)










lengt locus.breaker(t1d, p.label ='p', chr.label = 'CHR', pos.label = 'BP')

head(t1d)


