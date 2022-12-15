library(data.table)
library(dplyr)
library(MungeSumstats)
library(tidyr)

library(MungeSumstats) #Summary statistics here were prepared for being munged by the package MungeSumstats (https://neurogenomics.github.io/MungeSumstats/articles/MungeSumstats.html#overview).
library('BSgenome.Hsapiens.NCBI.GRCh38')
library('SNPlocs.Hsapiens.dbSNP144.GRCh38')
library("SNPlocs.Hsapiens.dbSNP144.GRCh37")
library("BSgenome.Hsapiens.1000genomes.hs37d5")


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
    sum_stats <- sum_stats  %>% dplyr::rename(c(SNP = all_of(rsID), EFFECT_ALLELE  = all_of(the_effect_allele), NON_EFFECT_ALLELE = all_of(the_non_effect_allele), p = all_of(pvalue) ))
    sum_stats$SNP <- tolower(sum_stats$SNP)
    sum_stats$p <- as.numeric(sum_stats$p)
    sum_stats$EFFECT_ALLELE <- toupper(as.character( sum_stats$EFFECT_ALLELE))
    sum_stats$NON_EFFECT_ALLELE <- toupper(as.character( sum_stats$NON_EFFECT_ALLELE))
    
    #conditional options
    #remove columns
    if(!is.na(to_remove[1])){sum_stats <- select(sum_stats,-(all_of(to_remove)))}
    
    #rename SE column
    if(!is.na(the_SE)){
      sum_stats <- dplyr::rename(sum_stats, SE=all_of(the_SE))
      sum_stats$SE <- as.numeric(sum_stats$SE)
    }
    
    #rename the effect column
    if(!is.na(the_OR)){
      sum_stats <-dplyr::rename(sum_stats, OR=all_of(the_OR))
      sum_stats$OR <- as.numeric(sum_stats$OR)
    }
    
    if (!is.na(the_Beta)){
      sum_stats <- dplyr::rename(sum_stats, Beta=all_of(the_Beta))
      sum_stats$Beta <- as.numeric(sum_stats$Beta)
    }
    
    if(is.na(the_OR) & is.na(the_Beta) ) {stop('Effect column not specified ')}
    
    
    #rename the CHR column
    if(!is.na(the_chr)){ sum_stats <-  dplyr::rename(sum_stats, CHR=all_of(the_chr))}
    
    #rename the BP column
    if(!is.na(the_bp)){ sum_stats <- dplyr::rename(sum_stats, BP=all_of(the_bp))}
    
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
                    to_remove = c('op', 'free', 'label', 'Z_Estimate', 'chisq', 'chisq_df' , 'chisq_pval', 'AIC', 'error', 'warning' , 'Q_chisq', 'Q_df'),
                    path = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/q_index_munged/f1_ready_for_munge_with_Q.txt'
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
                    to_remove = c('op', 'free', 'label', 'Z_Estimate', 'chisq', 'chisq_df' , 'chisq_pval', 'AIC', 'error', 'warning' , 'Q_chisq', 'Q_df'),
                    path = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/q_index_munged/f2_ready_for_munge_with_Q.txt'
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
                    to_remove = c('op', 'free', 'label', 'Z_Estimate', 'chisq', 'chisq_df' , 'chisq_pval', 'AIC', 'error', 'warning' , 'Q_chisq', 'Q_df'),
                    path = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/q_index_munged/f3_ready_for_munge_with_Q.txt'
)



#---------munge function -------------------------------------------------------


a1 <- list.files('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/q_index_munged/') #get the file names
my_paths <- as.list(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/q_index_munged/', a1)) #create paths



gwas_names <- list()
for( i in c(1:length(a1))){
  gwas_names[i] <- strsplit(a1, '_')[[i]][[1]]
  
}

names(my_paths) <- gwas_names #give names just for checking

data('sumstatsColHeaders') #load reference file coming from the package  MungeSumstats
output <- list()

for(i in c(1:3)){
  
  output[[i]] <-  format_sumstats(
    path= my_paths[[i]],
    convert_ref_genome = 'GRCH37',
    ref_genome = 'GRCH37',
    convert_small_p = F,
    compute_z = FALSE,
    force_new_z = FALSE,
    compute_n = 0L,
    convert_n_int = F,
    analysis_trait = NULL,
    INFO_filter = 0,
    FRQ_filter = 0,
    pos_se = T,   #importan for coloc
    effect_columns_nonzero = T, #important for coloc 
    N_std = 5,
    N_dropNA = F,
    rmv_chr = c("X", "Y", "MT"),
    rmv_chrPrefix = TRUE,
    on_ref_genome = TRUE,
    strand_ambig_filter = FALSE,
    allele_flip_check = TRUE,
    allele_flip_drop = TRUE,
    allele_flip_z = TRUE,
    allele_flip_frq = TRUE,
    bi_allelic_filter = TRUE,
    snp_ids_are_rs_ids = TRUE,
    remove_multi_rs_snp = T,    #important for my script 
    frq_is_maf = TRUE,
    sort_coordinates = TRUE,
    nThread = 12,
    save_path = paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/q_index_munged/', gwas_names[[i]], '_munged_q_index_build37.txt' ),
    write_vcf = FALSE,
    tabix_index = FALSE,
    return_data = FALSE,
    return_format = "data.table",
    ldsc_format = FALSE,
    log_folder_ind = FALSE,
    log_folder = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/q_index_munged',
    log_mungesumstats_msgs = TRUE,
    imputation_ind = FALSE,
    force_new = FALSE,
    mapping_file = sumstatsColHeaders
  )
  
}


saveRDS(output, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/q_index_munged/removed_SNP.RDS')




