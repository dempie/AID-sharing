 #in this script only one function will be run, since I requested to the function to output a log file
#when outputting a log file no warnings will be emitted, SO run this script as local job. 

#RESTART R AT THE END OF THE SCRIPT, ALWAYS! OTHERWISE YOU WILL NOT SEE WARNINGS!!!!

library(MungeSumstats) #Summary statistics here were prepared for being munged by the package MungeSumstats (https://neurogenomics.github.io/MungeSumstats/articles/MungeSumstats.html#overview).
library('BSgenome.Hsapiens.NCBI.GRCh38')
library('SNPlocs.Hsapiens.dbSNP144.GRCh38')
library("SNPlocs.Hsapiens.dbSNP144.GRCh37")
library("BSgenome.Hsapiens.1000genomes.hs37d5")



#-----------run the munging ---------------------------------------------------
#take the names of the traits and put the into a vector
#ATTENTION WHEN REDIRECTING THE LOG FILES OF FORMAT_SUMSTATS, NO WARNINGS WILL BE RETURNED TO THE CONSOLE!!!!
#RESTART R AFTER RUNNING format_sumstats TO TELL R TO SHOW THE ERRORS AND WARNINGS AGAIN !!!!!!!


#the function will infer the genome build. I tested and manualy checked it and it is accurate. 
#very important: the output sumstats will have the effect allele coded as A2 (a-two). . 

a1 <- list.files('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/ready_for_munge/') #get the file names
my_paths <- as.list(paste0('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/ready_for_munge/', a1)) #create paths


builds <- c("GRCH37" ,"GRCH37","GRCH38" ,"GRCH37", "GRCH37" ,"GRCH37" ,"GRCH37" ,"GRCH37" ,"GRCH37" ,"GRCH37" ,"GRCH38","GRCH37" )

gwas_names <- list()
for( i in c(1:length(a1))){
  gwas_names[i] <- strsplit(a1, '_')[[i]][[1]]
  
}

names(my_paths) <- gwas_names #give names just for checking

data('sumstatsColHeaders') #load reference file coming from the package  MungeSumstats
output <- list()

for(i in c(1:12)){
  
  output[[i]] <-  format_sumstats(
    path= my_paths[[i]],
    convert_ref_genome = 'GRCH37',
    ref_genome = builds[i],
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
    remove_multi_rs_snp = T,     
    frq_is_maf = TRUE,
    sort_coordinates = TRUE,
    nThread = 2,
    save_path = paste0('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/munged/', gwas_names[[i]], '_munged_build37.txt' ),
    write_vcf = FALSE,
    tabix_index = FALSE,
    return_data = FALSE,
    return_format = "data.table",
    ldsc_format = FALSE,
    log_folder_ind = FALSE,
    log_folder = 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/log/',
    log_mungesumstats_msgs = TRUE,
    imputation_ind = FALSE,
    force_new = FALSE,
    mapping_file = sumstatsColHeaders
  )
  
}


saveRDS(output, 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/removed_SNP/')
