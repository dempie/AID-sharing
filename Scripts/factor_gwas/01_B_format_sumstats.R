args<-commandArgs(trailingOnly = TRUE)
counter<-args[1]

#in this script only one function will be run, since I requested to the function to output a log file
#when outputting a log file no warnings will be emitted, SO run this script as local job. 

#RESTART R AT THE END OF THE SCRIPT, ALWAYS! OTHERWISE YOU WILL NOT SEE WARNINGS!!!!

library(MungeSumstats) #Summary statistics here were prepared for being munged by the package MungeSumstats (https://neurogenomics.github.io/MungeSumstats/articles/MungeSumstats.html#overview).
library('BSgenome.Hsapiens.NCBI.GRCh38')
library('SNPlocs.Hsapiens.dbSNP144.GRCh38')
library("SNPlocs.Hsapiens.dbSNP144.GRCh37")
library("BSgenome.Hsapiens.1000genomes.hs37d5")



#-----------run the mune function to add rsID and convert everything to the same build -------------
#take the names of the traits and put the into a vector
#ATTENTION WHEN REDIRECTING THE LOG FILES OF FORMAT_SUMSTATS, NO WARNINGS WILL BE RETURNED TO THE CONSOLE!!!!
#RESTART R AFTER RUNNING format_sumstats TO TELL R TO SHOW THE ERRORS AND WARNINGS AGAIN !!!!!!!

#very important: the output sumstats will have the effect allele coded as A2 (a-two). . 

setwd('/project/aid_sharing/AID_sharing')

my_paths <- c("outputs/rev_1/01_prepare/ready_format/asthma_han-2020.txt",
              "outputs/rev_1/01_prepare/ready_format/cd_build38_delange-2017.txt",
              "outputs/rev_1/01_prepare/ready_format/derma_sliz-2021.txt",
              ###
              "outputs/rev_1/01_prepare/ready_format/jia_beta_lopezisac-2020.txt",
              "outputs/rev_1/01_prepare/ready_format/psc_ji-2016.txt",
              "outputs/rev_1/01_prepare/ready_format/ra_eu_okada-2014.txt",
              ###
              "outputs/rev_1/01_prepare/ready_format/sle_beta_bentham-2015.txt",
              "outputs/rev_1/01_prepare/ready_format/t1d_chiou-2021.txt",
              "outputs/rev_1/01_prepare/ready_format/uc_build38_delange-2017.txt"
              ) 

#genome builds 
builds <- c("GRCH37" ,"GRCH38", "GRCH38",
            ###
            "GRCH37","GRCH37","GRCH37",
            ###
            "GRCH37", "GRCH38", "GRCH38" )

#gwas names
gwas_names <- c("asthma",
                "cd",
                "derma",
                ###
                "jia",
                "psc",
                "ra",
                ###
                "sle",
                "t1d",
                'uc')

data('sumstatsColHeaders') #load reference file coming from the package  MungeSumstats

#run function 
format_sumstats(
        path= my_paths[[counter]],
        convert_ref_genome = 'GRCH37',
        ref_genome = builds[counter],
        convert_small_p = F,
        compute_z = FALSE,
        force_new_z = FALSE,
        compute_n = 0L,
        convert_n_int = F,
        analysis_trait = NULL,
        INFO_filter = 0,
        FRQ_filter = 0,
        pos_se = F,   
        effect_columns_nonzero = F,
        N_std = 5,
        N_dropNA = F,
        rmv_chr = c("X", "Y", "MT"),
        rmv_chrPrefix = TRUE,
        on_ref_genome = TRUE,
        strand_ambig_filter = FALSE,
        allele_flip_check = F,
        allele_flip_drop = F,
        allele_flip_z = F,
        allele_flip_frq = F,
        bi_allelic_filter = F,
        snp_ids_are_rs_ids = TRUE,
        remove_multi_rs_snp = T,     
        frq_is_maf = TRUE,
        sort_coordinates = TRUE,
        nThread = 1,
        save_path = paste0('outputs/rev_1/01_prepare/ready_for_munge/', gwas_names[[counter]], '_munged_build37.txt' ),
        write_vcf = FALSE,
        tabix_index = FALSE,
        return_data = FALSE,
        return_format = "data.table",
        ldsc_format = FALSE,
        log_folder_ind = FALSE,
        log_folder = paste0('outputs/rev_1/01_prepare/log/',gwas_names[[counter]] ),
        log_mungesumstats_msgs = TRUE,
        imputation_ind = FALSE,
        force_new = T,
        mapping_file = sumstatsColHeaders
)
  
system('mv *log /project/aid_sharing/AID_sharing/outputs/rev_1/01_prepare/log/log_sbatch')



