#this scripts is for running a GWAS with GenomicSEM package in a fast way
#it takes the GWAS summary stats and splits into chunks so that it can be run as job array 

args<-commandArgs(trailingOnly = TRUE)
counter<-args[1]

#------------------------------------total number of chunks---------------------

# how_many_chunks <- function(summary_stats, chunk_size) {
#   require(data.table)
#   #give a number to the elements of the sum_states
#   my_index <- seq_along(1: nrow(summary_stats))
#   n <- nrow(summary_stats)
#   chunks <- rep(1:ceiling(n/chunk_size),each=chunk_size)[1:n] #the last piece [1:n] is necessary because otherwise the last chunk will become long as the chunk size and will replicate it thus exceeding the length of the sumstats
#   a <- split(my_index, chunks)
#   return(length(a))
# }

#----Function for splitting summary stats --------------------------------------

split_sum_stats <- function(summary_stats, chunk_size, which_chunk) {
  require(data.table)
  #give a number to the elements of the sum_states
  my_index <- seq_along(1: nrow(summary_stats))
  n <- nrow(summary_stats) 
  chunks <- rep(1:ceiling(n/chunk_size),each=chunk_size)[1:n] #the last piece [1:n] is necessary because otherwise the last chunk will become long as the chunk size and will replicate it thus exceeding the length of the sumstats
  a <- split(my_index, chunks)
  summary_stats[ (a[[which_chunk]]) , ]
}

#-------------------------------------------------------------------------------

library(GenomicSEM)
library(data.table)

aid_sumstats <- readRDS('/project/aid_sharing/AID_sharing/outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/02_sumstats_function/sumstats_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_3.RDS')
ldsc_model <- readRDS('/project/aid_sharing/AID_sharing/outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/02_sumstats_function/ldsc_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS')

aid_model <- 'F1 =~ NA*crohn + uc  + psc  
F2 =~ NA*jia + sle + ra+ t1d 
F3 =~ NA*asthma+ derma 

F1~~F2
F1~~F3
F2~~F3
F1 ~~ 1*F1
F2 ~~ 1*F2
F3 ~~ 1*F3

F1 ~ SNP
F2 ~ SNP
F3 ~ SNP

derma~~a*derma
a>0.001'


#dim(aid_sumstats) 3344158      24
#how_many_chunks(aid_sumstats, 10000) #335

chunk_to_use <- split_sum_stats(aid_sumstats, 10000, 203 )
chunk_to_use <- split_sum_stats(chunk_to_use, 1000, counter)

output <- userGWAS(covstruc = ldsc_model, 
                   SNPs = chunk_to_use, 
                   model = aid_model, 
                   sub = c("F1~SNP", "F2~SNP", 'F3~SNP'), 
                   parallel = F , 
                   cores = 1)


saveRDS(output, file.path('/project/aid_sharing/AID_sharing/outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/03_estimation/missing_chunks/chunk_203/', paste0('203_',counter, '_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_3.RDS')))















