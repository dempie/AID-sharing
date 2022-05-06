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
aid_sumstats <- readRDS('outputs/version3/04_output_sumstats-function/aid_sumstats_2F.RDS')
ldsc_model <- readRDS('outputs/version3/04_output_sumstats-function/ldsc_output_04')
model_ok <- aid_model <-'F1 =~ NA*croh + uc  + psc 
            F2 =~ NA*ra + t1d + jia + sle
F1~~F2
F1 ~~ 1*F1
F2 ~~ 1*F2
F1 ~ SNP
F2 ~ SNP
'

#dim(aid_sumstats)
#how_many_chunks(aid_sumstats, 10000) #429

chunk_to_use <- split_sum_stats(aid_sumstats, 5000, 3 )

output <- userGWAS(covstruc = ldsc_model, 
         SNPs = aid_sumstats[10, ], 
         model = model_ok, 
         sub = c("F1~SNP", "F2~SNP") )

saveRDS(output, file.path('outputs/version3/05_GWAS_results/GWAS_05-05-2022/', paste0(counter, '.RDS')))





















