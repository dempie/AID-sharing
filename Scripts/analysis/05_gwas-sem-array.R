#this scripts is for running a GWAS with GenomicSEM package in a fast way
#it takes the GWAS summary stats and splits into chunks so that it can be run as job array 

args<-commandArgs(trailingOnly = TRUE)
counter<-args[1]

#------------------------------------total number of chunks---------------------

how_many_chunks <- function(summary_stats, chunk_size) {
  require(data.table)
  #give a number to the elements of the sum_states
  my_index <- seq_along(1: nrow(summary_stats))
  n <- nrow(summary_stats) 
  chunks <- rep(1:ceiling(n/chunk_size),each=chunk_size)[1:n] #the last piece [1:n] is necessary because otherwise the last chunk will become long as the chunk size and will replicate it thus exceeding the length of the sumstats
  a <- split(my_index, chunks)
  return(length(a))
}

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

aid_sumstats <- readRDS('outputs/version3/04_output_sumstats-function/aid_sumstats')
ldsc_step4 <- readRDS('outputs/version3/04_output_sumstats-function/ldsc_V3_step4')
model_ok <- aid_model <-'F1 =~ NA*croh + uc  + psc 
                         F2 =~ NA*jia + pbc + sle + ra
F1~~F2
F1 ~~ 1*F1
F2 ~~ 1*F2
F1 ~ SNP
F2 ~ SNP
'

#dim(aid_sumstats)
#how_many_chunks(aid_sumstats, 5000)

#chunk_to_use <- split_sum_stats(aid_sumstats, 5000, counter )

output <- userGWAS(covstruc = ldsc_step4, 
         SNPs = aid_sumstats[9898:10000, ], 
         model = model_ok, 
         sub = c("F1~SNP", "F2~SNP"))

aid_factor_NA <-usermodel(ldsc_step4, estimation = "DWLS", model = model_ok, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)


#output<-list()
#output[['A']]<-1
#output[['B']]<-2
#counter<-20
save(output, file = file.path('/project/aid_sharing/AID_sharing/outputs/version3/', paste0(counter, ".Robj")))
#load(file.path('/project/aid_sharing/AID_sharing/outputs/version3/', paste0(counter, ".Robj")))



head(aid_sumstats)


















