#this scripts is for running a GWAS with GenomicSEM package in a fast way
#it takes the GWAS summary stats and splits into chunks so that it can be run as job array 

counter <-commandArgs(trailingOnly = TRUE)
counter <- as.numeric(counter)

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

aid_sumstats <- readRDS('outputs/version3/04_output_sumstats-function/aid_sumstats')

dim(aid_sumstats)

how_many_chunks(aid_sumstats, 5000)

chunk_to_use <- split_sum_stats(aid_sumstats, 5000, counter ) 

head(chunk_to_use)