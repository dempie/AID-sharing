#this scripts is for running a GWAS with GenomicSEM package in a fast way
#it takes the GWAS summary stats and splits into chunks so that it can be run as job array 

counter <-commandArgs(trailingOnly = TRUE)
counter <- as.numeric(counter)


library(data.table)
my_sum_stats <- as.data.table(c(seq_along(1:10000822)))



#give a number to the elements of the sum_states
index_total <- seq_along(1:nrow(my_sum_stats))
chunk_size <- 10000
n <- nrow(my_sum_stats)

chunks <- rep(1:ceiling(n/chunk_size),each=chunk_size)[1:n]
a <- split(my_sum_stats, chunks)
a[counter]


