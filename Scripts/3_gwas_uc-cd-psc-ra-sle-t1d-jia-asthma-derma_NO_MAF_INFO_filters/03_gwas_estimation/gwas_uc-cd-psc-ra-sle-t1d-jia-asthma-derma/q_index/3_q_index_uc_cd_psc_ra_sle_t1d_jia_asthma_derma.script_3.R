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

chunk_to_use <- split_sum_stats(aid_sumstats, 10000, counter )

aid_model_f1 <- 'F1 =~ NA*crohn + uc  + psc  
F2 =~ NA*jia + sle + ra+ t1d 
F3 =~ NA*asthma+ derma 

F1~~F2
F1~~F3
F2~~F3
F1 ~~ 1*F1
F2 ~~ 1*F2
F3 ~~ 1*F3



F2 ~ SNP
F3 ~ SNP

crohn + uc  + psc ~ SNP


crohn~~a*crohn
uc~~b*uc
psc~~c*psc
jia~~d*jia
sle~~e*sle
ra~~f*ra
t1d~~g*t1d
asthma~~h*asthma
derma~~i*derma

a>0.001
b>0.001
c>0.001
d>0.001
e>0.001
f>0.001
g>0.001
h>0.001
i>0.001
'


aid_model_f2 <- 'F1 =~ NA*crohn + uc  + psc  
F2 =~ NA*jia + sle + ra+ t1d 
F3 =~ NA*asthma+ derma 

F1~~F2
F1~~F3
F2~~F3
F1 ~~ 1*F1
F2 ~~ 1*F2
F3 ~~ 1*F3



F1 ~ SNP
F3 ~ SNP

jia + sle + ra+ t1d ~ SNP

crohn~~a*crohn
uc~~b*uc
psc~~c*psc
jia~~d*jia
sle~~e*sle
ra~~f*ra
t1d~~g*t1d
asthma~~h*asthma
derma~~i*derma

a>0.001
b>0.001
c>0.001
d>0.001
e>0.001
f>0.001
g>0.001
h>0.001
i>0.001'


aid_model_f3 <- 'F1 =~ NA*crohn + uc  + psc  
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

asthma + derma~ SNP

crohn~~a*crohn
uc~~b*uc
psc~~c*psc
jia~~d*jia
sle~~e*sle
ra~~f*ra
t1d~~g*t1d
asthma~~h*asthma
derma~~i*derma

a>0.001
b>0.001
c>0.001
d>0.001
e>0.001
f>0.001
g>0.001
h>0.001
i>0.001'


output <- list()

for(i in 1:3){
  tt <-  c('f1', 'f2', 'f3')[i]
  output[[tt]] <- userGWAS(covstruc = ldsc_model, 
                           SNPs = chunk_to_use, 
                           model = c(aid_model_f1, aid_model_f2, aid_model_f3)[i], 
                           sub=c('F1~~F2'),
                           parallel = F,cores = 1 )
  
  saveRDS(output[[tt]], file.path('/project/aid_sharing/AID_sharing/outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/03_estimation/q_index/', paste0(tt,'/',counter, '_qindex_',tt,'_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_3.RDS')))
  
}



















