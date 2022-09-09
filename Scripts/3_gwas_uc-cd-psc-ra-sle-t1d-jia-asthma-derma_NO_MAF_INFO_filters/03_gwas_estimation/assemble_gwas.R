##  in this script I will assemble the chuks of the GWAS that has been produced

#The tutorial and info on the package and how to run the code are here:  
#https://github.com/GenomicSEM/GenomicSEM/wiki/5.-User-Specified-Models-with-SNP-Effects

library(data.table)
library(qqman)
library(dplyr)

#----- fucntion for loading chunks and assembly---------------------------------


#this function takes the path and name of the of chunks and returns a summary stats for each of the three factors

assemble_f <- function(n_expected_chunks, first_part_of_path, terminal_part_of_path, n_factors ){
  
  
  n_expected_chunks <- c(n_expected_chunks)
  
  #check if files exist and return the one that do not exist
  chunks_found <- lapply(c(n_expected_chunks), function(x)file.exists(paste0(first_part_of_path, x, terminal_part_of_path)))
  missing_chunks <- c(n_expected_chunks)[unlist(chunks_found)==F]
  
  #issue a warning indicating which chunks are missing and that will not be included in the final sumstats
  if((length(missing_chunks))>0) { 
    
    warning(    paste0('SOME CHUNKS ARE MISSING!', '\n',
                       length(n_expected_chunks) - length(missing_chunks), ' chunks were found!', '\n',
                       length(missing_chunks), ' chunks were NOT found!', '\n',
                       'The missing chunks are the ', paste0(missing_chunks, collapse = ' - ') , '\n', 
                       '\n')
    )
    
  } else {
    
    cat(paste0( length(chunks_found) , ' chunks were found!', '\n', 
                'There are not missing chunks :) ','\n' ,
                '\n' ))
  }
  
  
  #load the chunks that have been found
  chunks_to_load <-  c(n_expected_chunks)[unlist(chunks_found)==T]
  
  chunks <- lapply( chunks_to_load , function(x)readRDS(paste0(first_part_of_path, x, terminal_part_of_path )))
  
  cat(paste0('I have loaded ' , length(chunks), ' chunks!', '\n'))
  
  #for each chunk, separate the factors and merge the chunk for each factor
  #allocate the list
  ordered <- list()
  sum_stats <- list()
  
  #cycle through each of the factors (n_factors)
  for(k in c(1:n_factors)) {
    
    #allocate a list in the k element of ordered according to the number of factors
    ordered[[k]] <- list()
    #separate the Factors in the ordered list
    for(i in (1:length(chunks))){
      ordered[[k]][[i]] <- chunks[[i]][[k]]
    }
    
    #merge the chunks for each facotr, put into a list a name the element of the list
    sum_stats[[k]] <-  do.call(rbind, ordered[[k]]) 
    names(sum_stats)[k] <- paste0('factor', k)
    
    #issue a warning if the number of unique SNP is less then the number of rows
    if( length(unique(sum_stats[[k]]$SNP))!= nrow(sum_stats[[k]]) ){warning('The number of unique SNP is different than the number of rows in F',k ,'!!!' )}
    
    #issue a warning if the cumulative numnber of unique SNPs in the merged dataset is different from the sum of the individual unique SNPs per chunk
    if( nrow(sum_stats[[k]]) != sum(unlist(lapply(ordered[[k]], nrow)))  )  warning('The number of unique SNP is different between the merged and the sum of the individual chunks in F',k,  '!!!')
    
  }
  
  #calculate some useful qc information
  SNPs_unique <- list()
  SNP_error <- list()
  
  qc_info <- for( i in (1:length(sum_stats))){
    
    #number of SNP in total without error
    SNPs_unique[[i]]  <- cat(paste0('The number of unique SNPs in F',i, ' is ', length(unique(sum_stats[[i]]$SNP)), '\n'))
    
    #operator found
    cat('The lhs operators found in F',i, ' are ',  unique(sum_stats[[i]]$lhs), '\n')
    
    #number of SNP not estimated
    SNP_error[[i]] <-  cat(paste0('The number of not estimated SNPs in F', i , ' is ', nrow(sum_stats[[i]][sum_stats[[i]]$error != 0,]), '\n' , 
                                  '\n'))
    
  }
  
  
  invisible(sum_stats)
}

#------ append chunks to existing summary stats---------------------------------

#it takes the output of the aseemble_3f function,
#  
append_chunk <- function(list_complete, n_factors ){
  
  #allocate list
  ordered <- list()
  sum_stats <- list()
  
  #cycle through each of the factors (n_factors)
  for(k in c(1:n_factors)) {
    
    #allocate a list in the k element of ordered according to the number of factors
    ordered[[k]] <- list()
    #separate the Factors in the ordered list
    for(i in (1:length(list_complete))){
      ordered[[k]][[i]] <- list_complete[[i]][[k]]
    }
    
    #merge the chunks for each facotr, put into a list a name the element of the list
    sum_stats[[k]] <-  do.call(rbind, ordered[[k]])
    names(sum_stats)[k] <- paste0('factor', k)
    
    #issue a warning if the number of unique SNP is less then the number of rows
    if( length(unique(sum_stats[[k]]$SNP))!= nrow(sum_stats[[k]]) ){warning('The number of unique SNP is different than the number of rows in F',k ,'!!!' )}
    
    #issue a warning if the cumulative numnber of unique SNPs in the merged dataset is different from the sum of the individual unique SNPs per chunk
    if( nrow(sum_stats[[k]]) != sum(unlist(lapply(ordered[[k]], nrow)))  )  warning('The number of unique SNP is different between the merged and the sum of the individual chunks in F',k,  '!!!')
    
  }
  
  #calculate some useful qc information
  SNPs_unique <- list()
  SNP_error <- list()
  
  for( i in (1:length(sum_stats))){
    
    #number of SNP in total without error
    SNPs_unique[[i]]  <- cat(paste0('The number of unique SNPs in F',i, ' is ', length(unique(sum_stats[[i]]$SNP)), '\n'))
    
    #operator found
    cat('The lhs operators found in F',i, ' are ',  unique(sum_stats[[i]]$lhs), '\n')
    
    #number of SNP not estimated
    SNP_error[[i]] <-  cat(paste0('The number of not estimated SNPs in F',i, ' is ', nrow(sum_stats[[i]][sum_stats[[i]]$error != 0,]), '\n' , 
                                  '\n'))
    
  }
  
  
  invisible(sum_stats)
  
}


#------ Assemble, what chunks are missing ---------------------------------------

#there are 335 chunks in the folder, I decided this number in the job array
chunks_1_335<- assemble_f(1:335, 'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/03_estimation/chunks/', '_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_3.RDS' , 3)

# 330 chunks were found!
#   5 chunks were NOT found!
#   The missing chunks are the 12 - 131 -196 - 203 - 250

#--------missing chunks---------------------------------------------------------

chunk_12 <- assemble_f(1:10,  'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/03_estimation/missing_chunks/chunk_12/12_', '_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_3.RDS', 3 )
chunk_131 <- assemble_f(1:10,  'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/03_estimation/missing_chunks/chunk_131/131_', '_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_3.RDS', 3 )
chunk_196 <- assemble_f(1:10,  'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/03_estimation/missing_chunks/chunk_196/196_', '_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_3.RDS', 3 )
chunk_203 <- assemble_f(1:10,  'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/03_estimation/missing_chunks/chunk_203/203_', '_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_3.RDS', 3 )
chunk_250 <- assemble_f(1:10,  'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/03_estimation/missing_chunks/chunk_250/250_', '_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_3.RDS', 3 )


#put all the chunks together and save the output 
all_chunks <- append_chunk(list(chunks_1_335, chunk_12, chunk_131, chunk_196, chunk_203, chunk_250), 3) #The number of unique SNPs in F123 is 3343158

saveRDS(all_chunks, 'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/03_estimation/factors_summary_stats.RDS')

all_chunks <- readRDS( 'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/03_estimation/factors_summary_stats.RDS')


#-------- save the unique summary stata-----------------------------------------------------

F1 <- all_chunks$factor1
F2 <- all_chunks$factor2
F3 <- all_chunks$factor3

F1_noNA <- F1[(which(!is.na(F1$Pval_Estimate))),]
dim(F1_noNA)

F2_noNA <- F2[(which(!is.na(F2$Pval_Estimate))),]
dim(F2_noNA)

F3_noNA <- F3[(which(!is.na(F3$Pval_Estimate))),]
dim(F2_noNA)


#save the factor gwas individually (clean format)
for( i in c(1:length( list(F1_noNA, F2_noNA, F3_noNA)))){
  to_save <- list(F1_noNA, F2_noNA, F3_noNA)[[i]]
  to_save <- to_save[,c('SNP', 'CHR', 'BP', 'A1', 'A2', 'est', 'SE', 'Pval_Estimate')]
  to_save <- rename(to_save, Beta ='est', P = 'Pval_Estimate')
  fwrite(to_save, paste0('outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/03_estimation/summarystats_f',i,'_.txt'), 
         col.names = T, row.names = F, sep = '\t', quote = F)
  rm(to_save)
}

#------- calculate effective sample size of the individual factors -------------

#restrict to MAF of 40% and 10%
N_hat <- vector()
for( i in c(1:length( list(F1_noNA, F2_noNA, F3_noNA)))){
  to_use <- list(F1_noNA, F2_noNA, F3_noNA)[[i]]
  maf_filtered <-subset(to_use , to_use$MAF <= .4 & to_use$MAF >= .1)
  N_hat[i]<-mean(1/((2*to_use$MAF*(1-to_use$MAF))*to_use$SE^2))
  names(N_hat)[i] <- paste0('N_hat F', i)
}

saveRDS(N_hat,'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/03_estimation/effective_sample_size_factors.RDS')




