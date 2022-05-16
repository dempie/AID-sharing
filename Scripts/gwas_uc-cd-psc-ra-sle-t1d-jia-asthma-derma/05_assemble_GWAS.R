#  LD score regression on auotoimmunity GWAS
##  in this script I will assemble the chuks of the GWAS that has been produced
#in /sbatch_scripts/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma

#The tutorial and info on the package and how to run the code are here:  
#https://github.com/GenomicSEM/GenomicSEM/wiki/5.-User-Specified-Models-with-SNP-Effects

library(data.table)
library(qqman)


#----- fucntion for loading chunks and assembly---------------------------------


#this function takes the path and name of the of chunks and returns a summary stats for each of the three factors

assemble_3f <- function(n_expected_chunks, first_part_of_path, terminal_part_of_path){
  
  
      n_expected_chunks <- c(n_expected_chunks)
        
      #check if files exist and return the one that do not exist
      chunks_found <- lapply(c(n_expected_chunks), function(x)file.exists(paste0(first_part_of_path, x, terminal_part_of_path )))
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
      
      #for each chunk, separate F1 from F2 from F3
      
        #allocate three lists
        all_F1 <- list() #allocate the list
        all_F2 <- list() #allocate the list
        all_F3 <- list() #allocate the list
        
      
      #create a list of lists, in which each element of the list is a chunk, three lists one per factor
      for(i in (1:length(chunks))){
        all_F1[[i]] <- chunks[[i]][[1]]
        all_F2[[i]] <- chunks[[i]][[2]]
        all_F3[[i]] <- chunks[[i]][[3]]
      }
      
      #create the dataframe for F1 and F2 and F3
      F1 <-do.call(rbind, all_F1) 
      F2 <- do.call(rbind, all_F2)
      F3 <- do.call(rbind, all_F3)
      
      #create a list of summary stats, each element of the list are the summary stats
      sum_stats <- list(sumstats_F1 = F1, sumstats_F2 = F2, sumstats_F3 = F3)
      
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
      
      
      #create the list of objects to be returned
      return_object <- list(sumstats_F1 = F1, sumstats_F2 = F2, sumstats_F3 = F3, SNPs_unique = SNPs_unique, SNP_error = SNP_error)
      
      invisible(return_object)
}

#------ append chunks to existing summary stats---------------------------------

#it takes the output of the aseemble_3f function,
#  
# REMBER TO remove all the accesory outputs from assemble_3f

append_chunk <- function(list_complete ){
      
      
    #allocate list
    list_F1 <- list()
    list_F2 <- list()
    list_F3 <- list()
    
    
      for( i in c(1:length(list_complete))){
        list_F1[[i]] <- list_complete[[i]][[1]]
        list_F2[[i]] <- list_complete[[i]][[2]]
        list_F3[[i]] <- list_complete[[i]][[3]]
      }
        
        
      #merge the respective dataframes
      F1 <-do.call(rbind, list_F1) 
      F2 <- do.call(rbind, list_F2)
      F3 <- do.call(rbind, list_F3)
      
      sum_stats <- list(sumstats_F1 = F1, sumstats_F2 = F2, sumstats_F3 = F3)
      
      
      #calculate some useful qc information
      SNPs_unique <- list()
      SNP_error <- list()
      
      qc_info <- for( i in (1:length(sum_stats))){
            
            #number of SNP in total without error
            SNPs_unique[[i]]  <- cat(paste0('The number of unique SNPs in F',i, ' is ', length(unique(sum_stats[[i]]$SNP)), '\n'))
            
            #operator found
            cat('The lhs operators found in F',i, ' are ',  unique(sum_stats[[i]]$lhs), '\n')
            
            #number of SNP not estimated
            SNP_error[[i]] <-  cat(paste0('The number of not estimated SNPs in F',i, ' is ', nrow(sum_stats[[i]][sum_stats[[i]]$error != 0,]), '\n' , 
                                          '\n'))
        
      }
      
             #issue a warning if the cumulative numnber of unique SNPs in the merged dataset is different from the sum of the individual unique SNPs per chunk
              row_F1 <- vector()
              row_F2 <- vector()
              row_F3 <- vector()
            for( i in c(1:length(list_complete))){
                row_F1[i] <- nrow(list_complete[[i]][[1]])
                row_F2[i] <- nrow(list_complete[[i]][[2]])
                row_F3[i] <- nrow(list_complete[[i]][[3]])
              }  
              
              
              if(!identical(sum(row_F1),nrow(F1),sum(row_F2),nrow(F2),sum(row_F3),nrow(F3))){
                
                warning('The number of unique SNP is different between the merged and the sum of the individual chunks!!!')
                
              }
      
      #create a list of summary stats, each element of the list are the summary stats
      sum_stats <- list(sumstats_F1 = F1, sumstats_F2 = F2, sumstats_F3 = F3)
      
      return(sum_stats)
      
}

#------ Assemble, what chunks are missing ---------------------------------------

#there are 335 chunks in the folder, I decided this number in the job array
chunks_1_335<- assemble_3f(1:335, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/chunks/', '_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS' )

# 329 chunks were found!
# 6 chunks were NOT found!
# The missing chunks are 12 - 25 - 87 - 131 - 194 - 203


#chunk 12 : task 144 failed - "system is computationally singular: reciprocal condition number = 3.33818e-38"
#chunk25 : Error in unserialize(socklist[[n]]) : error reading from connection. for this just sent the job again
#chunk 87 : Error in unserialize(socklist[[n]]) : error reading from connection. for this just sent the job again
#chunk 131: task 211 failed - "system is computationally singular: reciprocal condition number = 1.50212e-24"
#chunk 194: Error in unserialize(socklist[[n]]) : error reading from connection. for this just sent the job again
#chunk 203:  task 148 failed - "system is computationally singular: reciprocal condition number = 2.86735e-37"


#--------missing chunks---------------------------------------------------------

#chunk 12 was divided in 10 subchuks and re-estimated 

chunk_12 <- assemble_3f(1:10, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/chunks/chunk_12/chunk_12_1-10/12_', '_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS' )
chunk_25 <- assemble_3f(25, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/chunks/chunk_25/', '_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS')
chunk_87 <- assemble_3f(87, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/chunks/chunk_87/', '_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS') 
chunk_131 <- assemble_3f(1:10, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/chunks/chunk_131/chunk_131_1-10/131_', '_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS') 
#  SOME CHUNKS ARE MISSING!
# 9 chunks were found!
# 1 chunks were NOT found!
# The missing chunks are the 7

chunk_194 <-  assemble_3f(194, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/chunks/chunk_194/', '_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS') 
chunk_203 <-  assemble_3f(1:10, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/chunks/chunk_203/chunk_203_1-10/203_', '_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS') 


#put all the chunks together and save the output 
all_chunks <- append_chunk(list(chunks_1_335[1:3], chunk_12[1:3], chunk_25[1:3], chunk_87[1:3], chunk_131[1:3], chunk_194[1:3], chunk_203[1:3]))

saveRDS(all_chunks, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/factors_summary_stats.RDS')

#all_chunks <- readRDS( 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/factors_summary_stats.RDS')

#-------- plot the results -----------------------------------------------------
F1 <- all_chunks$sumstats_F1
F2 <- all_chunks$sumstats_F2
F3 <- all_chunks$sumstats_F3
#prepare F1 fro plotting
F1_noNA <- F1[(which(!is.na(F1$Pval_Estimate))),]
dim(F1_noNA)
F1_plot <- F1_noNA[F1_noNA$Pval_Estimate<0.005,]

#prepare F2 for plotting
F2_noNA <- F2[(which(!is.na(F2$Pval_Estimate))),]
dim(F2_noNA)
F2_plot <- F2_noNA[F2_noNA$Pval_Estimate<0.005,]

#prepare F3 for plotting
F3_noNA <- F3[(which(!is.na(F3$Pval_Estimate))),]
dim(F2_noNA)
F3_plot <- F3_noNA[F3_noNA$Pval_Estimate<0.005,]



jpeg(file='outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/miamiF1_F2_F3.jpeg', width = 800, height = 400, units='mm', res = 600)
par(mfrow=c(3,1))
par(mar=c(5,5,3,3))
manhattan(F1_plot, chr="CHR", bp="BP", snp="SNP", p="Pval_Estimate" ,ylim=c(0,30), main='F1 =~ crohn + uc  + psc', 
          col = c("darksalmon", "darkseagreen4"))  
par(mar=c(5,5,3,3))
manhattan(F2_plot, chr="CHR", bp="BP", snp="SNP", p="Pval_Estimate" ,ylim=c(0,30), main='F2 =~ jia + sle + ra + t1d', 
          col = c("darksalmon", "darkseagreen4"))
par(mar=c(5,5,3,3))
manhattan(F3_plot, chr="CHR", bp="BP", snp="SNP", p="Pval_Estimate" ,ylim=c(0,30), main='F3 =~ asthma + derma', 
          col = c("darksalmon", "darkseagreen4"))
dev.off()

#-------------------------------------------------------------------------------
#Q statistics calculation for Factor dependence of the SNPs

#check the chunks that are there and the ones that must be re-estimated 

chunks_found <- lapply(c(1:335), function(x)file.exists(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/qindex/', x,'_qindex_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS' ))) 
which(unlist(chunks_found)==F)# 53 and 93 were not estimated.  

#Error in serverSocket(port = port) : creation of server socket failed: port 11274 cannot be opened
#Error in serverSocket(port = port) : creation of server socket failed: port 11304 cannot be opened
#assmeble the Q statistic summary stats



