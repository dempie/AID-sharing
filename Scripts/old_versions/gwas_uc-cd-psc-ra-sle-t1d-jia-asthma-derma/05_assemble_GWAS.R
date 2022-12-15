#  LD score regression on auotoimmunity GWAS
##  in this script I will assemble the chuks of the GWAS that has been produced
#in /sbatch_scripts/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma

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
chunks_1_335<- assemble_f(1:335, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/chunks/', '_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS' , 3)

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

chunk_12 <- assemble_f(1:10, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/chunks/chunk_12/chunk_12_1-10/12_', '_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS', n_factors=3 )
chunk_25 <- assemble_f(25, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/chunks/chunk_25/', '_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS', 3)
chunk_87 <- assemble_f(87, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/chunks/chunk_87/', '_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS', 3) 
chunk_131 <- assemble_f(1:10, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/chunks/chunk_131/chunk_131_1-10/131_', '_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS', 3) 
#  SOME CHUNKS ARE MISSING!
# 9 chunks were found!
# 1 chunks were NOT found!
# The missing chunks are the 7

chunk_194 <-  assemble_f(194, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/chunks/chunk_194/', '_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS', 3) 
chunk_203 <-  assemble_f(1:10, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/chunks/chunk_203/chunk_203_1-10/203_', '_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS', 3) 


#put all the chunks together and save the output 
all_chunks <- append_chunk(list(chunks_1_335, chunk_12, chunk_25, chunk_87, chunk_131, chunk_194, chunk_203), 3) #The number of unique SNPs in F123 is 3343158

saveRDS(all_chunks, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/factors_summary_stats.RDS')

all_chunks <- readRDS( 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/factors_summary_stats.RDS')

#-------- plot the results -----------------------------------------------------
F1 <- all_chunks$factor1
F2 <- all_chunks$factor2
F3 <- all_chunks$factor3
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


#save the factor gwas individually for fuma 
for( i in c(1:length( list(F1_noNA, F2_noNA, F3_noNA)))){
  to_save <- list(F1_noNA, F2_noNA, F3_noNA)[[i]]
  to_save <- to_save[,c('SNP', 'CHR', 'BP', 'A1', 'A2', 'est', 'SE', 'Pval_Estimate')]
  to_save <- rename(to_save, Beta ='est', P = 'Pval_Estimate')
  fwrite(to_save, paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/gwas_f',i,'_only_fuma.txt'), 
         col.names = T, row.names = F, sep = '\t', quote = F)
  rm(to_save)
}


#output form

#------- calculate effective sample size of the individual factors -------------


#restrict to MAF of 40% and 10%
N_hat <- vector()
for( i in c(1:length( list(F1_noNA, F2_noNA, F3_noNA)))){
      to_use <- list(F1_noNA, F2_noNA, F3_noNA)[[i]]
      maf_filtered <-subset(to_use , to_use$MAF <= .4 & to_use$MAF >= .1)
      N_hat[i]<-mean(1/((2*to_use$MAF*(1-to_use$MAF))*to_use$SE^2))
      names(N_hat)[i] <- paste0('N_hat F', i)
}

saveRDS(N_hat,'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/effective_sample_size_factors.RDS')


#---plot the manatthan plot 
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

q_1_335 <- assemble_f(1:335, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/qindex/', '_qindex_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS', 1) #The number of unique SNPs is 3324158

#Error in serverSocket(port = port) : creation of server socket failed: port 11274 cannot be opened
#Error in serverSocket(port = port) : creation of server socket failed: port 11304 cannot be opened
#assmeble the Q statistic summary stats

q_53 <- assemble_f(53, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/qindex/chunk_53/', '_qindex_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS', 1)
q_93 <- assemble_f(93, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/qindex/chunk_93/', '_qindex_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS', 1)

#assemble evetyhing together
qindex <- append_chunk(list(q_1_335, q_53, q_93),1) #3344158

saveRDS(qindex, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/qindex_summary_stats.RDS')
qindex <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/qindex_summary_stats.RDS')


#---------calculate the proper qindex-------------------------------------------
#there is a difference of 1000 SNPs between the qindex and the GWAS because 1000 SNP could not be estimated for singulatiry of the matrix in the gwas.
#first create a qindex sumstats that matches the SNPs in sumstats

qindex_f1 <- qindex[[1]]
f1 <- all_chunks[[1]]

dim(qindex_f1) #3 344 158      21
dim(f1) #3 343 158      21

qindex_to_remove <- chmatch(qindex_f1$SNP, f1$SNP, nomatch = NA)  

qindex_ok <- qindex_f1[ which(!is.na(qindex_to_remove)), ]
dim(qindex_ok) #3 343 158 , 1000 less than the original qindex


length(qindex_ok$SNP)==length(f1$SNP) #TRUE

#are there any NA left?
sum(is.na(chmatch(qindex_ok$SNP, f1$SNP, nomatch = NA)))  #0


#chmatch returns a vector of the positions of (first) matches of its first argument in its second.
#reorder qindex as the order of rows in f1
qindex_ok <- qindex_ok[ chmatch( f1$SNP, qindex_ok$SNP,nomatch  = NA) , ]
dim(qindex_ok) #3343158

#check if the order is the same between qindex and SNP
sum(qindex_ok$SNP == f1$SNP) #3 343 158

#calculate Q index
Q_chisq <- f1$chisq - qindex_ok$chisq
Q_df <- f1$chisq_df - qindex_ok$chisq_df
Q_chisq_pval <- pchisq(Q_chisq,Q_df,lower.tail=FALSE)

#add the columns to the sumstats
for(i in c(1:3)){
  all_chunks[[i]]$Q_chisq <- Q_chisq
  all_chunks[[i]]$Q_df <- Q_df
  all_chunks[[i]]$Q_chisq_pval <- Q_chisq_pval 
}


saveRDS(all_chunks, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/gwas_final_withQindex.RDS')

#save FACTOR GWAS individually 
for(i in c(1:3)){
  fwrite(all_chunks[[i]], paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/factor',i,'_gwas_final_withQindex.txt' ),
         sep = '\t', col.names = T, row.names = F, quote = F)
}




#QSNP is a χ2-distributed test statistic, with larger values indexing a violation of the null hypothesis that the SNP acts entirely through the common factor(s)


#---------- plotting -----------------------------------------------------------
gw <- list()
g <- list()
q <- list()
for(i in c(1:3)){
  tt <- c('f1', 'f2', 'f3')[i]
  gw[[tt]] <- fread(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/factor',c('1', '2', '3')[i],'_gwas_final_withQindex.txt' ),
         data.table = F)
  g[[tt]] <-  gw[[tt]][gw[[tt]]$Pval_Estimate < 0.05 & (!is.na(gw[[tt]]$Pval_Estimate)),   ]
  g[[tt]] <-  gw[[tt]][gw[[tt]]$Pval_Estimate < 0.05 & (!is.na(gw[[tt]]$Pval_Estimate)),   ]
  q[[tt]] <- g[[tt]][g[[tt]]$Q_chisq_pval < 0.05 & (!is.na(g[[tt]]$Q_chisq_pval)),   ]
  
}




#pvalue of the SNP
jpeg(file='outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/manhattan_pvalue.jpeg', width = 800, height = 400, units='mm', res = 600)
par(mfrow=c(3,1))
for(i in 1:length(g)){
  tt <- c('f1', 'f2', 'f3')[i]
          par(mar=c(5,5,3,3))
          manhattan(g[[tt]], chr="CHR", bp="BP", snp="SNP", p="Pval_Estimate" ,ylim=c(0,30), main=paste0(tt,'=~ crohn + uc  + psc'), 
                    col = c("cornflowerblue", "coral2"))  
}
dev.off()


#Q value of the SNP
jpeg(file='outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/manhattan_qpvaule.jpeg', width = 800, height = 400, units='mm', res = 600)
par(mfrow=c(3,1))
for(i in 1:length(g)){
  tt <- c('f1', 'f2', 'f3')[i]
  par(mar=c(5,5,3,3))
  manhattan(q[[tt]], chr="CHR", bp="BP", snp="SNP", p="Q_chisq_pval" ,ylim=c(0,30), main=paste0(tt,'=~ crohn + uc  + psc'), 
            col = c("cornflowerblue", "coral2"))  
}
dev.off()




which(is.na(q[[tt]]$Pval_Estimate))
par(mfrow=c(2,1))
par(mar=c(0,5,3,3))
manhattan(gwasResults,ylim=c(0,10),cex=2.2,cex.lab=2.5,font.lab=2,font.axis=2,cex.axis=1.6,las=2,font=4)
par(mar=c(5,5,3,3))
manhattan(gwasResults,ylim=c(10,0),cex=2.2,cex.lab=2.5,font.lab=2,font.axis=2,cex.axis=1.6,las=2,font=4,xlab="",xaxt="n")
dev.off()



#--------circular plot---------------------------------------------------------
library("CMplot")


gw <- list()
g <- list()
q <- list()
for(i in c(1:3)){
  tt <- c('f1', 'f2', 'f3')[i]
  gw[[tt]] <- fread(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/factor',c('1', '2', '3')[i],'_gwas_final_withQindex.txt' ),
                    data.table = F)
  g[[tt]] <-  gw[[tt]][gw[[tt]]$Pval_Estimate < 0.5 & (!is.na(gw[[tt]]$Pval_Estimate)),   ]
  g[[tt]] <-  gw[[tt]][gw[[tt]]$Pval_Estimate < 0.05 & (!is.na(gw[[tt]]$Pval_Estimate)),   ]
}


SNPs <- unique(c(gw[['f1']]$SNP, gw[['f2']]$SNP, gw[['f3']]$SNP ))
circular <- data.frame('SNP'=SNPs)

circular[,c('CHR')] <- gw[['f1']][match(circular$SNP,gw[['f1']]$SNP ),]$CHR
circular[,c('BP')] <- gw[['f1']][match(circular$SNP,gw[['f1']]$SNP ),]$BP

circular[, 'P_f1'] <- gw[['f1']][match(circular$SNP,gw[['f1']]$SNP ),]$Pval_Estimate
circular[, 'P_f2'] <- gw[['f2']][match(circular$SNP,gw[['f2']]$SNP ),]$Pval_Estimate
circular[, 'P_f3'] <- gw[['f3']][match(circular$SNP,gw[['f3']]$SNP ),]$Pval_Estimate

CMplot(circular,type="p",plot.type="c",chr.labels=paste("Chr",c(1:22),sep=""),r=4,cir.legend=TRUE,
       outward=T,cir.legend.col="black",cir.chr.h=1.3,chr.den.col="black",file="jpg",
       memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=20,height=20, amplify = T, cir.band = 1, threshold=5e-10, band=0.2)




CMplot(circular[1:20000], plot.type="m",multracks=TRUE,threshold=c(5e-10),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","green","blue"),
       signal.cex=1, file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE, LOG10=T,  width=40, height=5)

#browseVignettes("karyoploteR")

