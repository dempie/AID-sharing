#  LD score regression on auotoimmunity GWAS
##  in this script I will assemble the chuks of the GWAS that has been produced
#in /sbatch_scripts/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma

#The tutorial and info on the package and how to run the code are here:  
#https://github.com/GenomicSEM/GenomicSEM/wiki/5.-User-Specified-Models-with-SNP-Effects

library(data.table)
library(qqman)
library(ggplot2)

#------ check what chunks are missing ------------------------------------------

chunks_found <- lapply(c(1:335), function(x)file.exists(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/chunks/', x,'_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS' ))) 
#which are not there
chunks_found <- which(chunks_found==F) #12  25  87 131 194 203 chunks are not there, inspect the log files. 

#chunk 12 : task 144 failed - "system is computationally singular: reciprocal condition number = 3.33818e-38"
#chunk25 : Error in unserialize(socklist[[n]]) : error reading from connection. for this just sent the job again
#chunk 87 : Error in unserialize(socklist[[n]]) : error reading from connection. for this just sent the job again
#chunk 131: task 211 failed - "system is computationally singular: reciprocal condition number = 1.50212e-24"
#chunk 194: Error in unserialize(socklist[[n]]) : error reading from connection. for this just sent the job again
#chunk 203:  task 148 failed - "system is computationally singular: reciprocal condition number = 2.86735e-37"


numbers <- rep(1:335)
#remobe chunk 16 and chunk 169
remove <- c(12, 25, 87, 131, 194, 203)
numbers <- numbers[! (numbers %in% remove)]

all_chunks <- lapply(numbers, function(x)readRDS(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/chunks/', x,'_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS' )))


#for each chunk, separate F1 from F2
all_F1 <- list() #allocate the list
all_F2 <- list() #allocate the list
all_F3 <- list() #allocate the list
#create a list of lists, in which each element of the lsit is a chunk
for(i in (1:length(all_chunks))){
  all_F1[[i]] <- all_chunks[[i]][[1]]
  all_F2[[i]] <- all_chunks[[i]][[2]]
  all_F3[[i]] <- all_chunks[[i]][[3]]
}

#create the dataframe for F1 and F2
F1 <-do.call(rbind, all_F1) 
F2 <- do.call(rbind, all_F2)
F3 <- do.call(rbind, all_F3)



length(unique(F1$SNP)) #3 284 158
length(unique(F2$SNP)) #3 284 158
length(unique(F3$SNP))# 3 284 158

head(F1)
head(F2)
head(F3)

unique(F1$lhs)
unique(F2$lhs)
unique(F3$lhs)

nrow(F1[F1$error != 0,]) #46917 SNP with an error 
nrow(F2[F2$error != 0,]) #46917 SNP with an error 
nrow(F3[F3$error != 0,]) #46917 SNP with an error 

nrow(F1[F1$warning != 0,]) #46917 SNP with an error
nrow(F2[F2$warning != 0,]) ##46917 SNP with an error
nrow(F3[F3$warning != 0,]) #46917 SNP with an error


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


pdf(file = '/project/aid_sharing/AID_sharing/F1_manhattan.pdf',width = 4, height = 4 )
#png(filename = '/project/aid_sharing/AID_sharing/F1_manhattan.png', width = 4, height = 4  )
manhattan(F2_plot, chr="CHR", bp="BP", snp="SNP", p="Pval_Estimate") + ggtitle('due')
dev.off()

#-------------------------------------------------------------------------------
#Q statistics calculation for Factor dependence of the SNPs
#assmeble the Q statistic summary stats

assemblechunks <- function(path, n_chunks, chunk_name)
  chunks <- lapply(c(1:n_chunks), function(x)readRDS(paste0('outputs/version3/05_GWAS_results/GWAS_05-05-2022/q_index/chunks_qindex/', x, chunk_name ))) 
  q_summary_stats <- do.call(rbind, chunksQ)


?manhattan
  

