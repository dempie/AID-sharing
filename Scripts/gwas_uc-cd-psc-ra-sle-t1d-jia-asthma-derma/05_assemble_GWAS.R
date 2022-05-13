#  LD score regression on auotoimmunity GWAS
##  in this script I will assemble the chuks of the GWAS that has been produced
#in /sbatch_scripts/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma

#The tutorial and info on the package and how to run the code are here:  
#https://github.com/GenomicSEM/GenomicSEM/wiki/5.-User-Specified-Models-with-SNP-Effects

library(data.table)
library(qqman)
getwd()




numbers <- rep(1:432)
#remobe chunk 16 and chunk 169
remove <- c(16, 169)
numbers <- numbers[! (numbers %in% remove)]

all_chunks <- lapply(numbers, function(x)readRDS(paste0(x, '_gwas_twof.RDS')))

#for each chunk, separate F1 from F2
all_F1 <- list() #allocate the list
all_F2 <- list() #allocate the list
#create a list of lists, in which each element of the lsit is a chunk
for(i in (1:length(all_chunks))){
  all_F1[[i]] <- all_chunks[[i]][[1]]
  all_F2[[i]] <- all_chunks[[i]][[2]]
}

#create the dataframe for F1 and F2
F1 <-do.call(rbind, all_F1) 
F2 <- do.call(rbind, all_F2)



length(unique(F1$SNP)) #4291927
length(unique(F2$SNP)) #4291927

head(F1)
head(F2)

unique(F1$lhs)
unique(F2$lhs)

nrow(F1[F1$error != 0,]) #27363 SNP with an error 
nrow(F2[F2$error != 0,]) #27363 SNP with an error

nrow(F1[F1$warning != 0,]) #27400 SNP with an error 
nrow(F2[F2$warning != 0,]) #27400 SNP with an error

F1_noNA <- F1[(which(!is.na(F1$Pval_Estimate))),]

?manhattan



pdf(file = '/project/aid_sharing/AID_sharing/F1_manhattan.pdf',width = 4, height = 4 )
#png(filename = '/project/aid_sharing/AID_sharing/F1_manhattan.png', width = 4, height = 4  )
manhattan(F1_noNA, chr="CHR", bp="BP", snp="SNP", p="Pval_Estimate", annotateTop = T )
dev.off()

#-------------------------------------------------------------------------------
#Q statistics calculation for Factor dependence of the SNPs
#assmeble the Q statistic summary stats

assemblechunks <- function(path, n_chunks, chunk_name)
  chunks <- lapply(c(1:n_chunks), function(x)readRDS(paste0('outputs/version3/05_GWAS_results/GWAS_05-05-2022/q_index/chunks_qindex/', x, chunk_name ))) 
  q_summary_stats <- do.call(rbind, chunksQ)




