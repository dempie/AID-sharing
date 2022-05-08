#  LD score regression on auotoimmunity GWAS
##  in this script I will assemble the chuks of the GWAS that has been produced
#in sbatch_scripts

#The tutorial and info on the package and how to run the code are here:  
#https://github.com/GenomicSEM/GenomicSEM/wiki/5.-User-Specified-Models-with-SNP-Effects

library(data.table)
library(qqman)
getwd()
setwd('outputs/version3/05_GWAS_results/GWAS_05-05-2022/chunks/')



numbers <- rep(1:432)
#remobe chunk 16 and chunk 169
remove <- c(16, 169)
numbers <- numbers[! (numbers %in% remove)]
#-------------------------------------------------------------------------------

  
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

F1_NA <- F1[(which(!is.na(F1$Pval_Estimate))),]
F2_NA <- F2[(which(!is.na(F2$Pval_Estimate))),]

#miami plot of the SNPs with p<0.005
F1_filtered <- F1_NA[F1_NA$Pval_Estimate < (0.005), ]
F2_filtered <- F2_NA[F2_NA$Pval_Estimate < (0.005), ]


jpeg('outputs/version3/05_GWAS_results/GWAS_05-05-2022/plots/plotF1_F2.jpeg', units="in", width=14, height=7, res=300)
par(mfrow=c(2,1))
par(mar=c(0,5,3,3))
manhattan(F1_filtered, chr="CHR", bp="BP", snp="SNP", p="Pval_Estimate" )
par(mar=c(5,5,3,3))
manhattan(F2_filtered, chr="CHR", bp="BP", snp="SNP", p="Pval_Estimate" ,ylim=c(30,0), xlab="",xaxt="n")
dev.off()







