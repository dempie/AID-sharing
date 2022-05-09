#  LD score regression on auotoimmunity GWAS
##  in this script I will assemble the chuks of the GWAS that has been produced
#in sbatch_scripts

#The tutorial and info on the package and how to run the code are here:  
#https://github.com/GenomicSEM/GenomicSEM/wiki/5.-User-Specified-Models-with-SNP-Effects

library(data.table)
library(qqman)


numbers <- rep(1:432)
#remobe chunk 16 and chunk 169
remove <- c(16, 169)
numbers <- numbers[! (numbers %in% remove)]
#-------------------------------------------------------------------------------

#all chunk exept 16 and 169  
all_chunks <- lapply(numbers, function(x)readRDS(paste0('outputs/version3/05_GWAS_results/GWAS_05-05-2022/chunks/', x, '_gwas_twof.RDS')))
#chunk 16
chunk_16 <- list()
chunk_16[[1]] <- readRDS('outputs/version3/05_GWAS_results/GWAS_05-05-2022/chunk_16/6_gwas_twof.RDS')
#what works of chunk 169
chunk_169 <- list()
chunk_169 <- lapply(c(1, 2, 3, 4, 7, 8, 9, 10), function(x)readRDS(paste0('outputs/version3/05_GWAS_results/GWAS_05-05-2022/chunk_169/chunk_169_1-10/', x, '_gwas_twof.RDS')))

#combine the three list
chunks_complete <- c(all_chunks, chunk_16, chunk_169)
length(chunks_complete) #439

#for each chunk, separate F1 from F2
all_F1 <- list() #allocate the list
all_F2 <- list() #allocate the list
#create a list of lists, in which each element of the lsit is a chunk
for(i in (1:length(chunks_complete))){
  all_F1[[i]] <- chunks_complete[[i]][[1]]
  all_F2[[i]] <- chunks_complete[[i]][[2]]
}

#create the dataframe for F1 and F2
F1 <-do.call(rbind, all_F1) 
F2 <- do.call(rbind, all_F2)

#---- save the files -----------------------------------------------------------

saveRDS(F1,'outputs/version3/05_GWAS_results/F1_completeGWAS.RDS')
saveRDS(F2,'outputs/version3/05_GWAS_results/F2_completeGWAS.RDS')

#------------------------------------------------------------------------------

length(F1$SNP) #4309927
length(F2$SNP) #4309927


length(unique(F1$SNP)) #4309927
length(unique(F2$SNP)) #4309927

head(F1)
head(F2)

unique(F1$lhs) #"F1" NA  
unique(F2$lhs) #"F2" NA  

nrow(F1[F1$error != 0,]) # 28444 SNP with an error 
nrow(F2[F2$error != 0,]) # 28444 SNP with an error

nrow(F1[F1$warning != 0,]) #28574 SNP with an error 
nrow(F2[F2$warning != 0,]) #28574 SNP with an error

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







