# LD score regression on auotoimmunity GWAS
## Step 2: check the columns, as for munge fucntion 
#The tutorial and info on the package and how to run the code are here:  
#https://github.com/GenomicSEM/GenomicSEM/wiki/3.-Models-without-Individual-SNP-effects

library(GenomicSEM)
library(data.table)
library(dplyr)

#read the data
crohn <- fread('Outputs/Version1/V1_step1/crohn_delange_rsID_noNA.txt', data.table = F)
ulcol <-fread('Outputs/Version1/V1_step1/ulccol_delange_rsID_noNA.txt', data.table = F)
alzh <- fread('Outputs/Version1/V1_step1/alzheimer_kunkle_rsID.txt', data.table = F)
allerg <- fread('Outputs/Version1/V1_step1/allergies_ferreira_rsID.txt', data.table = F)
reart <- fread('Outputs/Version1/V1_step1/ra_okada_rsID.txt', data.table = F)
sscler <- fread('Outputs/Version1/V1_step1/ssc_lopez_rsID.txt', data.table = F)
asthma <- fread('Outputs/Version1/V1_step1/asthma_demeanis_rsID.txt', data.table = F)
#alzh_wi <- fread('Outputs/Version1/V1_step1/alzheimer_wightman_rsID_noNA.txt', data.table = F)

## remove useless columns (Extra_1 and Extra_2) and 
## rearrange the columns for easy comparison


# crohn
head(crohn)
dim(crohn)  # 8361641      20
crohn <- crohn %>% select(-c(Extra_1, Extra_2, A1, A2, HetISq, HetChiSq, HetDf, HetPVal, Pval_IBDseq, Pval_IIBDGC, Pval_GWAS3, Min_single_cohort_pval  )) 
dim(crohn) # 8361641      8
head(crohn)

# ulcol
dim(ulcol) #  8375870      20
ulcol <- ulcol %>% select(- c(Extra_1, Extra_2, A1, A2, HetISq, HetChiSq, HetDf, HetPVal, Pval_IBDseq, Pval_IIBDGC, Pval_GWAS3, Min_single_cohort_pval))
dim(ulcol) # 8375870      8

#asthma
head(asthma)
asthma_euro <-  select(asthma, -c( 'Multiancestry_beta_fix', 'Multiancestry_se_fix',
                                   'Multiancestry_pval_fix',"Multiancestry_beta_rand",
                                   "Multiancestry_se_rand",          
                                   "Multiancestry_pval_rand", "Multiancestry_HetQtest",         
                                   "Multiancestry_df_HetQtest","Multiancestry_pval_HetQtest",
                                    "European_ancestry_beta_rand",    
                                   "European_ancestry_se_rand" , "European_ancestry_pval_rand" ,  
                                   "European_ancestry_HetQtest","European_ancestry_df_HetQtest", 
                                   "European_ancestry_pval_HetQtest"))
colnames(asthma_euro)[2] #"rsid"
colnames(asthma_euro)[2] <- 'rsID'

colnames(asthma_euro)[6] #"European_ancestry_beta_fix"
colnames(asthma_euro)[6] <- 'Beta'

colnames(asthma_euro)[8] #"European_ancestry_pval_fix"
colnames(asthma_euro)[8] <- 'P'

colnames(asthma_euro)[4] #"reference_allele" is usually the A2 the non-effect
colnames(asthma_euro)[4] <- 'A2'

colnames(asthma_euro)[5] #"alternate_allele" is the effect Allele A1
colnames(asthma_euro)[5] <- 'A1'

head(asthma_euro)

#allerg
head(allerg)
colnames(allerg)[9] <- 'rsID'

#alzh
head(alzh)
#reart
head(reart)
colnames(reart)[7] <- 'P'
colnames(reart)[1] <- 'rsID' 
colnames(reart)[6] <- 'OR'

#sscler
head(sscler)
sscler <- fread('Outputs/Version1/V1_step2/ssc_lopez_s2.txt')


#save the txt in Version1 > V1_step2
write.table(crohn, file = 'Outputs/Version1/V1_step2/crohn_delange_s2.txt',
            sep = "\t", quote = F, row.names = F, col.names=T )
write.table(ulcol, file = 'Outputs/Version1/V1_step2/ulccol_delange_s2.txt',
            sep = "\t", quote = F, row.names = F, col.names=T )
write.table(allerg, file ='Outputs/Version1/V1_step2/allergies_ferreira_s2.txt',
            sep = "\t", quote = F, row.names = F, col.names=T )
write.table(alzh, file = 'Outputs/Version1/V1_step2/alzheimer_kunkle_s2.txt',
            sep = "\t", quote = F, row.names = F, col.names=T )
write.table(sscler, file= 'Outputs/Version1/V1_step2/ssc_lopez_s2.txt',
            sep = "\t", quote = F, row.names = F, col.names=T )
write.table(reart, file = 'Outputs/Version1/V1_step2/reart_okada_s2.txt',
            sep = "\t", quote = F, row.names = F, col.names=T )
write.table(asthma_euro, file= 'Outputs/Version1/V1_step2/asthma_demeanis_euro_s2.txt',
            sep = "\t", quote = F, row.names = F, col.names=T )









