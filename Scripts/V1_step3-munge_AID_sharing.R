# LD score regression on auotoimmunity GWAS
## Step 3: Munge the summary statistics
#The tutorial and info on the package and how to run the code are here:  
#https://github.com/GenomicSEM/GenomicSEM/wiki/3.-Models-without-Individual-SNP-effects

library(GenomicSEM)
library(data.table)
library(pheatmap)
library(corrplot)

#munge
file_vector <- c('Outputs/Version1/V1_step2/crohn_delange_s2.txt', 
                 'Outputs/Version1/V1_step2/ulccol_delange_s2.txt',
                 'Outputs/Version1/V1_step2/allergies_ferreira_s2.txt', 
                 'Outputs/Version1/V1_step2/alzheimer_kunkle_s2.txt',
                 'Outputs/Version1/V1_step2/ssc_lopez_s2.txt',
                 'Outputs/Version1/V1_step2/reart_okada_s2.txt',
                 'Outputs/Version1/V1_step2/asthma_demeanis_euro_s2.txt')
munge(file_vector,
      hm3 = 'SNP/w_hm3.snplist',
      trait.names = c('crohn', 'ulccol', 'allergies' ,'alzheimer', 'ssc',  'ra', 'asthma_euro'),
      N= c( 40266, 45975, 360838 , 63926 , 26679 , 103638 , 127669 )) 

#the munged files are in Outuputs/Version1/V1_step3
#ssc_lopez does not report A2, or at least they are not read

#read the munged files and inspect them 
allergies <- fread('Outputs/Version1/V1_step3/Munged_files/allergies.sumstats.gz', data.table = F)
alzheimer <- fread('Outputs/Version1/V1_step3/Munged_files/alzheimer.sumstats.gz', data.table = F)
asthma_euro <- fread('Outputs/Version1/V1_step3/Munged_files/asthma_euro.sumstats.gz', data.table = F)
crohn <- fread('Outputs/Version1/V1_step3/Munged_files/crohn.sumstats.gz', data.table = F)
ra <- fread('Outputs/Version1/V1_step3/Munged_files/ra.sumstats.gz', data.table = F)
ulcol <- fread('Outputs/Version1/V1_step3/Munged_files/ulccol.sumstats.gz', data.table = F)

#alleregies
head(allergies)
dim(allergies) # 1174535       5
#alzheimer 
head(alzheimer)
dim(alzheimer) #1207073       5
#asthma
head(asthma_euro)
dim(asthma_euro) #1045391       5
#crohn
head(crohn)
dim(crohn) #1152121       5
#ulcol
head(ulcol)
dim(ulcol) # 1152245       5

#---------------------------------------------------------------------------

#ldsc function 

##ldsc with asthma_euro
traits <- c('Outputs/Version1/V1_step3/Munged_files/crohn.sumstats.gz', 
            'Outputs/Version1/V1_step3/Munged_files/ulccol.sumstats.gz',
            'Outputs/Version1/V1_step3/Munged_files/allergies.sumstats.gz',
            'Outputs/Version1/V1_step3/Munged_files/alzheimer.sumstats.gz', 
            'Outputs/Version1/V1_step3/Munged_files/ra.sumstats.gz', 
            'Outputs/Version1/V1_step3/Munged_files/asthma_euro.sumstats.gz') 
trait.names <- c('crohn', 'ulccol', 'allergies' ,'alzheimer',  'reart', 'asthma')

sample.prev <- round(c( 12194/40266, 12366/45975, 180129/360838 , 21982/63926 , 19234/103638, 19954/107715 ), 2)

population.prev <- round(c((100/100000),(30/100000),(0.20),(0.059),(460/100000), 0.0357), 4)
 # poluation prevalence calculated 0.0010 0.0003 0.2000 0.0590 0.0046 0.0357
ld <- "ldscores/eur_w_ld_chr"
wld <- "ldscores/eur_w_ld_chr"
LDS_version1 <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, stand = T)
saveRDS(LDS_version1, file = 'Outputs/Version1/V1_step3/LDscores_V1')

#LDS_version1 <- readRDS('Outputs/Version1/V1_step3/LDscores_V1')

LDS_version1$S_Stand
cormatrix <-  LDS_version1$S_Stand
rownames(cormatrix) <- c('crohn', 'ulccol', 'allergies' ,'alzheimer',  'reart', 'asthma')
pheatmap(cormatrix)
a <- corrplot(cormatrix, order = 'hclust', addCoef.col = 'black')

?png()

#standard errors
k<-nrow(LDS_version1$S)
SE_corr <-matrix(0, k, k)
SE_corr[lower.tri(SE,diag=TRUE)] <-sqrt(diag(LDS_version1$V) #the SE of the correlation


#the heritability lies on the diagonal of the genetic covariance matrix 
diag(LDS_version1$S)
#---------------------------------------------------------------------------

















