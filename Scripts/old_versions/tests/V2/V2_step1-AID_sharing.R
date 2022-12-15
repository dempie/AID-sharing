#  LD score regression on auotoimmunity GWAS
## Version2 will explor more traits and go into more details
 
library(data.table)
library(GenomicSEM)
library(tidyr)
library(dplyr)
library(corrplot)

#-----------------------------------------------------------------------------

# load in the summary stats
arm_fat <- fread('Summary_Stats/continuous-23123-both_sexes-irnt.tsv.gz', data.table = F) 
head(arm_fat)

auto_thyro <- fread('Summary_Stats/saevarsdottir_auto-thyroid_build37_AITD2020', data.table = F)
head(auto_thyro)

sle <- fread('Summary_Stats/bentham-2015_sle_build37_26502338_sle_efo0002690_1_gwas.sumstats.tsv.gz', data.table = F)
head(sle)

celiac <- fread('Summary_Stats/dubois-2010_celiac_build37_20190752_cel_efo0001060_1_gwas.sumstats.tsv.gz', data.table = F)
head(celiac)

pbc <- fread('Summary_Stats/cordell-2015_pbc_build37_26394269_pbc_efo1001486_1_gwas.sumstats.tsv.gz', data.table = F)
head(pbc)

psc <- fread('Summary_Stats/ji-2016_psc_build37_ipscsg2016.result.combined.full.with_header.txt', data.table = T)
head(psc)

ms <- fread('Summary_Stats/andlauer-2011_ms_build36_imsgc_2011_21833088_ms_efo0003885_1_gwas.sumstats.tsv.gz', data.table = F)
head(ms)
#----------------------------------------------------------------------------
#for arm fat

#the column names must be A1 and A2 for genomicSEM
colnames(arm_fat)[3] #"ref"
colnames(arm_fat)[3] <- 'A2'

colnames(arm_fat)[4] #alt
colnames(arm_fat)[4] <- 'A1'

colnames(arm_fat)[8] # pval
colnames(arm_fat)[8] <- 'P'


colnames(arm_fat)[6] #'"beta_meta"'
colnames(arm_fat)[6] <- 'Beta'

# add the rsID to run ldsc
#load in and prepare the reference file 
referenceSNP <- fread('SNP/reference.1000G.maf.0.005.txt.gz')
referenceSNP <- referenceSNP %>% select(-c('MAF', 'A1', 'A2')) %>% unite(CHR, BP, sep= ':', na.rm = F, remove = T, col = 'chrPosition')
head(referenceSNP)

#prepare the file of arm_fat
arm_fat <- unite(arm_fat, chr, pos, na.rm = F, remove = T, col = 'chrPosition' , sep = ':')

#add rsID
arm_fat_rsID <- merge(x=arm_fat, 
                      y=referenceSNP,
                      by.x= 'chrPosition', 
                      by.y= 'chrPosition', 
                      all.x=T, 
                      all.y=F, 
                      sort=T
                      )  

arm_fat_rsID_NA <- arm_fat_rsID[which(!is.na(arm_fat_rsID$SNP)), ]


#------------------------------------------------------------------------------
#auto_thyr
#for auto_thyroid there are a lot of snps without rsID
dim(auto_thyro) #44690176       13
length(grep('rs', auto_thyro$rsID))

colnames(auto_thyro)

colnames(auto_thyro)[4] #"A0"
colnames(auto_thyro)[4] <- 'A2'

colnames(auto_thyro)[5] #A1"
colnames(auto_thyro)[5] <- 'A1'

colnames(auto_thyro)[10] #"OR-A1"
colnames(auto_thyro)[10] <- 'EFFECT'

#munge function does not read the P column 
#let's create a new matrix 

auto_thyro_ok <- data.frame(SNP= auto_thyro$rsID, 
                       A1 = auto_thyro$A1, 
                       A2= auto_thyro$A2, 
                       P = auto_thyro$P, 
                       EFFECT = auto_thyro$EFFECT, 
                       Pos = auto_thyro$Pos)

#------------------------------------------------------------------------------
#sle 
colnames(sle)
colnames(sle)[4] #"other_allele"
colnames(sle)[4] <- 'A2'

colnames(sle)[5] # effect_allele"
colnames(sle)[5] <- 'A1'

#remove OR and keep only beta
sle <- select(sle, -c('OR' , 'OR_lower', 'OR_upper'))

#-----------------------------------------------------------------------------
#celiac 
colnames(celiac)
colnames(celiac)[4] #"other_allele"
colnames(celiac)[4] <- 'A2'

colnames(celiac)[5] #"effect_allele"
colnames(celiac)[5] <- 'A1'

#remove OR and keep only beta
celiac <- select(celiac, -c('OR' , 'OR_lower', 'OR_upper'))


dim(celiac) #523398     11

#---------------------------------------------------------------------------
#pbc 
colnames(pbc)
colnames(pbc)[4] #"other_allele"
colnames(pbc)[4] <- 'A2'

colnames(pbc)[5] #"effect_allele"
colnames(pbc)[5] <- 'A1'

#remove OR and keep only beta
pbc <- select(pbc, -c('OR' , 'OR_lower', 'OR_upper'))

dim(pbc) # 1134141      11

#--------------------------------------------------------------------------
#psc
head(psc)

#munge does not accept psc format for some reason I do not know, 
# I try to make a matrix

psc_ok <- data.frame( SNP = psc$SNP , 
                A2 = psc$allele_0  , 
                A1 = psc$allele_1 , 
                Pos = psc$pos, 
                Effect = psc$or , 
                SE = psc$se ,
                P = psc$p )

dim(psc) #7891602      13
#--------------------------------------------------------------------------
#ms
head(ms)

colnames(ms)

colnames(ms)[4] #"other_allele"
colnames(ms)[4] <- 'A2'

colnames(ms)[5] #"effect_allele"
colnames(ms)[5] <- 'A1'

#remove OR and keep only beta
ms <- select(ms, -c('OR' , 'OR_lower', 'OR_upper'))


dim(ms) #472086     11

#--------------------------------------------------------------------------
#save the files with fwrite that is very fast
data.table::fwrite(arm_fat_rsID, 'Outputs/Version2/Step1/Sumstats/arm_fat-l_rsID_noNA.txt', 
            col.names = T, row.names = F, sep = '\t', quote = F)

data.table::fwrite(auto_thyro_ok, 'Outputs/Version2/Step1/Sumstats/auto_thyro_ok_rsID.txt',
            col.names= T, row.names =F, sep = '\t', quote = F)

data.table::fwrite(sle, 'Outputs/Version2/Step1/Sumstats/sle_rsID.txt', sep = '\t', 
            col.names = T, row.names = F, quote = F)

data.table::fwrite(celiac, 'Outputs/Version2/Step1/Sumstats/celiac_rsID.txt', sep='\t',
            col.names = T, row.names = F, quote = F)

data.table::fwrite(pbc, 'Outputs/Version2/Step1/Sumstats/pbc_rsID.txt', sep= '\t', 
            col.names = T, row.names = F, quote = F)

data.table::fwrite(psc_ok, 'Outputs/Version2/Step1/Sumstats/psc_ok_rsID.txt', sep = '\t',
            col.names = T, row.names = F, quote = F)

data.table::fwrite(ms, 'Outputs/Version2/Step1/Sumstats/ms_rsID.txt', sep='\t', 
            col.names=T, row.names=F, quote=F)


#----------------------------------------------------------------------------
#munge the data
vector_of_traits <- c('Outputs/Version2/Step1/Sumstats/arm_fat-l_rsID_noNA.txt',
                      'Outputs/Version2/Step1/Sumstats/auto_thyro_ok_rsID.txt',
                      'Outputs/Version2/Step1/Sumstats/sle_rsID.txt',
                      'Outputs/Version2/Step1/Sumstats/celiac_rsID.txt',
                      'Outputs/Version2/Step1/Sumstats/pbc_rsID.txt',
                      'Outputs/Version2/Step1/Sumstats/psc_ok_rsID.txt',
                      'Outputs/Version2/Step1/Sumstats/ms_rsID.txt'
                      )
munge(vector_of_traits, N=c(492874,755406 , 14267, 15283, 13239, 14890, 26621), 
      trait.names = c('arm_fat', 'thyroid ','sle', 'celiac', 'pbc', 'psc', 'ms'),  
      hm3 = 'SNP/w_hm3.snplist' )


#-----------------------------------------------------------------------------


#ldsc
traits <- c('Outputs/Version1/V1_step3/Munged_files/crohn.sumstats.gz', 
            'Outputs/Version1/V1_step3/Munged_files/ulccol.sumstats.gz',
            'Outputs/Version1/V1_step3/Munged_files/allergies.sumstats.gz',
            'Outputs/Version1/V1_step3/Munged_files/alzheimer.sumstats.gz', 
            'Outputs/Version1/V1_step3/Munged_files/ra.sumstats.gz', 
            'Outputs/Version1/V1_step3/Munged_files/asthma_euro.sumstats.gz',
            
            'Outputs/Version2/Step1/Munged_files/arm_fat.sumstats.gz',
            'Outputs/Version2/Step1/Munged_files/thyroid.sumstats.gz',
            'Outputs/Version2/Step1/Munged_files/sle.sumstats.gz',
            'Outputs/Version2/Step1/Munged_files/celiac.sumstats.gz',
            'Outputs/Version2/Step1/Munged_files/pbc.sumstats.gz',
            'Outputs/Version2/Step1/Munged_files/psc.sumstats.gz',
            'Outputs/Version2/Step1/Munged_files/ms.sumstats.gz') 

trait.names <- c( 'crohn', 'ulccol', 'allergies' ,'alzheimer',  'ra', 'asthma',  
                  'arm_fat', 'auto_thyro', 'sle', 'celiac', 'pbc', 'psc', 'ms'
                  )

sample.prev <- round(c(12194/40266, 12366/45975, 180129/360838 , 21982/63926 , 19234/103638, 19954/107715, 
                       NA, (30234/725172), (5201/9066), (4533/10750), (2764/10475 ), (2871/12019), (9772/16849)), 2)

population.prev <- round(c((100/100000),(30/100000),(0.20),(0.059),(460/100000), 0.0357, 
                           NA,(0.05), (50/100000), (0.014), (10/100000), (5/100000), (35.9/100000) ), 5)

ld <- "ldscores/eur_w_ld_chr"

wld <- "ldscores/eur_w_ld_chr"

LDS_output <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names, stand = T)

saveRDS(LDS_output, file = 'Outputs/Version2/Step1/LDS_output')
output <- readRDS('Outputs/Version2/Step1/LDS_output')

diag(output$S)
rownames(output$S_Stand) <- colnames(output$S_Stand)
corrplot(output$S_Stand, order = 'hclust', addCoef.col = 'black', is.corr = F,)








