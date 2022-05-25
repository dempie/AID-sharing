# in this script all the necessary steps for the colocalizzation will be performed 

library(GenomicSEM)
library(data.table)
loci_all_gwas <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/06_gwas_analysis/loci_all_gwas.RDS')

cd <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/psc_ji-2016.txt', data.table = F)


fwrite(head(cd, 10000), 'cd_test.txt', sep = '\t', col.names = T, row.names = F, quote = F)
fwrite(head(cd, 10000), 'cd_2_test.txt', sep = '\t', col.names = T, row.names = F, quote = F)

cc<- sumstats(files = c('cd_test.txt', 'cd_2_test.txt'), ref = 'SNP/reference.1000G.maf.0.005.txt.gz',
         trait.names = c('c1', 'c2') ,
         se.logit =  c(F,F),  
         OLS = c(T,T), 
         linprob = c(F,F), 
         N = c(NA, NA), 
         parallel = T,
         keep.indel= F 
)

c1 <- fread('cd_test.txt')
cc$beta.c1 == 

head(cc$beta.c1 )

y <- cc
x <- (c1[(c1$SNP %in% cc$SNP),  ])

x <- x[order(x$CHR, x$BP),]
y <- y[order(y$CHR, y$BP  ) ,]
 

x[which(!(x$Beta==y$beta.c1)),]
y[which(!(x$Beta==y$beta.c1)),]

mean(abs(x$Beta)==abs(y$beta.c1))
mean(x$SE==y$se.c1)













