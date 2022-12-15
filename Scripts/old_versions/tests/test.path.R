library(GenomicSEM)
library(data.table)
library(RColorBrewer)

aid_sumstats <- readRDS('/project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/04_sumstat_outputs/sumstats_output/sumstats_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS')
ldsc_model <- readRDS('/project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/04_sumstat_outputs/ldsc_output/ldsc_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS')

aid_model <- 'F1 =~ NA*crohn + uc  + psc  
F2 =~ NA*jia + sle + ra+ t1d 
F3 =~ NA*asthma+ derma 

F1~~F2
F1~~F3
F2~~F3
F1 ~~ 1*F1
F2 ~~ 1*F2
F3 ~~ 1*F3

F1 ~ SNP
F2 ~ SNP
F3 ~ SNP

derma~~a*derma
a>0.001'


#dim(aid_sumstats) 3344158      24
#how_many_chunks(aid_sumstats, 10000) #335



output <- userGWAS(covstruc = ldsc_model, 
                   SNPs = aid_sumstats[1,], 
                   model = aid_model)

out <- list()
out$results <- output[[1]]


semPaths(a )

aid_factor$results
out

a <- semPlotModel_GSEM(out, est.label = 'est')

semPaths(a, what = 'path' , 
         whatLabels= 'est',
         residuals = T, 
         sizeMan = 7, 
         sizeInt = 2,
         label.cex=1, 
         theme="colorblind", 
         rotation = 4, 
         layout = "tree2", 
         sizeLat = 7,
         edge.color = "black",
         edge.label.cex = 1,
         curve = 2,
         nodeLabels = c('CD', 'UC', 'PSC', 'JIA', 
                        'SLE', 'RA','T1D' ,'AST', 'ECZ' , 'SNP', 'F1', 'F2', 'F3'), 
         label.cex = 1.5, 
         color = list( lat=brewer.pal(12, 'Paired')[c(1,3,5)]),
         groups=list(c("crohn","uc","psc"),c("jia","sle","ra_eu","t1d"), c( "asthma" ,"derma")),
         height = 10, width = 16, 
         edge.label.position=0.5,
         asize=3,
         esize=1)
