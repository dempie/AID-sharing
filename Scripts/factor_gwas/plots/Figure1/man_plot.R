
library(data.table)
library(qqman)
library(dplyr)
library(CMplot)
library(rcartocolor)
library(RColorBrewer)

#manathan plots
gw <- list()
g <- list()
q <- list()
for(i in c(1:3)){
  tt <- c('f1', 'f2', 'f3')[i]
  gw[[tt]] <- fread(paste0('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_munge_for_nicola/munged/',tt,'_munged_build37.txt' ),
                    data.table = F)
  # g[[tt]] <-  gw[[tt]][gw[[tt]]$Pval_Estimate < 0.5 & (!is.na(gw[[tt]]$Pval_Estimate)),   ]
  # g[[tt]] <-  gw[[tt]][gw[[tt]]$Pval_Estimate < 0.05 & (!is.na(gw[[tt]]$Pval_Estimate)),   ]
}


SNPs <- unique(c(gw[['f1']]$SNP, gw[['f2']]$SNP, gw[['f3']]$SNP ))
circular <- data.frame('SNP'=SNPs)

circular[,c('CHR')] <- gw[['f1']][match(circular$SNP,gw[['f1']]$SNP ),]$CHR
circular[,c('BP')] <- gw[['f1']][match(circular$SNP,gw[['f1']]$SNP ),]$BP

circular[, 'f1'] <- gw[['f1']][match(circular$SNP,gw[['f1']]$SNP ),]$P
circular[, 'f2'] <- gw[['f2']][match(circular$SNP,gw[['f2']]$SNP ),]$P
circular[, 'f3'] <- gw[['f3']][match(circular$SNP,gw[['f3']]$SNP ),]$P




factor_loci <- fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/06_genomic_regions_all_traits_munged/factor.regions.txt', data.table = F) 

#v3

matrix(c(rcartocolor::carto_pal(3, 'Safe')[1],rcartocolor::carto_pal(3, 'Safe')[2],rcartocolor::carto_pal(3, 'Safe')[3],
         RColorBrewer::brewer.pal(4, 'Blues')[4], RColorBrewer::brewer.pal(3, 'Reds')[3], RColorBrewer::brewer.pal(4, 'YlOrRd')[2]),
       3,2,byrow=T)

#run twice, one pdf, the other jpg
for(i in 1:2){
  CMplot(circular, 
         plot.type="m",
         multracks=TRUE,
         threshold=c(5e-8),
         threshold.lty=c(1,2), 
         threshold.lwd=c(1,1), 
         threshold.col=c("black","grey"), 
         amplify=F,
         #signal.col=list(brewer.pal(12, 'Paired')[c(1)], brewer.pal(12, 'Paired')[3], brewer.pal(12, 'Paired')[c(5)]),
         #signal.pch=c(1),
         signal.cex=0.7,
         col=matrix(c(brewer.pal(12, 'Paired')[c(1,3,5)], brewer.pal(12, 'Paired')[c(2,4,6)]),
                    3,2,byrow=F),
         highlight=list(f1=factor_loci[factor_loci$trait=='f1', ]$SNP, f2=factor_loci[factor_loci$trait=='f2', ]$SNP, f3=factor_loci[factor_loci$trait=='f3', ]$SNP),
         highlight.col='black',
         highlight.pch = 18, 
         highlight.cex=1.6,
         file=c("jpg", 'pdf')[i],
         memo="v3",
         dpi=300,
         file.output=T,
         verbose=TRUE, 
         LOG10=T,
         chr.border=F,
         cex.lab=2,
         width=50, 
         height=10)
}
system('mv Mul* /project/aid_sharing/AID_sharing/outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure1/')










