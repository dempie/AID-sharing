library(data.table)
library(RColorBrewer)
library(CMplot)


#Q index 
gw_q <- list()
for(i in c(1:3)){
  tt <- c('f1', 'f2', 'f3')[i]
  gw_q[[tt]] <- fread(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/q_index_munged/',tt,'_munged_q_index_build37.txt' ),
                      data.table = F)
}

locus_q_index <- locus.breaker(gw_q$f1, p.label = 'Q_CHISQ_PVAL') 

for(i in 1:3){
        tt <- c('f1', 'f2', 'f3')[i]
        gw_q[[tt]]$P <- -log10(gw_q[[tt]]$P)
        gw_q[[tt]]$Q_CHISQ_PVAL <- log10(gw_q[[tt]]$Q_CHISQ_PVAL)
        a <- gw_q[[tt]][ ,c(1,2,3,11)]
        b <- gw_q[[tt]][ ,c(1,2,3,12)]
        
        
        jpeg(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure2/miami/',tt,'miami_Q_snp.jpg'), width=6000  , height= 3000, res = 300)
    
        par(mar=c(1, 0.5, 0.2, 0.2), mfrow=c(2,1),
            oma = c(4, 4, 0.2, 0.2))
        
        CMplot(a, 
               plot.type="m",
               multracks=F,
               threshold=c(-log10(5e-8)),
               ylim = c(0,30),
               threshold.lty=c(1), 
               band = 1,
               threshold.lwd=c(1), 
               threshold.col=c("black"), 
               amplify=F,
               cex = 0.5,
               #signal.col=list(brewer.pal(12, 'Paired')[c(1)], brewer.pal(12, 'Paired')[3], brewer.pal(12, 'Paired')[c(5)]),
               #signal.pch=c(1),
               signal.cex=0.5,
               col=matrix(c(brewer.pal(12, 'Paired')[c(1,3,5)], brewer.pal(12, 'Paired')[c(2,4,6)]),
                          3,2,byrow=F)[i,],
               highlight= locus_q_index$SNP ,
               highlight.col='black',
               highlight.pch = 13, 
               highlight.cex=1.6,
               dpi=300,
               file.output=F,
               verbose=TRUE, 
               LOG10=F,
               chr.border=F,
               cex.lab=2)
        
        
        
        CMplot(b, 
               plot.type="m",
               multracks=F,
               threshold=c(log10(5e-8)),
               ylim = c(-30,0),
               threshold.lty=c(1), 
               band = 1,
               threshold.lwd=c(1), 
               threshold.col=c("black"), 
               amplify=F,
               cex = 0.5,
               #signal.col=list(brewer.pal(12, 'Paired')[c(1)], brewer.pal(12, 'Paired')[3], brewer.pal(12, 'Paired')[c(5)]),
               #signal.pch=c(1),
               signal.cex=0.5,
               col=brewer.pal(9, 'Set2')[c(7,8)],
               highlight=locus_q_index$SNP,
               highlight.col=  'black',
               highlight.pch = 18, 
               highlight.cex=1.5,
               dpi=300,
               ylab = '-log10(Qsnp Pvalue)',
               file.output=F,
               verbose=TRUE, 
               LOG10=F,
               chr.border=F,
               cex.lab=2)
        
        dev.off()

}
