##------------

locus.breaker=function(res,p.sig=5e-8, p.limit=1e-5,hole.size=250000
                       ,p.label="p",chr.label="chr",pos.label="pos"){
  
  res=res[which(res[,p.label]<p.limit),]
  trait.res=c()
  for(j in 1:22){
    
    res.chr=res[which(res[,chr.label]==j),]
    if(nrow(res.chr)>1){
      holes=res.chr[,pos.label][-1]-res.chr[,pos.label][-length(res.chr[,pos.label])] 
      gaps=which(holes>hole.size)
      if(length(gaps)>0){
        for(k in 1:(length(gaps)+1)){
          
          if(k==1){
            res.loc=res.chr[1:(gaps[k]),]  
          }else if(k==(length(gaps)+1)){
            res.loc=res.chr[(gaps[k-1]+1):nrow(res.chr),]  
          }else{
            res.loc=res.chr[(gaps[k-1]+1):(gaps[k]),]
          }
          if(min(res.loc[,p.label])<p.sig){
            
            start.pos=min(res.loc[,pos.label],na.rm=T)
            end.pos=max(res.loc[,pos.label],na.rm=T)
            chr=j
            best.snp=res.loc[which.min(res.loc[,p.label]),]
            line.res=c(chr,start.pos,end.pos,unlist(best.snp))
            trait.res=rbind(trait.res,line.res)
          }
          
          
        }
      }else{
        res.loc=res.chr
        if(min(res.loc[,p.label])<p.sig)  {
          
          start.pos=min(res.loc[,pos.label],na.rm=T)
          end.pos=max(res.loc[,pos.label],na.rm=T)
          chr=j
          best.snp=res.loc[which.min(res.loc[,p.label]),]
          line.res=c(chr,start.pos,end.pos,unlist(best.snp))
          trait.res=rbind(trait.res,line.res)
        }
        
      }
      
    }else if(nrow(res.chr)==1){
      
      res.loc=res.chr
      if(min(res.loc[,p.label])<p.sig){
        start.pos=min(res.loc[,pos.label],na.rm=T)
        end.pos=max(res.loc[,pos.label],na.rm=T)
        chr=j
        best.snp=res.loc[which.min(res.loc[,p.label]),]
        line.res=c(chr,start.pos,end.pos,unlist(best.snp))
        trait.res=rbind(trait.res,line.res)
      }
      
      
    }
  }
  
  print(trait.res)
  trait.res=as.data.frame(trait.res,stringsAsFactors=FALSE)
  trait.res=trait.res[,-(which(names(trait.res)==chr.label))]
  names(trait.res)[1:3]=c("chr","start","end")
  trait.res
}

#---------------
#test hyperloc 

library(data.table)
library(qqman)
library(dplyr)



prepare_manhattan <- function(x, p_col, CHR) {
  #x <- x[!psc_plot$CHR=='X',]
  x$CHR <- as.numeric(x$CHR) 
  x <- x[(which(!is.na(select(x, all_of(p_col))))),]
  x <- x[(which(!is.na(select(x, all_of(CHR))))),]
  x<- x[select(x, all_of(p_col))<0.005,]
  
  
  return(x)
}


my_sumstats <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/gwas_final_withQindex.RDS')
f1_gwas <- my_sumstats$factor1


uc <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/uc_delange-2017.txt', data.table = F)
cd <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/crohn_delange-2017.txt', data.table = F)
psc <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/psc_ji-2016.txt', data.table = F)

f1_plot <- prepare_manhattan(f1_gwas, 'Pval_Estimate', 'CHR')
uc_plot <- prepare_manhattan(uc, 'p', 'CHR')
cd_plot <- prepare_manhattan(cd, 'p', 'CHR')
psc_plot <- prepare_manhattan(psc, 'P', 'CHR')

jpeg(file='outputs/tests/miami_f1_uc_psc_cd.jpeg', width = 800, height = 400, units='mm', res = 600)
par(mfrow=c(4,1))
par(mar=c(5,5,3,3))
manhattan(f1_plot, chr="CHR", bp="BP", snp="SNP", p="Pval_Estimate" ,ylim=c(0,30), main='F1', 
          col = c("darksalmon", "darkseagreen4"))  
par(mar=c(5,5,3,3))
manhattan(uc_plot, chr="CHR", bp="BP", snp="SNP", p="p" ,ylim=c(0,30), main='uc', 
          col = c("darksalmon", "darkseagreen4"))
par(mar=c(5,5,3,3))
manhattan(cd_plot, chr="CHR", bp="BP", snp="SNP", p="p" ,ylim=c(0,30), main='cd', 
          col = c("darksalmon", "darkseagreen4"))
manhattan(psc_plot, chr="CHR", bp="BP", snp="SNP", p="P" ,ylim=c(0,30), main='psc', 
          col = c("darksalmon", "darkseagreen4"))
dev.off()



#locus zoom chromosome 1

jpeg(file='outputs/tests/miami_f1_uc_psc_cd_CHR1.jpeg', width = 800, height = 400, units='mm', res = 600)
par(mfrow=c(4,1))
par(mar=c(5,5,3,3))
manhattan(subset(f1_plot, CHR==1), chr="CHR", bp="BP", snp="SNP", p="Pval_Estimate" ,ylim=c(0,30), main='F1', 
          col = c("darksalmon", "darkseagreen4"))  
par(mar=c(5,5,3,3))
manhattan(subset(uc_plot, CHR==1), chr="CHR", bp="BP", snp="SNP", p="p" ,ylim=c(0,30), main='uc', 
          col = c("darksalmon", "darkseagreen4"))
par(mar=c(5,5,3,3))
manhattan(subset(cd_plot, CHR==1), chr="CHR", bp="BP", snp="SNP", p="p" ,ylim=c(0,30), main='cd', 
          col = c("darksalmon", "darkseagreen4"))
manhattan(subset(psc_plot, CHR==1), chr="CHR", bp="BP", snp="SNP", p="P" ,ylim=c(0,30), main='psc', 
          col = c("darksalmon", "darkseagreen4"))
dev.off()


#-------------------------------------------------------------------------------

f1_lb_chr1 <- locus.breaker(res = subset(f1_gwas, CHR==1), p.label='Pval_Estimate', chr.label ='CHR',  pos.label = 'BP' )
cd_lb_chr1 <- locus.breaker(res = subset(cd, CHR==1), p.label='p', chr.label ='CHR',  pos.label = 'BP' )
 cd_lb_chr1[9, 'SNP'] <- "1-25888000-25890373"
uc_lb_chr1 <- locus.breaker(res = uc, p.label='p', chr.label ='CHR',  pos.label = 'BP' )
psc_lb_chr1 <-  locus.breaker(res = subset(psc, CHR==1), p.label='P', chr.label ='CHR',  pos.label = 'BP' )

library(GenomicRanges)

f1_gr <- GRanges(seqnames = c(cd_lb_chr1$chr) ,ranges = IRanges(start = c(as.numeric(f1_lb_chr1$start)), end = c(as.numeric(f1_lb_chr1$end))) )
cd_gr <-  GRanges(seqnames = c(cd_lb_chr1$chr), ranges = IRanges(start = c(as.numeric(cd_lb_chr1$start)), end = c(as.numeric(cd_lb_chr1$end))) )
uc_gr <-  GRanges(seqnames = c(uc_lb_chr1$chr), ranges = IRanges(start = c(as.numeric(uc_lb_chr1$start)), end = c(as.numeric(uc_lb_chr1$end))) )
psc_gr <-  GRanges(seqnames = c(psc_lb_chr1$chr), ranges = IRanges(start = c(as.numeric(psc_lb_chr1$start)), end = c(as.numeric(psc_lb_chr1$end))) )


list_f1 <- GRangesList('f1'=f1_gr, 'cd'=cd_gr, 'uc'=uc_gr,'psc'=psc_gr) 
reduce(list_f1, drop.empty.ranges = F)

a<- reduce(c(f1_gr, cd_gr, psc_gr))




qui <- c(6.6e+7,7.0e+7)
jpeg(file='outputs/tests/miami_f1_uc_psc_cd_CHR1_loucs.jpeg', width = 400, height = 800, units='mm', res = 600)
par(mfrow=c(4,1))
par(mar=c(5,5,3,3))
manhattan(subset(f1_plot, CHR==1), chr="CHR", bp="BP", snp="SNP", p="Pval_Estimate" ,ylim=c(0,30),xlim=c(qui), main='F1', 
          col = c('darkseagreen4'))  
par(mar=c(5,5,3,3))
manhattan(subset(uc_plot, CHR==1), chr="CHR", bp="BP", snp="SNP", p="p" ,ylim=c(0,30),xlim=c(qui), main='uc', 
          col = c('darkseagreen4'))
par(mar=c(5,5,3,3))
manhattan(subset(cd_plot, CHR==1), chr="CHR", bp="BP", snp="SNP", p="p" ,ylim=c(0,30),xlim=c(qui), main='cd', 
          col = c( 'darkseagreen4'))
manhattan(subset(psc_plot, CHR==1), chr="CHR", bp="BP", snp="SNP", p="P" ,ylim=c(0,30),xlim=c(qui), main='psc', 
          col = c('darkseagreen4'))
dev.off()








