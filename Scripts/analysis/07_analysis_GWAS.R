#in this scirpt I will perform downstrean analysis of the GWAS results.
library(data.table)






#----locus.breaker function by Nicola, it identifies loci-----------------------
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
#------ run locus breaker function----------------------------------------------

#load the complete summary stats
f1_gwas <- readRDS('outputs/version3/05_GWAS_results/F1_completeGWAS.RDS')
f2_gwas <- readRDS('outputs/version3/05_GWAS_results/F2_completeGWAS.RDS')

nrow(f1_gwas) #4309927
nrow(f2_gwas) #4309927

#remove NAs
f1_gwas_ok <- f1_gwas[which(!is.na(f1_gwas$Pval_Estimate)), ]
f2_gwas_ok <- f2_gwas[which(!is.na(f2_gwas$Pval_Estimate)), ]

nrow(f1_gwas_ok) #4281353
nrow(f2_gwas_ok) #4281353

#run the function

f1_lb_output <- locus.breaker(res= f1_gwas_ok, p.label = 'Pval_Estimate', 
                              chr.label = 'CHR', pos.label = 'BP'  )
f1_lb_output$trait <- rep('f1', nrow(f1_lb_output))

f2_lb_output <- locus.breaker(res= f2_gwas_ok, p.label = 'Pval_Estimate', 
                              chr.label = 'CHR', pos.label = 'BP'  )

f2_lb_output$trait <- rep('f2', nrow(f2_lb_output))

dim(f1_lb_output)
dim(f2_lb_output)


a<- f1_lb_output$SNP
b <- f2_lb_output$SNP

#---------------------------------------

cd_ok <- fread('outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/crohn_delange-2017.txt', data.table = F)
uc_ok <- fread('outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/uc_delange-2017.txt', data.table = F)
psc_ok <- fread('outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/psc_ji-2016.txt', data.table = F)
jia <- fread('outputs/version3/04_output_sumstats-function/Sumstats_ready_for_munge/jia_lopez-2016_beta_se.txt', data.table = F)
sle <- fread('outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/sle_beta_bentham-2015.txt', data.table = F)
t1d <- fread('outputs/version3/01_output_prepare-sumstats/Sumstats_ready_for_munge/t1d_chiou-2021.txt', data.table=F) 

head(cd_ok)

cd_lb_output <- locus.breaker(res = cd_ok, p.label='p', chr.label = 'CHR', pos.label = 'BP')
cd_lb_output$trait <- rep('cd', nrow(cd_lb_output))

uc_lb_output <- locus.breaker(res = uc_ok, p.label='p', chr.label = 'CHR', pos.label = 'BP')
uc_lb_output$trait <- rep('cd', nrow(uc_lb_output))

psc_lb_oputput <- locus.breaker(res = psc_ok, p.label='P', chr.label = 'CHR', pos.label = 'Pos')
psc_lb_oputput$trait <- rep('psc', nrow(psc_lb_oputput)) 

jia_lb_output <- locus.breaker(res = jia, p.label='p', chr.label = 'chromosome', pos.label = 'position')
jia_lb_output$trait<- rep('jia', nrow(jia_lb_output))

sle_lb_output <- locus.breaker(res = sle, p.label='p', chr.label = 'chrom', pos.label = 'pos')
sle_lb_output$trait <- rep('sle', nrow(sle_lb_output))

t1d_lb_output <- locus.breaker(res = t1d, p.label='p', chr.label = 'chr', pos.label = 'position')
t1d_lb_output$trait <- rep('t1d', nrow(t1d_lb_output))
t1d_lb_output<- rename(t1d_lb_output, 'SNP'= rsID)

loci <- list(f1_lb_output, f2_lb_output, cd_lb_output, uc_lb_output, psc_lb_oputput, jia_lb_output, sle_lb_output, t1d_lb_output)

combine <-  function(loci) {
  require(dplyr)
  to_keep <- list(0)
  for( i in c(1:length(loci))) {
    to_keep[[i]] <- select(loci[[i]], c('chr', 'start', 'end', 'SNP', 'trait'))
  }
  
  return(to_keep)
  
}

loci<- combine(loci)
loci <-  do.call(rbind, loci)
saveRDS(loci, 'loci.RDS')

loci=loci[order(loci$start),]
loci=loci[order(loci$chr),]
head(loci)
