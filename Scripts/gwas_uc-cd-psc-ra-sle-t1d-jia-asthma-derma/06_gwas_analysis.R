#in this scirpt I will perform downstrean analysis of the GWAS results.
library(data.table)
library(GenomicRanges)

#----locus.breaker function by Nicola, it identifies loci-----------------------
locus.breaker=function(res,p.sig=5e-8, p.limit=1e-5,hole.size=250000
                       ,p.label="p",chr.label="chr",pos.label="pos"){
  require(data.table)
  
  res = res[order(res[, chr.label], res[,pos.label]),]
  
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


#------- overlapp function -----------------------------------------------------
#a function that tests for overlap between two datasets and outputs a vector of overlapping rsIDs. 







#------ run locus breaker function----------------------------------------------

#load the complete summary stats
my_sumstats <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/gwas_final_withQindex.RDS')

f1_gwas <- my_sumstats$factor1
f2_gwas <- my_sumstats$factor2
f3_gwas <- my_sumstats$factor3

nrow(f1_gwas) #3343158
nrow(f2_gwas) #3343158
nrow(f3_gwas) #3343158

#remove NAs
f1_gwas_ok <- f1_gwas[which(!is.na(f1_gwas$Pval_Estimate)), ]
f2_gwas_ok <- f2_gwas[which(!is.na(f2_gwas$Pval_Estimate)), ]
f3_gwas_ok <- f3_gwas[which(!is.na(f3_gwas$Pval_Estimate)), ]


nrow(f1_gwas_ok) #3294140
nrow(f2_gwas_ok) #3294140
nrow(f3_gwas_ok) #3294140

#run the locus breaker function

f1_lb_output <- locus.breaker(res= f1_gwas_ok, p.label = 'Pval_Estimate', 
                              chr.label = 'CHR', pos.label = 'BP'  )



f2_lb_output <- locus.breaker(res= f2_gwas_ok, p.label = 'Pval_Estimate', 
                              chr.label = 'CHR', pos.label = 'BP'  )


f3_lb_output <- locus.breaker(res= f3_gwas_ok, p.label = 'Pval_Estimate', 
                              chr.label = 'CHR', pos.label = 'BP'  )

dim(f1_lb_output) #70
dim(f2_lb_output) #65
dim(f3_lb_output) #62



#--------crohn locus breaker----------------------------------------------------



#------------------------------------------------------------------------------
uc <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/uc_delange-2017.txt', data.table = F)
cd <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/crohn_delange-2017.txt', data.table = F)
psc <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/psc_ji-2016.txt', data.table = F)
jia <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/jia_beta_lopezisac-2020.txt', data.table = F)
sle <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/sle_beta_bentham-2015.txt', data.table = F)
t1d <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/t1d_chiou-2021_build37.txt', data.table=F) 
colnames(t1d)[1] <- 'SNP'
asthma <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/asthma_ban-2020.txt', data.table = F)
ra <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/ra_eu_okada-2014_chr_bp.txt', data.table = F)
derma <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/derma_sliz-2021_build37.txt', data.table = F)


cd_lb_output <- locus.breaker(res = cd, p.label='p', chr.label ='CHR',  pos.label = 'BP' )

uc_lb_output <- locus.breaker(res = uc, p.label='p', chr.label ='CHR',  pos.label = 'BP' )
psc_lb_oputput <-  locus.breaker(res = psc, p.label='P', chr.label ='CHR',  pos.label = 'BP' )
jia_lb_output <-  locus.breaker(res = jia, p.label='p', chr.label ='CHR',  pos.label = 'BP' )
sle_lb_output <-  locus.breaker(res = sle, p.label='p', chr.label ='CHR',  pos.label = 'BP' )
t1d_lb_output <-  locus.breaker(res =t1d,  p.label='p', chr.label ='CHR',  pos.label = 'BP_37' )
asthma_lb_output <-  locus.breaker(res =asthma,  p.label='p', chr.label ='CHR',  pos.label = 'BP' )
ra_lb_output <-  locus.breaker(res = ra,  p.label='p', chr.label ='CHR',  pos.label = 'BP' )
derma_lb_output <-  locus.breaker(res = derma,  p.label='p', chr.label ='CHR',  pos.label = 'BP' )
  
#-------------------------------------------------------------------------------

#find the overlapping loci!


create_range <- function(res, chr='CHR', start='start', end='end', w=F ) {
  
  obj <-GRanges( seqnames =  res$chr ,IRanges(names = res$chr ,start = as.numeric(res$start), end = as.numeric(res$end))) 
  ifelse(w, return(ranges(obj, use.names = T) ), return(data.frame(names= obj@ranges@NAMES ,start=obj@ranges@start, width=obj@ranges@width))) 
  
}



f1_range <- create_range(f1_lb_output, w = F ) 
f2_range <- create_range(f2_lb_output, w=F)
f3_range <- create_range(f3_lb_output, w=F)
cd_ranges <- create_range()
uc_ranges <- create_range(uc_lb_output, w=F)
#------------------something is wrong with locus breaker 












