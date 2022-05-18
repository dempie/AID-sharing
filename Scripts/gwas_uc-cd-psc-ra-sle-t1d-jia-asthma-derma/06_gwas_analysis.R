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



#------- overlapp function -----------------------------------------------------
#a function that tests for overlap between two datasets and outputs a vector of overlapping rsIDs. 


test_overlap <-function(input_a, input_b, n_chr){     
  require(data.table)
  
        a <-as.data.table(input_a)
        b <- as.data.table(input_b)
        
        a$CHR <- as.numeric(a$chr)
        a$start <- as.numeric(a$start)
        a$end <- as.numeric(a$end)
        
        b$CHR <- as.numeric(b$chr)
        b$start <- as.numeric(b$start)
        b$end <- as.numeric(b$end)
        
        output <- list()
        overlapping_rsIDs <- list()
        
                for(i in c(1:n_chr)){
                  setkey(b,start, end)
                  output[[i]] <- foverlaps(a[a$CHR==i,], b[b$CHR==i,], type = 'any', which = TRUE, nomatch =NULL)
                  names(output)[[i]] <- paste0('CHR',i)
                  
                  #return the ovelapping SNPs
                  overlapping_rsIDs[[i]] <- input_b[output[[i]]$yid, 'SNP' ]
                  
                }
                
        #return the overlapping loci defined by rsIDs
        invisible(unlist(overlapping_rsIDs))
        
}

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


a<- f1_lb_output$SNP
b <- f2_lb_output$SNP

#--------crohn locus breaker----------------------------------------------------



#------------------------------------------------------------------------------
uc <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/uc_delange-2017.txt', data.table = F)
cd <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/crohn_delange-2017.txt', data.table = F)
psc <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/psc_ji-2016.txt', data.table = F)
jia <- fread('outputs/version3/04_output_sumstats-function/Sumstats_ready_for_munge/jia_lopez-2016_beta_se.txt', data.table = F)
sle <- fread('utputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/sle_beta_bentham-2015.txt', data.table = F)
t1d <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/t1d_chiou-2021_build37.txt', data.table=F) 
asthma <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/asthma_ban-2020.txt', data.table = F)
ra <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/ra_eu_okada-2014.txt', data.table = F)
derma <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/derma_sliz-2021_build37.txt', data.table = F)


# cd_lb_output <- locus.breaker(res = cd, p.label='p', chr.label = )
# uc_lb_output <- locus.breaker()
# ji_lb_oputput <- locus.breaker()
# jia_lb_output <- 
# sle_lb_output <- 
# t1d_lb_output <- 
#   
  
#---- hypercoloc ---------------------------------------------------------------




library(hyprcoloc)
vignette('hyprcoloc')



















