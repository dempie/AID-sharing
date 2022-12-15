#in this scirpt I will perform downstrean analysis of the GWAS results.
library(data.table)
library(GenomicRanges)
library(dplyr)
library(qqman)
library(ChIPpeakAnno)
library(RColorBrewer)

#----locus.breaker function by Nicola, it identifies loci-----------------------
locus.breaker=function(res,p.sig=5e-8, p.limit=1e-5,hole.size=250000
                       ,p.label="p",chr.label="chr",pos.label="pos"){
  
  res = res[order(as.numeric(res[, chr.label]), as.numeric(res[,pos.label])),]
  
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
prepare_plot <-  function(res, Pval_col,plimit=0.005, BP='BP', CHR='CHR') {
  require(dplyr)
  res[,Pval_col] <- as.numeric(res[, Pval_col])
  res[, CHR] <- as.numeric(res[, CHR])
  res <- res[res[, CHR] %in% c(1:22),]
  output <- res[(which(!is.na(res[,Pval_col]))),]
  output <- output[which(output[, Pval_col] <  plimit),]
  output <- rename(output, BP=all_of(BP), CHR=all_of(CHR), p=all_of(Pval_col) )
  return(output)
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


#-----------------------
#save a copy for running cheers, all the SNPS
for( i in c(1:length(list(f1_gwas_ok, f2_gwas_ok, f3_gwas_ok)))){
  to_save <- list(f1_gwas_ok, f2_gwas_ok, f3_gwas_ok)[[i]]
  to_save <- to_save[,c('SNP', 'CHR', 'BP')]
  to_save <- to_save[base::order(to_save [, 'CHR'], to_save[,'BP']), ]
  to_save$CHR <- paste0('chr', to_save$CHR)
  to_save <- rename(to_save, SNPID ='SNP')
  fwrite(to_save, paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/06_gwas_analysis/for_CHEERS/f',i,'_full_cheers.txt'), 
         col.names = T, row.names = F, sep = '\t', quote = F)
  rm(to_save)
}

#save only the signficant ones  
for( i in c(1:length( list(f1_gwas_ok, f2_gwas_ok, f3_gwas_ok)))){
  to_save <- list(f1_gwas_ok, f2_gwas_ok, f3_gwas_ok)[[i]]
  to_save <- to_save[to_save$Pval_Estimate<5*10^-8,]
  to_save <- to_save[,c('SNP', 'CHR', 'BP')]
  to_save <- to_save[base::order(to_save [, 'CHR'], to_save[,'BP']), ]
  to_save$CHR <- paste0('chr', to_save$CHR)
  to_save <- rename(to_save, SNPID ='SNP')
  fwrite(to_save, paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/06_gwas_analysis/for_CHEERS/f',i,'_all_significant_cheers.txt'), 
         col.names = T, row.names = F, sep = '\t', quote = F)
  rm(to_save)
}


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


f1_lb_output <- locus.breaker(res= f1_gwas_ok, p.label = 'Pval_Estimate', chr.label = 'CHR', pos.label = 'BP'  )
f2_lb_output <- locus.breaker(res= f2_gwas_ok, p.label = 'Pval_Estimate', chr.label = 'CHR', pos.label = 'BP'  )
f3_lb_output <- locus.breaker(res= f3_gwas_ok, p.label = 'Pval_Estimate',  chr.label = 'CHR', pos.label = 'BP'  )
cd_lb_output <- locus.breaker(res = cd, p.label='p', chr.label ='CHR',  pos.label = 'BP' )
uc_lb_output <- locus.breaker(res = uc, p.label='p', chr.label ='CHR',  pos.label = 'BP' )
psc_lb_oputput <-  locus.breaker(res = psc, p.label='P', chr.label ='CHR',  pos.label = 'BP' )
jia_lb_output <-  locus.breaker(res = jia, p.label='p', chr.label ='CHR',  pos.label = 'BP' )
sle_lb_output <-  locus.breaker(res = sle, p.label='p', chr.label ='CHR',  pos.label = 'BP' )
t1d_lb_output <-  locus.breaker(res =t1d,  p.label='p', chr.label ='CHR_37',  pos.label = 'BP_37' )
asthma_lb_output <-  locus.breaker(res =asthma,  p.label='p', chr.label ='CHR',  pos.label = 'BP' )
ra_lb_output <-  locus.breaker(res = ra,  p.label='p', chr.label ='CHR',  pos.label = 'BP' )
derma_lb_output <-  locus.breaker(res = derma,  p.label='p', chr.label ='CHR_37',  pos.label = 'BP_37' )
  

loci_all_gwas <- list(f1=f1_lb_output, f2=f2_lb_output, f3=f3_lb_output, cd=cd_lb_output, uc=uc_lb_output, 
                         psc= psc_lb_oputput, jia=jia_lb_output, sle=sle_lb_output, t1d= t1d_lb_output, asthma= asthma_lb_output, 
                          ra=ra_lb_output,derma= derma_lb_output ) 

saveRDS(loci_all_gwas, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/06_gwas_analysis/loci_all_gwas.RDS' )

loci_all_gwas <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/06_gwas_analysis/loci_all_gwas.RDS')


#save only the index one for running CHEERS
for( i in c(1:length( list(loci_all_gwas$f1, loci_all_gwas$f2, loci_all_gwas$f3)))){
  to_save <- list(loci_all_gwas$f1, loci_all_gwas$f2, loci_all_gwas$f3)[[i]]
  to_save <- to_save[,c('SNP', 'chr', 'BP')]
  to_save$chr <- paste0('chr', to_save$chr)
  to_save <- rename(to_save, SNPID ='SNP')
  fwrite(to_save, paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/06_gwas_analysis/for_CHEERS/f',i,'_index_cheers.txt'), 
         col.names = T, row.names = F, sep = '\t', quote = F)
  rm(to_save)
}



#--------find the overlapp between the loci-------------------------------------

#function to create a genomic ranges object useful to find overlapps between loci 

create_range <- function(res, chr='CHR', start='start', end='end', w=T, tag=NA ) {
  obj <-GRanges( seqnames =  res$chr ,IRanges(names = res$chr ,start = as.numeric(res$start), end = as.numeric(res$end))) 
  if( sum(obj@ranges@width> 2000000) >0 ) { warning(paste0( sum(obj@ranges@width> 2000000)) ,' ranges are wider than 2 MB') 
          print(obj[which(obj@ranges@width> 2000000),])
    }
  ifelse(w, return(obj), return(data.frame(names= obj@ranges@NAMES ,start=obj@ranges@start, width=obj@ranges@width))) 
  
}


#reduce function collapses the overlapping intervals into one single range
f1_range <- reduce(create_range(loci_all_gwas$f1)) 
f2_range <- reduce(create_range(loci_all_gwas$f2))
f3_range <- reduce(create_range(loci_all_gwas$f3))
cd_ranges <- reduce(create_range(loci_all_gwas$cd))
uc_ranges <- reduce(create_range(loci_all_gwas$uc))
psc_ranges <- reduce(create_range(loci_all_gwas$psc))
jia_ranges <- reduce(create_range(loci_all_gwas$jia))
sle_ranges <- reduce(create_range(loci_all_gwas$sle))
t1d_ranges <- reduce(create_range(loci_all_gwas$t1d))
asthma_ranges <- reduce(create_range(loci_all_gwas$asthma))
ra_ranges <- reduce(create_range(loci_all_gwas$ra))
derma_ranges <- reduce(create_range(loci_all_gwas$derma))


#---------------- unique new loci ----------------------------------------------

#search for loci that are present only in the factors and not in the unique traits

#f1
f1_traits <- reduce(c(cd_ranges, uc_ranges, psc_ranges))
f1_not_unique <- findOverlaps(f1_traits, f1_range)
f1_new_unique <- f1_range[-f1_not_unique@to, ]

#which are the most significative SNP?
f1_not_reduce <- create_range(loci_all_gwas$f1)
filter1 <- findOverlaps(f1_not_reduce, f1_new_unique)
loci_all_gwas$f1[filter1@from,][,c('SNP', 'BP', "chr", "start", "end", "est", "Pval_Estimate")] #new SNPs F1 "rs395157"  "rs9367849"


#f2
f2_traits <- reduce(c(jia_ranges, sle_ranges, t1d_ranges, ra_ranges, sle_ranges))
f2_not_unique <- findOverlaps(f2_traits, f2_range)
f2_new_unique <- f2_range[-f2_not_unique@to,]

#which are the most significative SNP in F"?
f2_not_reduce <- create_range(loci_all_gwas$f2)
filter2<- findOverlaps(f2_not_reduce, f2_new_unique)
loci_all_gwas$f2[filter2@from,][,c('SNP', 'BP', "chr", "start", "end", "est", "Pval_Estimate")]  #new SNPs F2 "rs7511656" "rs2366644" "rs4269168" "rs2746432" "rs6582578"


#f3
f3_traits <- reduce(c(derma_ranges, asthma_ranges))
f3_not_unique <- findOverlaps(f3_traits, f3_range) 
f3_new_unique <- f3_range[-f3_not_unique@to,]

#which are the most significative SNP in F3?
f3_not_reduce <- create_range(loci_all_gwas$f3)
filter3 <- findOverlaps(f3_not_reduce, f3_new_unique)
loci_all_gwas$f3[filter3@from,][,c('SNP', 'BP', "chr", "start", "end", "est", "Pval_Estimate")]  #new SNPs F3  "rs11265449"


#-------------------------------------------------------------------------------
#venn diagram of factor loci, to see which overlapp anche which do not 


myCol <- brewer.pal(1:3, "Set2")



pdf(width = 14, height = 14, file = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/06_gwas_analysis/venn_f1_f2_f3_loci.pdf')
venn_f1_f2_f3_loci<- makeVennDiagram(Peaks=list(f1=f1_range, f2=f2_range, f3=f3_range),connectedPeaks = 'keepAll',
                      NameOfPeaks=c("F1", "F2", "F3"),
                      imagetype="png" ,
                      height = 480 , 
                      width = 480 , 
                      resolution = 300,
                      compression = "lzw",
                      lwd = 1,
                      col=myCol,
                      fill = myCol,
                      cex = 1.5,
                      fontfamily = "sans",
                      cat.cex = 2,
                      cat.default.pos = "outer",
                      cat.pos = c(-27, 27, 135),
                      cat.dist = c(0.055, 0.055, 0.085),
                      cat.fontfamily = "sans",
                      cat.col = myCol,
                      rotation = 1
                      )

dev.off()


#-------loci that are shared by all factors------------------------------------
ol <- findOverlapsOfPeaks(f1_range, f2_range, f3_range )

#extrract the info on the 4 loci that are overlapping with everyhting 

ol$peaklist$`f1_range///f2_range///f3_range` 
# F1 SHARED 20, 34, 54, 56
#F2 SHARED 13, 24, 47, 52
#F3 SHARED 16, 25, 48, 50

ol$all.peaks$f1_range

F1_shared <- loci_all_gwas$f1[c(20, 34,54, 56),]$SNP
F2_shared <- loci_all_gwas$f2[c(13, 24, 47, 52),]$SNP
F3_shared <- loci_all_gwas$f3[c(16,25,48,50),]$SNP



?hyprcoloc



dim(loci_all_gwas$f3)







