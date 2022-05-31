#in this script plotting and analysis of colocalization will be done

library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(tidyr)
library(ChIPpeakAnno)
library(hyprcoloc)


#function to create a genomic ranges object useful to find overlapps between loci 

create_range <- function(res, chr='CHR', start='start', end='end', w=T, tag=NA ) {
  obj <-GRanges( seqnames =  res$chr ,IRanges(names = res$chr ,start = as.numeric(res$start), end = as.numeric(res$end))) 
  if( sum(obj@ranges@width> 2000000) >0 ) { warning(paste0( sum(obj@ranges@width> 2000000)) ,' ranges are wider than 2 MB') 
    print(obj[which(obj@ranges@width> 2000000),])
  }
  ifelse(w, return(obj), return(data.frame(names= obj@ranges@NAMES ,start=obj@ranges@start, width=obj@ranges@width))) 
  
}



#----------analysis of the factor gwas------------------------------------------
coloc_factors <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/hyprcoloc_only_factors/colocalize_factors.RDS')
loci_factors <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/hyprcoloc_only_factors/factor_loci.txt', data.table=F)
f1_range <- reduce(create_range(loci_factors[loci_factors$trait=='f1',])) 
f2_range <- reduce(create_range(loci_factors[loci_factors$trait=='f2',]))
f3_range <- reduce(create_range(loci_factors[loci_factors$trait=='f3',]))


#find the loci that are physically overlapping among all the three factors 
ovl_f <- findOverlapsOfPeaks(f1_range, f2_range, f3_range )

st <- ovl_f$peaklist$`f1_range///f2_range///f3_range`@ranges@start
en <- ovl_f$peaklist$`f1_range///f2_range///f3_range`@ranges@start + ovl_f$peaklist$`f1_range///f2_range///f3_range`@ranges@width -1
chr <-  ovl_f$peaklist$`f1_range///f2_range///f3_range`@seqnames@values
loc_index <- sort(unique(loci_factors$pan_locus))




report <- coloc_factors$report
dim(report[ report$unique==F,])== length(coloc_factors$res_coloc) #TRUE, we were able to run colocalization on all of the loci that showed overlapp 
report$unique <- rep(NA, nrow(report))

#a loop for identifying which are the loci that do not overlap and therefore do not colocalize for sure. 
a <- strsplit(coloc_factors$report$traits, split = ',') # trait elements
strsplit(rownames(coloc_factors$report), split = '_') #row names

for(i in 1: length(a)){
  if(length(a[[i]])==1) {
    report[i,'unique'] <- T  
  } else {
    report[i,'unique'] <- F
  }
  
}

#a loop to create as many list as the number of traits that have at least one unique locus
unique_locus <- rownames(report[ report$unique==T,])
unq <- report[ report$unique==T,]
f <- list(list())

for(i in 1:length(unique(unq$traits))){
  tt <- unique(unq$traits)[i]
  f[tt] <- list(rownames(unq[unq$traits==tt,])). 
}



apply(coloc_factors$report, MARGIN = 1, FUN = function(x)strsplit(x[2],split = ','))




length(coloc_factors$res_coloc)








