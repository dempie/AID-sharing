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







colocalize_2 <- function(list_of_paths, trait_names, loci){
  #packages 
  # require(data.table)
  # require(dplyr)
  # require(GenomicRanges)
  # require(hyprcoloc)
  
    SNPs <- list()
    list_of_files <- list()
    
    for( i in c(1:length(trait_names))){
      #load the sumstats and put them into a list
      list_of_files[[i]] <- fread(list_of_paths[[i]], data.table = F)
      SNPs[[i]]<- list_of_files[[i]]$SNP
    }
    
    names(list_of_files) <- unlist(trait_names)
    #start by searching the commong SNPs between all the gwas and generate a referenc set 
    shared_SNPs <- Reduce(intersect, SNPs) 
    cat(paste0('The number of shared SNPs is ', length(shared_SNPs )))
    
    reference_file <- fread('SNP/reference.1000G.maf.0.005.txt.gz', data.table = F) #load the reference file
    ref_set <- reference_file[reference_file$SNP   %in%  shared_SNPs, ] #create a reference set of the positions of the shared SNPs
    loc_index <- sort(unique(loci$pan_locus))
    
    betas <- list()
    SE <- list()
    res_coloc <- list()
    report <- data.frame(row.names = paste0('locus_', loc_index )  ) #allocate space to produce a report of the number of SNPs per locus
    
    for(i in c(1:length(unique(loci$pan_locus)))){
      
      #select the active locus, the beginning and the end of the locus
      start_loc <- min(loci[loci$pan_locus==loc_index[i], 'start'] )
      end_loc <- max(loci[loci$pan_locus==loc_index[i], 'end'] )
      chr_loc <- unique(loci[loci$pan_locus==loc_index[i], 'chr'])
      
      #the active SNP for locus i
      SNP_active <- ref_set[c( between(ref_set$BP,  start_loc ,  end_loc)  & ref_set$CHR== chr_loc ), ]$SNP
      
      #if there are less than 10 SNPs raise the window of SNPs to take of +-100 kb 
      if(length(SNP_active)<10) {
        start_loc_2 <-  as.numeric(ifelse((as.numeric(start_loc) - 100000)<0, 0, as.numeric(start_loc) - 100000 )) #not below zero 
        end_loc_2 <- as.numeric(as.numeric(end_loc) + 100000 )
        #select the active SNPs
        SNP_active <- ref_set[c( between(ref_set$BP,  start_loc_2 ,  end_loc_2)  & ref_set$CHR== chr_loc ), ]$SNP
      } else {
        #select the active SNPs
        SNP_active <- ref_set[c( between(ref_set$BP,  start_loc ,  end_loc)  & ref_set$CHR== chr_loc ), ]$SNP
      }
      
      #allocate the space
      beta_locus <- data.frame(row.names =SNP_active )
      se_locus <- data.frame(row.names =SNP_active )
      
      #for each GWAS take out the BETAs and the SE
      for(k in c(1:length(list_of_files))){
        #take all of them
        tt<- unlist(trait_names)[k] 
        tr_SNP_active <- list_of_files[[tt]][ list_of_files[[tt]]$SNP %in% SNP_active ,]
        tr_SNP_active <- tr_SNP_active[!duplicated(tr_SNP_active ), ] #remove duplicated SNPs 
        #exctract the betas and se
        beta_locus[, paste0(tt,'-',k,'_locus', i)] <- select(tr_SNP_active, 'BETA')
        se_locus[, paste0(tt,'-',k,'_locus', i)] <- select(tr_SNP_active, 'SE')
        
      }
      
      #run coloc only  if there are at least 10 SNPs
      if( nrow(beta_locus)>10){
        traits <- colnames(beta_locus)
        rsid <- rownames(beta_locus)
        res_coloc[[paste0('locus_', i)]] <- hyprcoloc( effect.est = as.matrix(beta_locus), effect.se = as.matrix(se_locus), trait.names=traits, snp.id=rsid, binary.outcomes = rep(1, length(traits)))
  
    } 
      
      betas[[paste0('locus_', i)]] <- beta_locus
      SE[[paste0('locus_', i)]] <- se_locus
      report[i,'SNP_active'] <-  length(SNP_active) 
      report[i,'traits'] <- paste0(loci[loci$pan_locus==i,]$trait, collapse = ',')
      print(i)
  }  
    
  
  output <- list(betas=betas, SE=SE, res_coloc=res_coloc, report=report) 
  return(output)
  }








