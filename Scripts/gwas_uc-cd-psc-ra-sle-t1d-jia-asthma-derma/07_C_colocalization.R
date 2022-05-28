
library(data.table)
library(dplyr)
library(GenomicRanges)
library(tidyr)
library(hyprcoloc) #https://github.com/jrs95/hyprcoloc

#function to create a genomic ranges object useful to find overlapps between loci 

create_range <- function(res, chr='CHR', start='start', end='end', w=T, tag=NA ) {
  obj <-GRanges( seqnames =  res$chr ,IRanges(names = res$chr ,start = as.numeric(res$start), end = as.numeric(res$end))) 
  if( sum(obj@ranges@width> 2000000) >0 ) { warning(paste0( sum(obj@ranges@width> 2000000)) ,' ranges are wider than 2 MB') 
    print(obj[which(obj@ranges@width> 2000000),])
  }
  ifelse(w, return(obj), return(data.frame(names= obj@ranges@NAMES ,start=obj@ranges@start, width=obj@ranges@width))) 
  
}
#-------------------------------------------------------------------------------


sstats_names <- list.files('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/')
sstats_names  <- sstats_names [c(-8, -9)]
mypaths <- as.list(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/', sstats_names))

trait_names <- list()
for( i in c(1:length(sstats_names ))){
  trait_names[i] <- strsplit(sstats_names , '_')[[i]][[1]]
  
}

#--------------------------------------------------------------------------------

#this function identifies significant loci with the locus.breaker fuction
#then creates a table of all the loci for all traits and puts them togehter.
#then looks at the overlapps between the loci and creates a macro loci that do not overlapp among them and assigns them a unique number (column pan_loci)
#the function requires the sumstatsto have the columns: SNP, CHR,BP,EFFECT_ALLELE, NON_EFFECT_ALLELE, BETA,SE, P

locus_lister <- function(my_paths, gwas_names) {
  require(data.table)
  require(dplyr)
  require(GenomicRanges)
  loci <- list()
  for(i in c(1:length(my_paths))){
    sstats <- fread(my_paths[[i]], data.table = F)
    colnames(sstats) <- toupper(colnames(sstats))
    
    sstats <- select(sstats, c(SNP, CHR,BP, A2, A1, BETA,SE, P))
    loci[[i]] <-  locus.breaker(sstats)  #locus breaker function is in my R profile
    names(loci[i]) <- gwas_names[[i]]
    loci[[i]]$trait <- rep(gwas_names[i], nrow(loci[[i]]))
  }
  
  all_loci <- do.call(rbind, loci) #create the listof all the loci
  
  pan_loci <- reduce(create_range(all_loci)) #create the genomic ranges object with my function create_range
  pan_loci_non_reduced <- create_range(all_loci)
  
  overlapping <- findOverlaps(pan_loci_non_reduced, pan_loci) #find overlaps between all the loci and use it as an index of the unique non overlapping loci
  overlapping
  all_loci$pan_locus <- rep(0, nrow(all_loci)) #allocate the column 
  all_loci[overlapping@from,]$pan_locus <- overlapping@to  #assing the number as index of which macro loci is overlapping 
  
  return(all_loci)
  
}


#---------generate the list of all loci------------------------------------------

sstats_names <- list.files('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/')
sstats_names  <- sstats_names[1:5]
mypaths <- as.list(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/', sstats_names))

trait_names <- list()
for( i in c(1:length(sstats_names ))){
  trait_names[i] <- strsplit(sstats_names , '_')[[i]][[1]]
  
}



loci <- locus_lister(mypaths, trait_names)


trait_1 <- fread(mypaths[[1]], data.table = F)
trait_2 <- fread(mypaths[[2]], data.table = F)
trait_3 <- fread(mypaths[[3]], data.table = F)
traits <- list(trait_1, trait_2, trait_3)  

  SNPs <- list()
  ss <- vector()
for(i in c(1:3)){
  ss <- traits[[i]]
  SNPs[[i]]<- ss$SNP
}


i=16
k=1

#start by searching the commong SNPs between all the gwas and generate a referenc set 
shared_SNPs <- Reduce(intersect, SNPs)  #shared SNPs among all the GWAS
reference_file <- fread('SNP/reference.1000G.maf.0.005.txt.gz', data.table = F) #load the reference file
ref_set <- reference_file[ref_set$SNP   %in%  shared_SNPs, ] #create a reference set of the positions of the shared SNPs
loc_index <- sort(unique(loci$pan_locus))


        #select the active locus, the beginning and the end of the locus
        start_loc <- min(loci[loci$pan_locus==loc_index[i], 'start'] )
        end_loc <- max(loci[loci$pan_locus==loc_index[i], 'end'] )
        chr_loc <- unique(loci[loci$pan_locus==loc_index[i], 'chr'])

        #the active SNP for locus i
        SNP_active <- ref_set[c( between(ref_set$BP,  start_loc ,  end_loc)  & ref_set$CHR== chr_loc ), ]$SNP
        
        betas <- list()
        se <- list()
                #for each GWAS take out the BETAs and the SE
                
        
                tr_SNP_active <- traits[[k]][traits[[k]]$SNP %in% SNP_active ,]
                betas[[k]]  <- select(tr_SNP_active, 'BETA')
                se[[k]] <- select(tr_SNP_active, 'SE')
                
                betas[[locus]]
                SE[[locus]]
                app
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        



