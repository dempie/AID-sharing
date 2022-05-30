
library(data.table)
library(dplyr)
library(GenomicRanges)
library(tidyr)
library(hyprcoloc) #https://github.com/jrs95/hyprcoloc

#vignette('hyprcoloc')

#function to create a genomic ranges object useful to find overlapps between loci 

create_range <- function(res, chr='CHR', start='start', end='end', w=T, tag=NA ) {
  obj <-GRanges( seqnames =  res$chr ,IRanges(names = res$chr ,start = as.numeric(res$start), end = as.numeric(res$end))) 
  if( sum(obj@ranges@width> 2000000) >0 ) { warning(paste0( sum(obj@ranges@width> 2000000)) ,' ranges are wider than 2 MB') 
    print(obj[which(obj@ranges@width> 2000000),])
  }
  ifelse(w, return(obj), return(data.frame(names= obj@ranges@NAMES ,start=obj@ranges@start, width=obj@ranges@width))) 
  
}


#--------------------------------------------------------------------------------
sstats_names <- list.files('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/')
sstats_names  <- sstats_names[c(-8, -9, -12)]
mypaths <- as.list(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/', sstats_names))

#loop for loading all the files, create a vector of names and find the SNPs that are present in all the GWAS
trait_names <- list()
SNPs <- list()
list_of_files <- list()

for( i in c(1:length(sstats_names ))){
  #createa lit of trait names
  trait_names[i] <- strsplit(sstats_names , '_')[[i]][[1]]
  
  #load the sumstats and put them into a list
  list_of_files[[i]] <- fread(mypaths[[i]], data.table = F)
  SNPs[[i]]<- list_of_files[[i]]$SNP
}

names(list_of_files) <- unlist(trait_names)

shared_SNPs <- Reduce(intersect, SNPs)

#this function identifies significant loci with the locus.breaker fuction
#then creates a table of all the loci for all traits and puts them togehter.
#then looks at the overlapps between the loci and creates a macro loci that do not overlapp among them and assigns them a unique number (column pan_loci)
#the function requires the sumstatsto have the columns: SNP, CHR,BP, A2 (EFFECT_ALLELE), A1 (NON_EFFECT_ALLELE), BETA,SE, P

locus_lister <- function(my_paths, gwas_names) {
  require(data.table)
  require(dplyr)
  require(GenomicRanges)
  loci <- list()
  
  for(i in c(1:length(my_paths))){
            sstats <- fread(my_paths[[i]], data.table = F)
            colnames(sstats) <- toupper(colnames(sstats))
            
            sstats <- select(sstats, c(SNP, CHR,BP, A2, A1, BETA,SE, P))
            sstats <- sstats[sstats$SNP %in% shared_SNPs,] #select only the SNPs that are shared among all the traits
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
  all_loci[overlapping@from,]$pan_locus <- overlapping@to  #assinging the number as index of which macro loci is overlapping 
  
  return(all_loci)
  
}

#337

#---------generate the list of all loci------------------------------------------

loci <- locus_lister(mypaths, trait_names)
head(loci)
dim(loci) #573 12

# fwrite(loci, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/all_loci_all_traits.txt', sep = '\t', col.names = T, row.names = F, quote = F)
# saveRDS(shared_SNPs, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/shared_SNPs_all_traits.RDS')

  

#start by searching the commong SNPs between all the gwas and generate a referenc set 
#shared_SNPs <- Reduce(intersect, SNPs)  #shared SNPs among all the GWAS
reference_file <- fread('SNP/reference.1000G.maf.0.005.txt.gz', data.table = F) #load the reference file
ref_set <- reference_file[reference_file$SNP   %in%  shared_SNPs, ] #create a reference set of the positions of the shared SNPs
loc_index <- sort(unique(loci$pan_locus))

betas <- list()
SE <- list()
res_coloc <- list()

for(i in c(1:length(unique(loci$pan_locus)))){
  
        #select the active locus, the beginning and the end of the locus
        start_loc <- min(loci[loci$pan_locus==loc_index[i], 'start'] )
        end_loc <- max(loci[loci$pan_locus==loc_index[i], 'end'] )
        chr_loc <- unique(loci[loci$pan_locus==loc_index[i], 'chr'])

        #the active SNP for locus i
        SNP_active <- ref_set[c( between(ref_set$BP,  start_loc ,  end_loc)  & ref_set$CHR== chr_loc ), ]$SNP
        
      
        beta_locus <- data.frame(row.names =SNP_active )
        se_locus <- data.frame(row.names =SNP_active )
        #             
                    
                          #for each GWAS take out the BETAs and the SE
                          for(k in c(1:length(unlist( loci[loci$pan_locus==i,]$trait)))){
                                    #take only the ones that are significant fot that loci
                                    tt<- unlist( loci[loci$pan_locus==i,]$trait)[k] 
                                    tr_SNP_active <- list_of_files[[tt]][ list_of_files[[tt]]$SNP %in% SNP_active ,]
                                    tr_SNP_active <- tr_SNP_active[!duplicated(tr_SNP_active ), ] #remove duplicated SNPs 
                                    #exctract the betas and se
                                    beta_locus[, paste0(tt,'-',k,'_locus', i)] <- select(tr_SNP_active, 'BETA')
                                    se_locus[, paste0(tt,'-',k,'_locus', i)] <- select(tr_SNP_active, 'SE')
                                  
                          }
                    
        
        #run coloc only if the number of traits for that locus is more than one (otherwise you will have an error) and if there are at least 10 SNPs
        if(ncol(beta_locus)>1 & nrow(beta_locus)>10){
          traits <- colnames(beta_locus)
          rsid <- rownames(beta_locus)
          res_coloc[[paste0('locus_', i)]] <- hyprcoloc( effect.est = as.matrix(beta_locus), effect.se = as.matrix(se_locus), trait.names=traits, snp.id=rsid)

        }

        betas[[paste0('locus_', i)]] <- beta_locus
        SE[[paste0('locus_', i)]] <- se_locus
             
        print(i)
     
}      

loci

length(res_coloc)

saveRDS(res_coloc, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/hyprcoloc_outputs.RDS')



  
  