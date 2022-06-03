
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

#this function identifies significant loci with the locus.breaker fuction
#then creates a table of all the loci for all traits and puts them togehter.
#then looks at the overlapps between the loci and creates a macro loci that do not overlapp among them and assigns them a unique number (column pan_loci)
#the function requires the sumstatsto have the columns: SNP, CHR,BP, A2 (EFFECT_ALLELE), A1 (NON_EFFECT_ALLELE), BETA,SE, P
#restar R to clear memory before and after running locus_lister, as it sometimes too large objects are generated
locus_lister <- function(my_paths, gwas_names) {
  # require(data.table)
  # require(dplyr)
  # require(GenomicRanges)
  
  SNPs <- list()
  list_of_files <- list()

  for( i in c(1:length(my_paths))){
    #load the sumstats and put them into a list
    list_of_files[[i]] <- fread(my_paths[[i]], data.table = F)
    SNPs[[i]]<- list_of_files[[i]]$SNP
  }
  
  shared_SNPs <- Reduce(intersect, SNPs) #compute the SNPs that are present in all the summary stats
  names(list_of_files) <- unlist(gwas_names)
  loci <- list()
  
  for(i in c(1:length(my_paths))){
            to_take <- gwas_names[[i]]
            sstats <- list_of_files[[to_take]]
            colnames(sstats) <- toupper(colnames(sstats))
            
            sstats <- select(sstats, c(SNP, CHR,BP, A2, A1, BETA,SE, P))
            sstats <- sstats[sstats$SNP %in% shared_SNPs,] #select only the SNPs that are shared among all the traits
            loci[[i]] <-  locus.breaker(sstats)  #locus breaker function is in my R profile
            loci[[i]]$trait <- rep(to_take, nrow(loci[[i]]))
  }
  
  all_loci <- do.call(rbind, loci) #create the listof all the loci
  
  pan_loci <- reduce(create_range(all_loci)) #create the genomic ranges object with my function create_range
  pan_loci_non_reduced <- create_range(all_loci)
  
  overlapping <- findOverlaps(pan_loci_non_reduced, pan_loci) #find overlaps between all the loci and use it as an index of the unique non overlapping loci
  all_loci$pan_locus <- rep(0, nrow(all_loci)) #allocate the column 
  all_loci[overlapping@from,]$pan_locus <- overlapping@to  #assinging the number as index of which macro loci is overlapping 
  
  return(all_loci)
  
}


#-------------function for colocalization -------------------------------------
#this function takes as input a list of paths of munged summmary stats, thier respective names in list format, the output from locus_breaker. 
#it returns an object with the betas and se involved in each of the loci (the loci argument), the outputs from colocalization and how many snps are involved for each locus. 
#restar R to clear memory before and after running colocalize, as it sometimes too large objects are generated
colocalize <- function(list_of_paths, trait_names, loci){
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

#------------------colocazlize 2------------------------------------------------

#Runs hyprcoloc on all the loci withouth checking if they are physically overlapping
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

#----------- run the function---------------------------------------------------

#gather te info about the paths and the gwas names
sstats_names <- list.files('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/')
sstats_names  <- sstats_names[c(-8, -9, -12)]
the_paths<- as.list(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/', sstats_names))

#loop for loading all the files, create a vector of names and find the SNPs that are present in all the GWAS
gwas_names <- list()


for( i in c(1:length(sstats_names ))){
  #createa lit of trait names
  gwas_names[i] <- strsplit(sstats_names , '_')[[i]][[1]]
}


#--------- all against all -----------------------------------------------------

all_coli <- locus_lister(the_paths, gwas_names)
fwrite(all_coli, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/hyprcoloc_all_against_all/loci_all_traits.txt')
head(all_coli)
dim(all_coli) #573  12

all_coli<- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/hyprcoloc_all_against_all/loci_all_traits.txt',data.table = F)
#restar R to clear memory before and after running, as it sometimes too large objects are generated


all_coloc <- colocalize(list_of_paths = the_paths, trait_names =  gwas_names ,loci = all_coli)
saveRDS(all_coloc, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/hyprcoloc_all_against_all/colocalize_all_loci_all_traits.RDS')

#------only factors ------------------------------------------------------------

factor_loci <- locus_lister(the_paths[c(4,5,6)], gwas_names[c(4,5,6)])
fwrite(factor_loci , 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/hyprcoloc_only_factors/factor_loci.txt', sep = '\t', col.names = T, row.names = F, quote = F)
#factor_loci <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/hyprcoloc_only_factors/factor_loci.txt', data.table = F)

factor_coloc <- colocalize(list_of_paths = the_paths[c(4,5,6)], trait_names = gwas_names[c(4,5,6)],loci = factor_loci)
saveRDS(factor_coloc, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/hyprcoloc_only_factors/colocalize_factors.RDS')









































  