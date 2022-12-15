locus_lister_all <- function(my_paths, gwas_names) {
  require(data.table)
  require(dplyr)
  require(GenomicRanges)
  loci <- list()
  
  for(i in c(1:length(my_paths))){
    sstats <- fread(my_paths[[i]], data.table = F)
    colnames(sstats) <- toupper(colnames(sstats))
    
    sstats <- select(sstats, c(SNP, CHR,BP, A2, A1, BETA,SE, P))
    #sstats <- sstats[sstats$SNP %in% shared_SNPs,] #select only the SNPs that are shared among all the traits
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

locus_lister_shared <- function(my_paths, gwas_names) {
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

loci_all <- locus_lister_all(mypaths, trait_names)
loci_shared <- locus_lister_shared(mypaths, trait_names)

dim(loci_all)
dim(loci_shared)



data
for(n in c(1:length(trait_names)))
  i=trait_names[[n]]
a <- nrow(loci_all[loci_all$trait==i,  ])
q <- reduce(create_range(loci_all[loci_all$trait==i,  ]))
b <- length(a@ranges@width)
