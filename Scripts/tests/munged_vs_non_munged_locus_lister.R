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
    
    sstats <- select(sstats, c(SNP, CHR,BP, A2, A1, EST ,SE, PVAL_ESTIMATE))
    sstats <- sstats[sstats$SNP %in% shared_SNPs,] #select only the SNPs that are shared among all the traits
    loci[[i]] <-  locus.breaker(sstats, p.label = 'PVAL_ESTIMATE')  #locus breaker function is in my R profile
    loci[[i]]$trait <- rep(to_take, nrow(loci[[i]]))
  }
  
  all_loci <- do.call(rbind, loci) #create the listof all the loci
  
  pan_loci <- reduce(create_range(all_loci)) #create the genomic ranges object with my function create_range
  pan_loci_non_reduced <- create_range(all_loci)
  
  overlapping <- findOverlaps(pan_loci_non_reduced, pan_loci) #find overlaps between all the loci and use it as an index of the unique non overlapping loci
  all_loci$pan_locus <- rep(0, nrow(all_loci)) #allocate the column 
  all_loci[overlapping@from,]$pan_locus <- overlapping@to  #assinging the number as index of which macro loci is overlapping 
  all_loci$pan_locus_name <- rep(0, nrow(all_loci))
  
  #assign a name refereed to the position of each pan_locus
  for(k in 1:length(unique(all_loci$pan_locus))){
    all_loci[all_loci$pan_locus==k, ]$pan_locus_name <- paste0(all_loci[all_loci$pan_locus==k, ]$chr,'_',  min(all_loci[all_loci$pan_locus==k, ]$start), '_',  max(all_loci[all_loci$pan_locus==k, ]$end))
  }
  
  return(all_loci)
  
}


create_range <- function(res, chr='CHR', start='start', end='end', w=T, tag=NA ) {
  obj <-GRanges( seqnames =  res$chr ,IRanges(names = res$chr ,start = as.numeric(res$start), end = as.numeric(res$end))) 
  if( sum(obj@ranges@width> 2000000) >0 ) { warning(paste0( sum(obj@ranges@width> 2000000)) ,' ranges are wider than 2 MB') 
    print(obj[which(obj@ranges@width> 2000000),])
  }
  ifelse(w, return(obj), return(data.frame(names= obj@ranges@NAMES ,start=obj@ranges@start, width=obj@ranges@width))) 
  
}



#non munged
the_paths <- list(f1='outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/factor1_gwas_final_withQindex.txt', 
                    f2='outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/factor2_gwas_final_withQindex.txt',
                    f3='outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/factor3_gwas_final_withQindex.txt')


#non munged
loc_non_munged <- locus_lister(the_paths, list('f1', 'f2', 'f3'))

#munged
factor_loci <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/factor_loci_moloc_closest_gene_info.txt', data.table = F)


table(loc_non_munged$SNP %in% factor_loci$SNP)
table(factor_loci$SNP %in% loc_non_munged$SNP)





