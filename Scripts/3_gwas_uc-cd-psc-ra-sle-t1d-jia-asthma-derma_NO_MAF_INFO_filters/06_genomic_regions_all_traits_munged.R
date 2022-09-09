
library(data.table)
library(GenomicRanges)
library(dplyr)
library(qqman)
library(ChIPpeakAnno)
library(RColorBrewer)
library(ComplexHeatmap)
library(stringr)

#----locus.breaker function by Nicola, it identifies loci-----------------------


#this function identifies significant genomic regions with the locus.breaker fuction
#then creates a table of all the loci for all traits and puts them togehter.
#then looks at the overlapps between the loci and creates a macro loci that do not overlapp among them and assigns them a unique number (column pan_loci)
#the function requires the sumstatsto have the columns: SNP, CHR,BP, A2 (EFFECT_ALLELE), A1 (NON_EFFECT_ALLELE), BETA,SE, P
#restar R to clear memory before and after running locus_lister, as it sometimes too large objects are generated

locus_lister <- function(my_paths, gwas_names) {
  
      require(data.table)
      require(dplyr)
      require(GenomicRanges)
      
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
      
      pan_loci <- reduce(GRanges(seqnames = all_loci$chr, IRanges(as.numeric(all_loci$start), as.numeric(all_loci$end)))) 
      pan_loci_non_reduced <- GRanges(seqnames = all_loci$chr, IRanges(as.numeric(all_loci$start), as.numeric(all_loci$end)))
      
      overlapping <- findOverlaps(pan_loci_non_reduced, pan_loci) #find overlaps between all the loci and use it as an index of the unique non overlapping loci
      all_loci$pan_locus <- rep(0, nrow(all_loci)) #allocate the column 
      all_loci[overlapping@from,]$pan_locus <- overlapping@to  #assinging the number as index of which macro loci is overlapping 
      all_loci$pan_locus_name <- rep(0, nrow(all_loci))
      
      #assign a name refereed to the position of each pan_locus
      for(k in 1:length(unique(all_loci$pan_locus))){
              all_loci[all_loci$pan_locus==k, ]$pan_locus_name <- paste0(all_loci[all_loci$pan_locus==k, ]$chr,'_',  min(all_loci[all_loci$pan_locus==k, ]$start), '_',  max(all_loci[all_loci$pan_locus==k, ]$end))
      }
      
      rownames(all_loci) <- NULL
      return(all_loci)
  
}


factors <- list('outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/05_munge_for_nicola/munged/f1_munged_build37.txt',
          'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/05_munge_for_nicola/munged/f2_munged_build37.txt',
          'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/05_munge_for_nicola/munged/f3_munged_build37.txt'
          )



factor_regions <- locus_lister(factors, list('f1', 'f2', 'f3')) 
fwrite(factor_regions, 'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/06_genomic_regions_all_traits_munged/factor.regions.txt', sep='\t', col.names = T, row.names = F)


#regions all gwas

files <- list.files('outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/05_munge_for_nicola/munged/') #get the file names
my_paths <- as.list(paste0('outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/05_munge_for_nicola/munged/', files)) #create paths
my_paths <- my_paths[str_ends(my_paths, 'txt')]



gwas_names <- list()
for( i in c(1:length(files))){
  gwas_names[i] <- strsplit(files, '_')[[i]][[1]]
  
}

gwas_names <- gwas_names[c(-11, -8)]
names(my_paths) <- gwas_names

all_regions <- locus_lister(my_paths, gwas_names)
fwrite(all_regions, 'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/06_genomic_regions_all_traits_munged/all_regions.txt', sep='\t', col.names = T, row.names = F)




#------- search for ovelapping genomic regions between factors and their respective traits ---------
regions_factors <- fread('outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/06_genomic_regions_all_traits_munged/factor.regions.txt', data.table = F)
all_regions <- fread('outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/06_genomic_regions_all_traits_munged/all_regions.txt', data.table=F)

f_ranges <- list()
for(i in c('f1', 'f2', 'f3')){
  
  f_ranges[[i]] <- GRanges(seqnames = regions_factors[regions_factors$trait==i,]$chr, IRanges(as.numeric(regions_factors[regions_factors$trait==i,]$start),as.numeric(regions_factors[regions_factors$trait==i,]$end)))
  
}

#plot venndiagram
pdf(width = 14, height = 14, file = 'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/06_genomic_regions_all_traits_munged/venn_f1_f2_f3_loci.pdf')
venn_f1_f2_f3_loci<- makeVennDiagram(Peaks=f_ranges,connectedPeaks = 'keepAll',
                                     NameOfPeaks=c("F1", "F2", "F3"),
                                     imagetype="png" ,
                                     height = 480 , 
                                     width = 480 , 
                                     resolution = 300,
                                     compression = "lzw",
                                     lwd = 1,
                                     cex = 1.5,
                                     cat.cex = 2,
                                     cat.default.pos = "outer",
                                     cat.pos = c(-27, 27, 135),
                                     cat.dist = c(0.055, 0.055, 0.085),
                                     cat.fontfamily = "sans",
                                     rotation = 1
)

dev.off()

#----- plot UpSet plot ---------------------------------------------------------

#function for obtaining a combination matrix for plotting and UpSet plot from Granges names
lists_forupset <- function(peaks, traits){ 
  #allocate lists
  require(ComplexHeatmap)
  require(ChIPpeakAnno)
  require(stringr)
  
  ovl_f <- findOverlapsOfPeaks(peaks)
  
  un <- list()
  final <- list(list())
  sh <- list()
  #create the lists and put it the elements that are unique for each trait
  for( i in 1:length(traits)){
    tt<-traits[i]
    filt <- str_starts(names(ovl_f$uniquePeaks@ranges),tt, negate = FALSE)
    un[[tt]] <- names(ovl_f$uniquePeaks@ranges)[filt]
  }
  
  #put into the lists the element that are shared among the lists
  for(k in 1:length(ovl_f$mergedPeaks$peakNames)){
    
    chek_in <- traits %in% substring(ovl_f$mergedPeaks$peakNames[[k]], first = 1, last = 2) #check if there is f1, f2 or f3
    
    for(u in 1:length(chek_in)){
      tt<-traits[u]
      if(chek_in[u]==T){
        un[[tt]][[length(un[[tt]])+1]] <- paste0(ovl_f$mergedPeaks$peakNames[[k]], collapse = '-')
      }
      
    }
  }
  
  #generate the output
  
  cm <- make_comb_mat(un, mode = 'distinct')
  return(cm)
  
}


#plot

cm <- lists_forupset(f_ranges, traits = c('f1', 'f2', 'f3')) 

pdf(width = 14, height = 14, file = 'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/06_genomic_regions_all_traits_munged/upset_f1_f2_f3_loci.pdf')
UpSet(cm, set_order = c("f1", "f2", "f3"), comb_order = order(comb_size(cm), decreasing = F),
      comb_col = c(brewer.pal(8, 'Set2')[c(7,6,6,6)],brewer.pal(12, 'Paired')[c(2,4,6)]),
      top_annotation = upset_top_annotation(cm, add_numbers = TRUE, height = unit(6, "cm")),
      left_annotation = upset_left_annotation(cm, add_numbers = TRUE, width = unit(3,'cm'), gp = gpar(fill = brewer.pal(5, 'Greys')[1])  ),
      row_title = "Factor", 
      column_title = " Physical intersection of factor loci"
)
dev.off()


#---------------- f1 and its loci -----------------------------------------------------
# 
# 
# f1_ranges <- list()
# for(i in c('f1', 'cd', 'uc', 'psc')){
#   
#   f1_ranges[[i]] <- GRanges(seqnames = all_regions[all_regions$trait==i,]$chr, IRanges(as.numeric(all_regions[all_regions$trait==i,]$start),as.numeric(all_regions[all_regions$trait==i,]$end)))
#   
# }
# 
# 
# cm<- lists_forupset(f1_ranges, c('f1', 'cd', 'uc', 'ps'))
# UpSet(cm)
# 
# 
# 
# UpSet(cm, set_order = c('f1', 'cd', 'uc', 'ps'),
#       top_annotation = upset_top_annotation(cm, add_numbers = TRUE, height = unit(6, "cm")),
#       left_annotation = upset_left_annotation(cm, add_numbers = TRUE, width = unit(3,'cm'), gp = gpar(fill = brewer.pal(5, 'Greys')[1])  ),
#       row_title = "Factor", 
#       column_title = " Physical intersection of F1 with its trait loci"
# )
# #plot venndiagram
# #pdf(width = 14, height = 14, file = 'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/06_genomic_regions_all_traits_munged/venn_f1_uc_cd_psc.pdf')
# venn_f1_f2_f3_loci<- makeVennDiagram(Peaks=f1_ranges,connectedPeaks = 'keepAll',
#                                      NameOfPeaks=c('f1', 'cd', 'uc', 'psc'),
#                                      imagetype="png" ,
#                                  
# )
# 
# 
# 
# #-----f2 and its loci --------------------------------------------------------
# 
# f2_ranges <- list()
# for(i in c('f2', 't1d', 'ra', 'sle', 'jia')){
#   
#   f2_ranges[[i]] <- GRanges(seqnames = all_regions[all_regions$trait==i,]$chr, IRanges(as.numeric(all_regions[all_regions$trait==i,]$start),as.numeric(all_regions[all_regions$trait==i,]$end)))
#   
# }
# 
# 
# cm<- lists_forupset(f2_ranges, c('f2', 't1', 'ra', 'sl', 'ji'))
# UpSet(cm)
# 
# 
# 
# UpSet(cm, set_order =  c('f2', 't1', 'ra', 'sl', 'ji'),
#       top_annotation = upset_top_annotation(cm, add_numbers = TRUE, height = unit(6, "cm")),
#       left_annotation = upset_left_annotation(cm, add_numbers = TRUE, width = unit(3,'cm'), gp = gpar(fill = brewer.pal(5, 'Greys')[1])  ),
#       row_title = "Factor", 
#       column_title = " Physical intersection of F2 with its trait loci"
# )
# #plot venndiagram
# #pdf(width = 14, height = 14, file = 'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/06_genomic_regions_all_traits_munged/venn_f1_uc_cd_psc.pdf')
# venn_f1_f2_f3_loci<- makeVennDiagram(Peaks=f1_ranges,connectedPeaks = 'keepAll',
#                                      NameOfPeaks=c('f1', 'cd', 'uc', 'psc'),
#                                      imagetype="png" ,
#                                      
# )


