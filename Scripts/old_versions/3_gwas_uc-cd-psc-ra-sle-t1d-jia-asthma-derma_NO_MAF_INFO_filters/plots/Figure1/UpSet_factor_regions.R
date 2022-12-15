library(data.table)
library(ComplexHeatmap)
library(RColorBrewer)
library(ChIPpeakAnno)
library(stringr)
#------ load the dataset -------------------------------------------------------
regions.factors <- fread('outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/06_genomic_regions_all_traits_munged/factor.regions.txt', data.table = F) #take the original factor loci list


#------ upse tplot of the genomic regions of the plot --------------------------
f_ranges <- list()
for(i in c('f1','f2', 'f3')){
f_ranges[[i]] <- GRanges(seqnames = regions.factors[regions.factors$trait==i,]$chr, IRanges(regions.factors[regions.factors$trait==i,]$start, regions.factors[regions.factors$trait==i,]$end))
}

#find the loci that are physically overlapping among all the three factors 
ovl_f <- findOverlapsOfPeaks(f_ranges , connectedPeaks = 'keepAll')

#a loop to create as many list as the number of traits that have at least one unique locus
#function for obtaining a combination matrix for plotting and UpSet plot from Granges names, only works in this specific case. 

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



up_n <- lists_forupset( f_ranges, traits = c('f1', 'f2', 'f3')) 


pdf(width = 6, height = 5, file = 'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/plots/figure1/upset_plot_loci_ranges_factors.pdf')
UpSet(up_n, set_order = c("f1", "f2", "f3"), comb_order = order(comb_size(up_n), decreasing = F),
      comb_col = c(brewer.pal(8, 'Set2')[c(7,6,6,6)],brewer.pal(12, 'Paired')[c(2,4,6)]),
      top_annotation = upset_top_annotation(up_n, add_numbers = TRUE, height = unit(6, "cm")),
      left_annotation = upset_left_annotation(up_n, add_numbers = TRUE, width = unit(3,'cm'), gp = gpar(fill = brewer.pal(5, 'Greys')[1])  ),
      row_title = "Factor", 
      column_title = " Physical intersection of factor loci"
)
dev.off()


#to check if all of this is true 
pdf(width = 6, height = 5, file = 'outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/plots/figure1/venn_plot_loci_ranges_factors.pdf')
makeVennDiagram(ovl_f)
dev.off()
system('mv *log /project/aid_sharing/AID_sharing/outputs/3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO_MAF_INFO/plots/figure1/logs')




