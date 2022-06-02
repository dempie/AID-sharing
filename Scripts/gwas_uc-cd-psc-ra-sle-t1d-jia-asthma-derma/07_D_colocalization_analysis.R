#in this script plotting and analysis of colocalization will be done

library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(tidyr)
library(ChIPpeakAnno)
library(hyprcoloc)
library(stringr)
library(RColorBrewer)



#function to create a genomic ranges object useful to find overlapps between loci 

create_range <- function(res, chr='CHR', start='start', end='end' ) {
  cnd<- length(GRanges( seqnames =  res$chr ,IRanges(names = rownames(res) ,start = as.numeric(res$start), end = as.numeric(res$end)))@ranges) != (length(reduce(GRanges( seqnames =  res$chr ,IRanges(names = rownames(res) ,start = as.numeric(res$start), end = as.numeric(res$end))))))
  if(cnd==T) {warning('There are overlaps inside this thing') }

  return(GRanges( seqnames =  res$chr ,IRanges(names = rownames(res) ,start = as.numeric(res$start), end = as.numeric(res$end))) )
  #check if there are overlapping ranges
}

#----------analysis of the factor gwas------------------------------------------

coloc_factors <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/hyprcoloc_only_factors/colocalize_factors.RDS')
loci_factors <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/hyprcoloc_only_factors/factor_loci.txt', data.table=F)
f1_range <- create_range(loci_factors[loci_factors$trait=='f1',])
f2_range <- create_range(loci_factors[loci_factors$trait=='f2',])
f3_range <- create_range(loci_factors[loci_factors$trait=='f3',])

dim(loci_factors[loci_factors$trait=='f3',])
#find the loci that are physically overlapping among all the three factors 
ovl_f <- findOverlapsOfPeaks(f1_range, f2_range, f3_range , connectedPeaks = 'keepAll')

#a loop to create as many list as the number of traits that have at least one unique locus
#it takes the output from findoverlapps function and the names of the traits as in findoverlapps!!
lists_forupset <- function(ovl_f, traits){ 
    #allocate lists
    require(ComplexHeatmap)
    un <- list()
    final <- list(list())
    sh <- list()
    #create the lists and put it the eleemtnts that are unique
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
                un[[tt]][[length(un[[tt]])+1]]<- paste0(ovl_f$mergedPeaks$peakNames[[k]], collapse = '-')
              }
              
            }
    }
    
    
    #-generate a matrix showing which locus is present in each disease
    output <- list()
    for(k in 1:length(un)){
      tt <- c(names(un))[k]
      to_get <- strsplit(unlist(strsplit(un[[tt]], split = '-')), split='__')
      
      for(i in c(1:length(to_get))){
        output[[tt]][[i]]<-to_get[[i]][2]
        
      }
      output[[tt]]<- unlist(output[[tt]])
      
      
    }
    output <- list_to_matrix(output)
    output <- output[order(as.numeric(rownames(output))),]
    
    #generate the output
    return(list(loci=output, traits_overlap=un))
}



up_n <- lists_forupset(ovl_f = ovl_f, traits = c('f1', 'f2', 'f3')) 


cm <- make_comb_mat(up_n, mode = 'distinct')
pdf(width = 10, height = 5, file = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/upset_plot_factors.pdf')
UpSet(cm, set_order = c("f1", "f2", "f3"), comb_order = order(comb_size(cm), decreasing = T),
      comb_col = c((brewer.pal(6, 'Set3'))[4:6][comb_degree(cm)]),
       top_annotation = upset_top_annotation(cm, add_numbers = TRUE, height = unit(6, "cm")),
      right_annotation = upset_right_annotation(cm, add_numbers = TRUE, width = unit(5,'cm') )
       )
dev.off()


#----- function locus namer-----------------------------------------------------






