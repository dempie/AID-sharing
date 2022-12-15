library(data.table)
library(ComplexHeatmap)
library(RColorBrewer)
library(qgraph)
library(ChIPpeakAnno)
library(stringr)
#------ load the dataset -------------------------------------------------------
loci_factors <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/factor_loci_moloc_closest_gene_info.txt', data.table = F) #take the original factor loci list


#------ upse tplot of the genomic regions of the plot --------------------------

f1_range <- GRanges(seqnames = loci_factors[loci_factors$trait=='f1',]$chr, IRanges(loci_factors[loci_factors$trait=='f1',]$start, loci_factors[loci_factors$trait=='f1',]$end))
f2_range <- GRanges(seqnames = loci_factors[loci_factors$trait=='f2',]$chr, IRanges(loci_factors[loci_factors$trait=='f2',]$start, loci_factors[loci_factors$trait=='f2',]$end))
f3_range <- GRanges(seqnames = loci_factors[loci_factors$trait=='f3',]$chr, IRanges(loci_factors[loci_factors$trait=='f3',]$start, loci_factors[loci_factors$trait=='f3',]$end))


#find the loci that are physically overlapping among all the three factors 
ovl_f <- findOverlapsOfPeaks(f1_range, f2_range, f3_range , connectedPeaks = 'keepAll')

#a loop to create as many list as the number of traits that have at least one unique locus
#it takes the output from findoverlapps function and the names of the traits as in findoverlapps!!
lists_forupset <- function(ovl_f, traits){ 
          #allocate lists
          require(ComplexHeatmap)
          require(stringr)
          
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


cm <- make_comb_mat(up_n$traits_overlap, mode = 'distinct')
pdf(width = 6, height = 5, file = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure1/upset_plot_loci_ranges_factors.pdf')
UpSet(cm, set_order = c("f1", "f2", "f3"), comb_order = order(comb_size(cm), decreasing = F),
      comb_col = c(brewer.pal(8, 'Set2')[c(7,6,6,6)],brewer.pal(12, 'Paired')[c(2,4,6)]),
      top_annotation = upset_top_annotation(cm, add_numbers = TRUE, height = unit(6, "cm")),
      left_annotation = upset_left_annotation(cm, add_numbers = TRUE, width = unit(3,'cm'), gp = gpar(fill = brewer.pal(5, 'Greys')[1])  ),
      row_title = "Factor", 
      column_title = " Physical intersection of factor loci"
)
dev.off()


#to check if all of this is true makeVennDiagram(ovl_f)




