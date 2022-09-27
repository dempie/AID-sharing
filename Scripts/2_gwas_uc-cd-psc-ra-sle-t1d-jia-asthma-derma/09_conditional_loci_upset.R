library(data.table)
library(ComplexHeatmap)
library(RColorBrewer)
library(qgraph)

#------ load the dataset -------------------------------------------------------
loci.table <- fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/09_loci_definitions_Nicola/loci_definitions/final_locus_table.tsv', data.table = F) #take the original factor loci list
loci.table$final.locus=paste0(loci.table$Chr,":",loci.table$start,"-",loci.table$end,"_",loci.table$sub_locus)



# for(i in 1:3){
#       tt <- c('f1', 'f2', 'f3')[i]
#       actf <- loci_table[loci_table$trait==tt,]
# 
#      actf[actf$final.locus %in% unique(actf$final.locus[ duplicated(actf$final.locus)]), ]
# }







#----- make an upset plot-------------------------------------------------------
q <- make_comb_mat(list(f1=loci.table[loci.table$trait=='f1',]$final.locus, f2=loci.table[loci.table$trait=='f2',]$final.locus, f3=loci.table[loci.table$trait=='f3',]$final.locus), mode = 'distinct')

pdf(width = 10, height = 4, file = 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/09_loci_definitions_Nicola/upset_plot_final_loci_nicola.pdf')
UpSet(q, set_order = c("f1", "f2", "f3"), 
      comb_order = c(5,6,7,2,3,4,1),
      comb_col = c(brewer.pal(8, 'Set2')[c(7,6,6,6)],brewer.pal(12, 'Paired')[c(2,4,6)]),
      top_annotation = upset_top_annotation(q, add_numbers = TRUE, height = unit(6, "cm")),
      right_annotation = upset_right_annotation(q, add_numbers = TRUE, width = unit(5,'cm'),gp = gpar(fill = brewer.pal(5, 'Greys')[1])  ),
      row_title = "Factor", 
      column_title = "Intersection of all loci, colocalization")
dev.off()


head(loci.table)
coloc_tab <-  fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/09_loci_definitions_Nicola/loci_definitions/final_locus_table.tsv', data.table = F)
coloc_tab$final.locus=paste0(coloc_tab$Chr,":",coloc_tab$start,"-",coloc_tab$end,"_",coloc_tab$sub_locus)

topl_205 <- coloc_tab[coloc_tab$pan.locus==205, c('trait', 'sub_locus', 'SNP', 'final.locus')]



pdf('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/09_loci_definitions_Nicola/qgraph_locus_205_factors.pdf', width = 16, height = 16)
qgraph(topl_205[!topl_205$sub_locus %in% c(7, 10, 9, 8), c(1,4)], 
       directed=F, 
       layout='groups',
       groups=list(c(1,3,7,8,9), c(4,6,11,12), c(2,5,10)),
       border.color=c("#1F78B4" , "#E31A1C", "#1F78B4" ,"#33A02C", "#E31A1C" ,"#33A02C" ,
                      "#A6CEE3" ,"#A6CEE3", '#A6CEE3' ,"#FB9A99" ,"#B2DF8A", "#B2DF8A" ),
       # label.color= c("#1F78B4" , "#E31A1C", "#1F78B4" ,"#33A02C", "#E31A1C" ,"#33A02C" ,"#1F78B4" ,"#1F78B4", "#1F78B4" ,"#E31A1C" ,"#33A02C", "#33A02C" ),
       border.width= rep(2, 12),
       color=c("#1F78B4" , "#E31A1C", "#1F78B4" ,"#33A02C", "#E31A1C" ,"#33A02C" , "#A6CEE3" ,"#A6CEE3", '#A6CEE3' ,"#FB9A99" ,"#B2DF8A", "#B2DF8A" ),
       #edge.labels=topl_205[!topl_205$sub_locus %in% c(7, 10, 9, 8), 3], 
       edge.label.cex=0.5)
dev.off()





##### naming the nodes as the SNP that has the lower pvalue
topl_205 <- coloc_tab[coloc_tab$pan.locus==205, ]
topl_205$to_name <- rep(0, nrow(topl_205))


for(i in 1: length(unique(topl_205[!topl_205$sub_locus %in% c(7, 10, 9, 8),]$final.locus))){
  ll <- unique(topl_205[!topl_205$sub_locus %in% c(7, 10, 9, 8),]$final.locus)[i]
  act <- topl_205[!topl_205$sub_locus %in% c(7, 10, 9, 8) & topl_205$final.locus==ll, ]
  topl_205[!topl_205$sub_locus %in% c(7, 10, 9, 8) & topl_205$final.locus==ll, ]$to_name <- rep(act[order(act$pJ),]$SNP[1], length(topl_205[!topl_205$sub_locus %in% c(7, 10, 9, 8) & topl_205$final.locus==ll, ]$to_name))
}

topl_205

topl_205[!topl_205$sub_locus %in% c(7, 10, 9, 8), c('trait', 'to_name')]



qgraph(topl_205[!topl_205$sub_locus %in% c(7, 10, 9, 8), c('trait','to_name')], 
       directed=F, 
       layout='groups',
       groups=list(c(1,3,7,8,9), c(4,6,11,12), c(2,5,10)),
       border.color=c("#1F78B4" , "#E31A1C", "#1F78B4" ,"#33A02C", "#E31A1C" ,"#33A02C" ,
                      "#A6CEE3" ,"#A6CEE3", '#A6CEE3' ,"#FB9A99" ,"#B2DF8A", "#B2DF8A" ),
       # label.color= c("#1F78B4" , "#E31A1C", "#1F78B4" ,"#33A02C", "#E31A1C" ,"#33A02C" ,"#1F78B4" ,"#1F78B4", "#1F78B4" ,"#E31A1C" ,"#33A02C", "#33A02C" ),
       border.width= rep(2, 12),
       color=c("#1F78B4" , "#E31A1C", "#1F78B4" ,"#33A02C", "#E31A1C" ,"#33A02C" , "#A6CEE3" ,"#A6CEE3", '#A6CEE3' ,"#FB9A99" ,"#B2DF8A", "#B2DF8A" ),
       #edge.labels=topl_205[!topl_205$sub_locus %in% c(7, 10, 9, 8), 3], 
       edge.label.cex=0.5)



#---HLA multiple names problem --------------------
dim(loci.table)

length(loci.table[loci.table$trait=='f1',]$final.locus) #93
length(unique(loci.table[loci.table$trait=='f1',]$final.locus)) #92

length(loci.table[loci.table$trait=='f2',]$final.locus) #107
length(unique(loci.table[loci.table$trait=='f2',]$final.locus)) #103

length(loci.table[loci.table$trait=='f3',]$final.locus) #93
length(unique(loci.table[loci.table$trait=='f3',]$final.locus)) #91












