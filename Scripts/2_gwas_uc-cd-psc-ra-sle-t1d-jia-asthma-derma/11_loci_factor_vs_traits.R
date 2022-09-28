library(data.table)
library(ComplexHeatmap)
library(RColorBrewer)
library(qgraph)
library(ggplot2)

#------ load the dataset -------------------------------------------------------
loci.table <- fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/09_loci_definitions_Nicola/loci_definitions/final_locus_table.tsv', data.table = F) #take the original factor loci list
loci.table$final.locus=paste0(loci.table$Chr,":",loci.table$start,"-",loci.table$end,"_",loci.table$sub_locus)

#----- make an upset plot-------------------------------------------------------
q <- make_comb_mat(list(f1=loci.table[loci.table$trait=='f1',]$final.locus, cd=loci.table[loci.table$trait=='cd',]$final.locus, uc=loci.table[loci.table$trait=='uc',]$final.locus, psc=loci.table[loci.table$trait=='psc',]$final.locus), mode = 'distinct')

pdf(width = 9, height = 5, file = 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/09_loci_definitions_Nicola/upset_plot_final_loci_nicola.pdf')
UpSet(q, set_order = c("f1", "f2", "f3"), 
      comb_order = c(5,6,7,2,3,4,1),
      comb_col = c(brewer.pal(8, 'Set2')[c(7,6,6,6)],brewer.pal(12, 'Paired')[c(2,4,6)]),
      top_annotation = upset_top_annotation(q, add_numbers = TRUE, height = unit(6, "cm")),
      right_annotation = upset_right_annotation(q, add_numbers = TRUE, width = unit(5,'cm'),gp = gpar(fill = brewer.pal(5, 'Greys')[1])  ),
      row_title = "Factor", 
      column_title = "Intersection of all loci, colocalization")
dev.off()





###bar plot new loci compared to old version

#----f1-------------------------------------------------------------------------
UpSet(q, top_annotation = upset_top_annotation(q, add_numbers = TRUE, height = unit(6, "cm")),
      right_annotation = upset_right_annotation(q, add_numbers = TRUE, width = unit(5,'cm'),gp = gpar(fill = brewer.pal(5, 'Greys')[1])  ))


f1 <- unique(loci.table[loci.table$trait=='f1',]$final.locus)
cd <-  unique(loci.table[loci.table$trait=='cd',]$final.locus)
uc <- unique(loci.table[loci.table$trait=='uc',]$final.locus)
psc <- unique(loci.table[loci.table$trait=='psc',]$final.locus)




q <- make_comb_mat(list(f1=loci.table[loci.table$trait=='f1',]$final.locus, 
                        cd=loci.table[loci.table$trait=='cd',]$final.locus, 
                        uc=loci.table[loci.table$trait=='uc',]$final.locus, 
                        psc=loci.table[loci.table$trait=='psc',]$final.locus), 
                   mode = 'distinct')
    
q

f1_uc <- extract_comb(q, "1010") #1
f1_psc <- extract_comb(q, "1001")
f1_cd_uc <- extract_comb(q, '1110') #2
f1_cd <- extract_comb(q, '1100') #1
f1_cd_uc_psc <- extract_comb(q, '1111') #3
new_f1 <- extract_comb(q, '1000') #new


f1_1 <- length(unique(c(f1_cd, f1_uc, f1_psc)))
f1_2 <- length(c(f1_cd_uc))
f1_3 <- length(c(f1_cd_uc_psc))
f1_new <- length(new_f1)

prova <- data.frame(c(f1_new,f1_1, f1_2, f1_3))
rownames(prova) <- c('Fgut', '1', '2', '3')
colnames(prova) <- 'number'


pdf('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/11_loci_factors_vs_traits/fgut_barplot.pdf')
ggplot(prova,aes(rownames(prova), number, label=number)) +
  geom_bar(stat = 'identity') + 
  geom_text()+
  theme_classic()
dev.off()





#---------f2--------------------------------------------------------------------

f2 <- unique(loci.table[loci.table$trait=='f2',]$final.locus)
t1d <-  unique(loci.table[loci.table$trait=='t1d',]$final.locus)
sle <- unique(loci.table[loci.table$trait=='sle',]$final.locus)
jia <- unique(loci.table[loci.table$trait=='jia',]$final.locus)
ra <- unique(loci.table[loci.table$trait=='ra',]$final.locus)

f2_q <- make_comb_mat(f2=f2, t1d=t1d, sle=sle,jia=jia, ra=ra)

UpSet(f2_q, top_annotation = upset_top_annotation(f2_q, add_numbers = TRUE, height = unit(6, "cm")),
      right_annotation = upset_right_annotation(f2_q, add_numbers = TRUE, width = unit(5,'cm'),gp = gpar(fill = brewer.pal(5, 'Greys')[1])  ))



#new f2
f2 <- length(unique(extract_comb(f2_q, '10000')))

#1
f2_t1d <- extract_comb(f2_q, '11000')
f2_sle <- extract_comb(f2_q, '10100')
f2_jia <- extract_comb(f2_q, '10010') #empty
f2_ra <- extract_comb(f2_q, '10001')


f2_1 <-length(unique( c(f2_t1d, f2_sle, f2_jia, f2_ra))) #53

#2
f2_t1d_sle <- extract_comb(f2_q, '11100') #empty
f2_t1d_jia <- extract_comb(f2_q, '11010') #
f2_t1d_ra <- extract_comb(f2_q, '11001')
f2_sle_jia <- extract_comb(f2_q, '10110')
f2_sle_ra <- extract_comb(f2_q, '10101')
f2_ra_jia <- extract_comb(f2_q, '10011')

f2_2 <- length(unique(c(f2_t1d_sle, f2_t1d_jia,f2_t1d_ra ,f2_sle_jia ,f2_sle_ra ,f2_ra_jia))) #13

#3
f2_t1d_sle_jia <- extract_comb(f2_q, '11110')
f2_t1d_sle_ra <- extract_comb(f2_q, '11101')
f2_t1d_ra_jia <- extract_comb(f2_q, '11011')
f2_ra_jia_sle <- extract_comb(f2_q, '10111')

f2_3 <- length(unique(c(f2_t1d_sle_jia,f2_t1d_sle_ra, f2_t1d_ra_jia, f2_ra_jia_sle  ))) #none




faid <- data.frame(c(f2,f2_1, f2_2, f2_3))
rownames(faid) <- c('Faid', '1', '2','3')
colnames(faid) <- 'number'


pdf('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/11_loci_factors_vs_traits/faid_barplot.pdf')
ggplot(faid,aes(rownames(faid), number, label=number)) +
  geom_bar(stat = 'identity') + 
  geom_text()+
  theme_classic()
dev.off()



#------f3-----------------------------------------------------------------------

f3 <- unique(loci.table[loci.table$trait=='f3',]$final.locus)
asthma <-  unique(loci.table[loci.table$trait=='asthma',]$final.locus)
derma <- unique(loci.table[loci.table$trait=='derma',]$final.locus)



UpSet(f3_q, top_annotation = upset_top_annotation(f3_q, add_numbers = TRUE, height = unit(6, "cm")),
      right_annotation = upset_right_annotation(f3_q, add_numbers = TRUE, width = unit(5,'cm'),gp = gpar(fill = brewer.pal(5, 'Greys')[1])  ))




f3_q <- make_comb_mat(f3=f3, asthma=asthma, derma=derma)

#new f3
f3 <-length(unique( extract_comb(f3_q, '100')))

#f2 2
f3_asth <- extract_comb(f3_q, '110')
f3_derm<- extract_comb(f3_q, '101')

f3_2 <- length(unique(c(f3_asth, f3_derm)))

#f3 3
f3_3 <- length(unique(extract_comb(f3_q, '111')))





falrg <- data.frame(c(f3,f3_2, f3_3))
rownames(falrg) <- c('Falrg', '1', '2' )
colnames(falrg) <- 'number'


pdf('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/11_loci_factors_vs_traits/falrg_barplot.pdf')
ggplot(falrg,aes(rownames(falrg), number, label=number)) +
  geom_bar(stat = 'identity') + 
  geom_text()+
  theme_classic()
dev.off()














  
  
  

  
  
 