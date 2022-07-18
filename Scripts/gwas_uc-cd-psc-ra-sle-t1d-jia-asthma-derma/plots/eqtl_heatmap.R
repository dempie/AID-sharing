library(data.table)
library(stringr)
library(ComplexHeatmap)
library(RColorBrewer)

eqtl <- read.csv2('Nicola/sc_eqtl/All_colocalization_eqtl.tsv', sep = '\t')

pp_ok <- eqtl[as.numeric(eqtl$PP.H4.abf)>0.9,]
factor_eqtl<- pp_ok[pp_ok$t1%in%c('f1', 'f2', 'f3'), ]

barplot(table(factor_eqtl$t1), main = 'eQTL per Factor')





genes <- unique(factor_eqtl[order(factor_eqtl$t1),]$t2) #order by factor 
cell_types <- unique(factor_eqtl$cell_type)

semaph <- matrix(nrow = length(cell_types), ncol=length(genes), rep(0, length(cell_types)*length(genes)) )
colnames(semaph) <- genes
rownames(semaph) <- cell_types 


for(i in 1:nrow(semaph)){
  
  for(j in 1:ncol(semaph)){
    semaph[i,j] <- paste0(unique(factor_eqtl[factor_eqtl$cell_type==rownames(semaph)[i] &factor_eqtl$t2==colnames(semaph)[j], ]$t1) , collapse = '-')
    
  }
 
  
}

#how many different entries
unique(as.vector(semaph))

#order the cell types
semaph <- semaph[c('CD4NC', 'CD4ET', 'CD4SOX4', 'CD8NC', 'CD8ET', 'CD8S100B', 'BMem', 'BIN','Plasma', 'NK', 'NKR','MonoC','MonoNC','DC'),]
dim(semaph)


#prepare row annotation
a <- data.frame()
ccc <- matrix(nrow = length(cell_types), ncol=6, rep(0, length(cell_types)*6) )
colnames(ccc) <- unique(as.vector(semaph))
rownames(ccc) <- rownames(semaph)

for(i in 1:nrow(semaph)){
  a <- table(semaph[i,])
  k <- names(a)[str_detect(names(a), '')]
  ccc[i,k] <- a[k]
}

#select colors

colori<- ComplexHeatmap:::default_col(semaph)
colori[c('f1', 'f2', 'f3', 'f1-f2', 'f3-f1')] <-  rep('white', 5)

pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure_eqtl/semaphor_qtl.pdf', width = 24, height = 6.33)
Heatmap(semaph, 
           column_title = "eQTL genes and Cell types factors",
           rect_gp = gpar(col = brewer.pal(5, 'Greys')[2], lwd = 0.25), 
           column_title_gp = gpar(fontsize = 20, fontface = "bold"),
           col=colori, 
           border=T,
           row_names_gp = gpar( fontsize=10), row_names_side = 'left',
           column_names_gp = gpar(fontsize=10),
           column_names_side = 'bottom',
           show_heatmap_legend = F,
           row_title = "Cell types",
           na_col = 'white', 
           heatmap_width = unit(24, 'cm'), heatmap_height = unit(10, 'cm'),
        #barplots
           right_annotation =  rowAnnotation('eQTL Genes' = anno_barplot(t(apply(ccc, 1, function(x) x/sum(x)))[, c('f1', 'f2', 'f3', 'f1-f2', 'f3-f1')], 
                                                                         gp=gpar(fill= c(brewer.pal(10, 'Paired')[c(2,4,6)], 'yellow', 'black'), 
                                                                        col=rep('grey', 6)), 
                                                                         bar_width = .7, 
                                                                         width = unit(3, "cm"))) ,
           
           
           
           cell_fun = function(j, i, x, y, width, height, fill) {
             if(!is.na(semaph[i, j] )){
               
               if(semaph[i, j]== 'f1' ){
                 grid.points(x = x, y = y, pch = 22, size = unit(4, "mm"), gp = gpar(col = brewer.pal(6, 'Paired')[1], fill=brewer.pal(5, 'Paired')[2] ))
               }
               
               
               if(semaph[i, j]== 'f2' ){
                 grid.points(x = x, y = y, pch = 22, size = unit(4, "mm"), gp = gpar(col = brewer.pal(5, 'Paired')[3], fill=brewer.pal(5, 'Paired')[4]))
               }
               
               if(semaph[i, j]== 'f3' ){
                 grid.points(x = x, y = y, pch = 22, size = unit(4, "mm"), gp = gpar(col = brewer.pal(6, 'Paired')[5], fill=brewer.pal(6, 'Paired')[6]))
               }
               
               if(semaph[i, j]== 'f1-f2' ){
                 grid.points(x = x, y = y, pch=22, size = unit(4, "mm"), gp = gpar(col = 'grey', fill='yellow'))
               }
               
               if(semaph[i, j]== 'f3-f1' ){
                 grid.points(x = x, y = y, pch=22, size = unit(4, "mm"), gp = gpar(col = 'grey', fill='black'))
               }
             }
             
           }
)

dev.off()























