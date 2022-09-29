library(data.table)
library(ComplexHeatmap)
library(RColorBrewer)
mr_new <- fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/tables/Supplementary_table_MR_results.csv')

mr_f <- mr_new[mr_new$TRAIT %in% c('Falrg', 'Faid', 'Fgut') & mr_new$PVAL<0.05,]

mr_f$effect <- ifelse(mr_f$B>0, 1, -1)

plot_f <- mr_f[, c(5,6,7,9,8)]
duplicated(plot_f[,c(1,2,3,4)])

for(i in 1:nrow(plot_f)){
  plot_f$name[i] <-paste(plot_f$GENE[i], plot_f$`LOCUS_NAME_(CHR:START_END)`[i]) 
}


plot_f



mat <-matrix(nrow =length(unique(plot_f$name)), ncol = length(unique(plot_f$CELL_TYPE))) 

rownames(mat) <- unique(plot_f$name)
colnames(mat) <- unique(plot_f$CELL_TYPE)


for(i in rownames(mat)){
  
  for(k in colnames(mat)){
    
    mat[i,k] <- ifelse(length(plot_f[plot_f$name==i & plot_f$CELL_TYPE==k,]$effect)!=0, paste0(plot_f[plot_f$name==i & plot_f$CELL_TYPE==k , c( 'effect', 'TRAIT')], collapse = ''), NA)
    
  }
}

names(table(mat))
colori <- ComplexHeatmap:::default_col(mat)
colori[names(table(mat))] <-  rep('white', 7)

for(i in 1:nrow(mat)){
  rownames(mat)[i] <- strsplit(rownames(mat)[i],' ')[[1]][1]
}

colnames(mat)
mat<- mat[, c('CD4NC', 'CD4ET', 'CD4SOX4', 'CD8NC', 'CD8ET', 'CD8S100B','BIN', 'BMem', 'Plasma', 'MonoC', 'MonoNC', 'DC' ,'NK', 'NKR')]

pdf('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure5/mr_heatmap.pdf', height = 30, width=8)
Heatmap(mat, 
        #column_title = "Genes and pathways KEGG",
        rect_gp = gpar(col = brewer.pal(5, 'Greys')[2], lwd = 0.25), 
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        col=colori, 
        border=T,
        row_names_gp = gpar( fontsize=8),
        column_names_gp = gpar(fontsize=8),
        column_names_side = 'bottom',
        row_names_side = 'left',
        show_heatmap_legend = F,
        row_title = "Pathways",
        na_col = 'white', 
        height = unit(30, 'cm'), width = 8, 
        
        
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(!is.na(mat[i, j] )){
            
            if(mat[i, j]== '-1Faid' ){
              grid.points(x = x, y = y, pch = 25, size = unit(5, "mm"), gp = gpar(col = brewer.pal(6, 'Paired')[4], fill=brewer.pal(5, 'Paired')[4] ))
            }
            
            
            if(mat[i, j]== '1Faid' ){
              grid.points(x = x, y = y, pch = 24, size = unit(5, "mm"), gp = gpar(col = brewer.pal(5, 'Paired')[4], fill=brewer.pal(5, 'Paired')[4]))
            }
            
            if(mat[i, j]== '-1Falrg' ){
              grid.points(x = x, y = y, pch = 25, size = unit(5, "mm"), gp = gpar(col = brewer.pal(6, 'Paired')[6], fill=brewer.pal(6, 'Paired')[6]))
            }
            
            if(mat[i, j]== '1Falrg' ){
              grid.points(x = x, y = y, pch = 24, size = unit(5, "mm"), gp = gpar(col = brewer.pal(6, 'Paired')[6], fill=brewer.pal(6, 'Paired')[6]))
            }
            
            if(mat[i, j]== '1Fgut' ){
              grid.points(x = x, y = y, pch = 24, size = unit(5, "mm"), gp = gpar(col = brewer.pal(6, 'Paired')[2], fill=brewer.pal(6, 'Paired')[2]))
            }
            
            if(mat[i, j]== '-1Fgut' ){
              grid.points(x = x, y = y, pch = 25, size = unit(5, "mm"), gp = gpar(col = brewer.pal(6, 'Paired')[2], fill=brewer.pal(6, 'Paired')[2]))
            }
            
            if(mat[i, j]== 'c(-1, -1)c(\"Fgut\", \"Faid\")' ){
              grid.points(x = x, y = y, pch=25, size = unit(5, "mm"), gp = gpar(col = 'grey', fill='yellow'))
            }
            
          }
          
        }
)

dev.off()

