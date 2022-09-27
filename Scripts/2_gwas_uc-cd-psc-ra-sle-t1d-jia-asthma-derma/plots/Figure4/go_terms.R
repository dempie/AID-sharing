
library(RColorBrewer)
library(VennDiagram)
library(stringr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

go <- readRDS('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/10_genes_and_pathways/go_pathway_analysis.RDS')

#t cell activation venn diagram 

go_t <- go$result[go$result$term_name =='T cell activation',c('query', 'p_value', 'term_name', 'intersection')]

genes <- list()
for(i in 1:3){
  genes[[go_t$query[i]]] <- str_split(go_t$intersection[i], ',')[[1]]
}



#Make the plot
venn.diagram(
  x = genes,
  category.names = c('F1', 'F2', 'F3'),
  filename = 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure4/venn-diagram_go_t_cell_activationvenn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 1200 , 
  width = 1200 , 
  resolution = 300,
  compression = "lzw",
  lwd = 0.5,
  col=brewer.pal(6, 'Paired')[c(2,4,6)],
  fill = c(alpha(brewer.pal(6, 'Paired')[1],0.3), alpha(brewer.pal(6, 'Paired')[3],0.3), alpha(brewer.pal(6, 'Paired')[5],0.3)),
  cex = 0.4,
  fontfamily = "sans",
  cat.cex = 1,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 0),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = brewer.pal(6, 'Paired')[c(2,4,6)],
  rotation = 1
)


#lymphocyte activation venn diagram

go_t <- go$result[go$result$term_name =='lymphocyte activation',c('query', 'p_value', 'term_name', 'intersection')]

genes <- list()
for(i in 1:3){
  genes[[go_t$query[i]]] <- str_split(go_t$intersection[i], ',')[[1]]
}

#Make the plot
a <- venn.diagram(
  x = genes,
  category.names = c('F1', 'F2', 'F3'),
  filename = NULL,
  output = TRUE ,
  imagetype="png" ,
  height = 1200 , 
  width = 1200 , 
  resolution = 300,
  compression = "lzw",
  lwd = 0.5,
  col=brewer.pal(6, 'Paired')[c(2,4,6)],
  fill = c(brewer.pal(6, 'Paired')[1],brewer.pal(6, 'Paired')[3], brewer.pal(6, 'Paired')[5]),
  cex = 0.4,
  fontfamily = "sans",
  cat.cex = 1,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 0),
  cat.dist = c(0.055, 0.055, -0.45),
  cat.fontfamily = "sans",
  cat.col = brewer.pal(6, 'Paired')[c(2,4,6)],
  rotation = 1
)


pdf('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure4/venn-diagram_go_lymphocyte_activationvenn.pdf', width = 16, height = 16)
grid.draw(a)
dev.off()


#-------- barplot for shared genes top 10 go terms -----------------------------

#search the top 10 unique pathways and searcj the pvalues for all the factors for the top 10 pathways
top_10_pathways <- unique(go$result[order(go$result$p_value, decreasing = F),][,'term_name'])[(1:10)]
kg <- go$result[go$result$term_name %in% top_10_pathways,c('query', 'p_value', 'term_name', 'intersection')]
genes <- list()
to_merge <- list()
#loop for the pathways
for(i in 1:length(top_10_pathways)){
  pp <- top_10_pathways[i]
  
  #loop for the three facotors
  for(k in 1:3){
    tt <- c('f1','f2', 'f3')[k]
    genes[[pp]][[tt]]<- str_split(kg[kg$term_name==pp & kg$query==tt,]$intersection, ',')[[1]]
    
  }
  
  #make a table for each of the pathways
  shared_2 <- unique(c(intersect(genes[[pp]]$f1, genes[[pp]]$f2), intersect(genes[[pp]]$f1, genes[[pp]]$f3), intersect(genes[[pp]]$f2, genes[[pp]]$f3)))
  shared_3 <- Reduce(intersect,genes[[pp]]) 
  all_genes <- unlist(genes[[pp]])
  
  tabp <- data.frame('gene'=unique(all_genes), 'sharing'=rep(NA, length(unique(all_genes))), 'pathway'=rep(pp, length(unique(all_genes))))
  tabp[tabp$gene%in% shared_2,]$sharing <- 'two'
  
  if( length(shared_3)!=0 ){tabp[tabp$gene%in% shared_3,]$sharing <-  'three'} #the intersection between three could be 0 and cause an error
  
  tabp[!tabp$gene%in% c(shared_3, shared_2),]$sharing <- 'unique'
  to_merge[[pp]] <- tabp
  
  
}
final_t <- do.call(rbind, to_merge)



pdf('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure4/bar_plot_top10GOterm.pdf', width = 10, height = 7)
ggplot(final_t, aes(x=pathway, fill=sharing)) +
  geom_bar(stat='count',  color='black', position='fill')+
  scale_fill_manual(values = brewer.pal(8,'Set2')[c(2,1,8)]) +
  geom_text(aes(label = paste0("n=", ..count..)),position= 'fill',stat='count', color='Black')+
  labs(y = 'Number of genes', x = '')+
  theme_classic() + 
  theme(legend.position="bottom") + ggtitle('Top 10 GO Terms')
dev.off()



#test that the code is working
# a <- venn.diagram(
#   x = genes$`lymphocyte activation`  ,
#   category.names = c('F1', 'F2', 'F3'),
#   filename = NULL,
#   output = TRUE ,
#   imagetype="png" ,
#   height = 1200 ,
#   width = 1200 ,
#   resolution = 300,
#   compression = "lzw",
#   lwd = 0.5,
#   col=brewer.pal(6, 'Paired')[c(2,4,6)],
#   fill = c(brewer.pal(6, 'Paired')[1],brewer.pal(6, 'Paired')[3], brewer.pal(6, 'Paired')[5]),
#   cex = 0.4,
#   fontfamily = "sans",
#   cat.cex = 1,
#   cat.default.pos = "outer",
#   cat.pos = c(-27, 27, 0),
#   cat.dist = c(0.055, 0.055, -0.45),
#   cat.fontfamily = "sans",
#   cat.col = brewer.pal(6, 'Paired')[c(2,4,6)],
#   rotation = 1
# )
# 
# grid.draw(a)



#------------------------------------------------------------------

#GO TERMS top 5 each


top_10_pathways <- unique(go$result[order(go$result$p_value, decreasing = F),][,'term_name'])[(1:10)]
kg <- go$result[go$result$term_name %in% top_10_pathways,c('query', 'p_value', 'term_name', 'intersection')]


colnames(kg) <- c( 'trait', 'padj','pathway', 'genes')


kg <- kg %>%
  arrange(-desc(padj)) %>%
  group_by(trait) %>%
  slice(1:5)

kg$pathway_f <- paste(kg$pathway, kg$trait, sep = '_')
kg$pathway_f <- factor(kg$pathway_f )


pdf('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure4/barplot_top5eachF.pdf', height = 10, width = 12)
 
ggplot(kg, aes(x = factor(pathway_f, levels = rev(pathway_f)), y = -log10(padj))) +
    coord_flip() +
    geom_bar(stat = 'identity') +
    theme(axis.line.x = element_line(),
          axis.line.y = element_line(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          legend.position = "off")

dev.off()

