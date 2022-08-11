
library(RColorBrewer)
library(VennDiagram)
library(stringr)

go <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/go_pathway_analysis.RDS')

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
  filename = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure_enrichments/venn-diagram_go_t_cell_activationvenn.png',
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


pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure_enrichments/venn-diagram_go_lymphocyte_activationvenn.pdf', width = 16, height = 16)
grid.draw(a)
dev.off()



