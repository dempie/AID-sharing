library(data.table)
library(stringr)
library(gprofiler2)
library(biomaRt)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)


mart <- useDataset("hsapiens_gene_ensembl", useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")) #select the database to convert and maake sure it is build 37 


#----plot the GO terms heatmap -------------------------------------------------


go <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/go_pathway_analysis.RDS') 
he <- go$result
genes <- unique(strsplit(paste0(he$intersection, collapse = ','), split = ',')[[1]])
genes_symbol <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values=  genes,mart=mart )

genes
t_t <- data.frame(genes) 
for(i in 1:length(he$term_name)){
  t_t[, paste0(paste(strsplit(he$term_name[i], split=' ')[[1]], collapse = '_'),'_' , he$query[i])] <- as.numeric(genes %in% strsplit(he[he$term_name==he$term_name[i] & he$query==he$query[i], 'intersection'], split=',' )[[1]])
  
}


#add gene symbol and remove ENSEMBLE
rownames(t_t) <- genes_symbol$hgnc_symbol[match(t_t$genes, genes_symbol$ensembl_gene_id)]   
t_t$genes <- NULL


t_t['trait',] <- c(he$query) 
t_t['p_value',] <- c(he$p_value)
tt <- t(t_t)


#color selection
colori<- ComplexHeatmap:::default_col(x = tt)
colori[c('f1', 'f2', 'f3')] <-  RColorBrewer::brewer.pal(5, 'Paired')[c(1,3,5)]#colorblind safe
colori[c('0','1')] <- brewer.pal(9, 'Greys')[c(2,8)]

#column split
sep <- c(rep('a', nrow(t_t)-2), rep('b', 2))
names(sep) <- colnames(tt)





pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/supplementary/gene_enrichment/heatmaps/heatmap_GO_gene_and_pathways_only_path_genes.pdf', height = 20, width = 20)
Heatmap(tt, column_title = "Genes and pathways GO terms",
        rect_gp = gpar(col = "white", lwd = 0.25), 
        column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        left_annotation = rowAnnotation( '-log10(Padj)' = anno_barplot(   axis_param = list(direction = "reverse"),-log10(as.numeric(tt[,'p_value'])), width = unit(2, "cm"))),
        row_split = tt[,'trait'],
        border=T,
        column_split = sep,
        column_gap = unit(3, "mm"),
        row_gap = unit(5, "mm"), 
        row_names_gp = gpar(col =   RColorBrewer::brewer.pal(5, 'Paired')[c(1,3,5)], fontsize=2),
        column_names_gp = gpar(fontsize=3),
        column_names_side = 'bottom',
        show_heatmap_legend = F,
        row_title = "Pathways",
        na_col = 'red', 
        col = colori, 
        heatmap_width = unit(30, 'cm'), heatmap_height = unit(30, 'cm')
)
dev.off()


#----- plot go barplots --------------------------------------------------------

go <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/go_pathway_analysis.RDS')

top_10 <- go$result[order(go$result$p_value, decreasing = F),][1:20,c('query', 'p_value', 'term_name')]
top_10$p_log <- -log10(top_10$p_value)

pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/supplementary/gene_enrichment/heatmaps/go_top10_barplot.pdf', height = 16, width = 16)
ggplot(top_10, aes(x = term_name, y=p_log, fill=query)) + geom_col(position = position_dodge2(width = 0.9, preserve = "single"), col='grey') + 
  scale_fill_manual(values=brewer.pal(5,'Paired')[c(1,3,5)]) +
  ylab('-log10(P-adjusted)') + xlab('GO term')+ ggtitle('Top 20 GO terms')+
  coord_flip()  +theme_classic() 
dev.off()



#top 10 go terms circle plot

#search the top 10 unique pathways and searcj the pvalues for all the factors for the top 10 pathways
top_10_pathways <- unique(go$result[order(go$result$p_value, decreasing = F),][,'term_name'])[(1:10)]
got <- go$result[go$result$term_name %in% top_10_pathways,c('query', 'p_value', 'term_name')]

mtp <- matrix(nrow = length(unique(got$term_name)), ncol = 3)
colnames(mtp) <- c('f1', 'f2', 'f3')
rownames(mtp) <- unique(got$term_name)

for(i in c('f1', 'f2', 'f3')){
  ko <- got[ got$query==i,]$term_name
  mtp[got[ got$query==i, ]$term_name , i] <- -log10(got[ got$query==i, ]$p_value) #importat to set the radius of the circle
}


for(i in c('f1', 'f2', 'f3')){
  mtp[,i][is.na(mtp[,i])] <- 0
}

colnames(mtp) <- toupper(colnames(mtp))

pdf(height = 16, width = 16, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/supplementary/gene_enrichment/heatmaps/go_top_10_circle_heatmpat.pdf')
Heatmap(mtp,  
        rect_gp = gpar(col = brewer.pal(5, 'Greys')[2], lwd = 0.25), 
        column_title = 'Top 10 GO terms',
        cluster_rows = F,
        cluster_columns = F,
        col = c('white', 'white'), 
        border=T,
        show_column_dend = F,
        show_row_dend = F,
        row_names_gp = gpar(fontsize=10),
        row_names_side = 'left',
        width = unit(10, 'cm'),
        height = unit(30, 'cm'),
        show_heatmap_legend = F,
        column_names_gp = gpar(fontsize=25),
        column_names_side = 'bottom',
        column_names_rot = 360,
        top_annotation =  columnAnnotation(legend = anno_empty(border = F,
                                                               width =  unit(10, "cm"), height = unit(2, 'cm'))), #space for the legend annotation
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(mtp[i, j] > 0){
            grid.circle(x = x, y = y, r = unit(mtp[i,j], 'mm'), gp = gpar(fill = brewer.pal(7,'BrBG')[4], col = 'black'))
            #grid.text(signif(mtp[i,j],2), x=x,  y = y, gp=gpar(fontsize=10))
          }
        }
        
)


for(k in 1:5) {
  i <- c(0.05,0.30,0.50, 0.70, 0.80)[k] #this are the positions 
  
  decorate_annotation('legend', {
    
    grid.circle(x=i, r= unit(-log10(c(10^-15,10^-12,10^-9,10^-6,10^-3)), 'mm')[k],gp = gpar(fill = brewer.pal(7,'BrBG')[3], col = 'black')) #this are the pvalues that are shown in the legend, the dimension of the -log10(p) is in mm
    grid.text(x=i, round(-log10(c(10^-15,10^-12,10^-9,10^-6,10^-3)), 3)[k]) #show the p
  })
}
dev.off()





#----plot the GO terms heatmap -------------------------------------------------


react <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/react_pathway_analysis.RDS') 
he <- react$result
genes <- unique(strsplit(paste0(he$intersection, collapse = ','), split = ',')[[1]])
genes_symbol <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values=  genes,mart=mart )

t_t <- data.frame(genes) 
for(i in 1:length(he$term_name)){
  t_t[, paste0(paste(strsplit(he$term_name[i], split=' ')[[1]], collapse = '_'),'_' , he$query[i])] <- as.numeric(genes %in% strsplit(he[he$term_name==he$term_name[i] & he$query==he$query[i], 'intersection'], split=',' )[[1]])
  
}


#add gene symbol and remove ENSEMBLE
rownames(t_t) <- genes_symbol$hgnc_symbol[match(t_t$genes, genes_symbol$ensembl_gene_id)]   
t_t$genes <- NULL


t_t['trait',] <- c(he$query) 
t_t['p_value',] <- c(he$p_value)
tt <- t(t_t)


#color selection
colori<- ComplexHeatmap:::default_col(x = tt)
colori[c('f1', 'f2', 'f3')] <-  RColorBrewer::brewer.pal(5, 'Paired')[c(1,3,5)]#colorblind safe
colori[c('0','1')] <- brewer.pal(9, 'Greys')[c(2,8)]

#column split
sep <- c(rep('a', nrow(t_t)-2), rep('b', 2))
names(sep) <- colnames(tt)


pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/supplementary/gene_enrichment/heatmaps/heatmap_react_gene_and_pathways_only_path_genes.pdf', height = 20, width = 20)
Heatmap(tt, column_title = "Genes and pathways GO terms",
        rect_gp = gpar(col = "white", lwd = 0.25), 
        column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        left_annotation = rowAnnotation( '-log10(Padj)' = anno_barplot(   axis_param = list(direction = "reverse"),-log10(as.numeric(tt[,'p_value'])), width = unit(2, "cm"))),
        row_split = tt[,'trait'],
        border=T,
        column_split = sep,
        column_gap = unit(3, "mm"),
        row_gap = unit(5, "mm"), 
        row_names_gp = gpar(col =   RColorBrewer::brewer.pal(5, 'Paired')[c(1,3,5)], fontsize=8),
        column_names_gp = gpar(fontsize=8),
        column_names_side = 'bottom',
        show_heatmap_legend = F,
        row_title = "Pathways",
        na_col = 'red', 
        col = colori
)
dev.off()


#-----heatmap clusterprofiler---------------------------------------------------
# 
# go <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/go_pathway_analysis.RDS')
# 
# 
# 
# #----clusterprofiler go terms -------------------------------------------------
# loci.table <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/loci_table_names_nicola_genes.csv', data.table = F)
# mart <- useDataset("hsapiens_gene_ensembl", useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")) #select the database to convert and maake sure it is build 37
# 
# 
# 
# f_list_genes <- list()
# for(i in 1:3){
#   tt <- c('f1', 'f2', 'f3')[i]
#   f_list_genes[[tt]] <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values= loci.table[loci.table$trait==tt, ]$closest_gene ,mart= mart)
# 
# }
# 
# 
# aaa <-list(f1=as.character(unique(f_list_genes$f1$entrezgene_id)), f2=as.character(unique(f_list_genes$f2$entrezgene_id)), f3=as.character(unique(f_list_genes$f3$entrezgene_id)))
# #bbb <-list(f1=as.character(unique(f_list_genes$f1$ensembl_gene_id)), f2=as.character(unique(f_list_genes$f2$ensembl_gene_id)), f3=as.character(unique(f_list_genes$f3$ensembl_gene_id))) ensemble gene name
# 
# 
# test_enrich <- compareCluster(geneCluster = aaa,
#                              ont           = "BP",
#                               OrgDb = org.Hs.eg.db,
#                               #pAdjustMethod = "BH",
#                               pvalueCutoff  = 0.01,
#                               qvalueCutoff  = 0.05,
#                               readable      = TRUE,
#                               fun =  enrichGO)
# 
# 
# 
# dotplot(test_enrich, showCategory=5)
# 
# 
# 
# 
# 
# 
# go_order <- go$result[order(go$result$p_value), c('query', 'p_value', 'term_name')]
# go_order$logp <- -log10(go_order$p_value)
# 
# 
# a <- list()
# for(i in 1:3){
#   tt <- c('f1', 'f2', 'f3')[i]
#   a[[tt]]<- go_order[go_order$query==tt,][1:5,]
#   a[[tt]]$top <- rep('top', 5)
# }
# toplot<- Reduce(rbind, a)
# 
# ggplot(toplot,aes(x=query,y=term_name,color=logp)) +geom_point(size=10)+scale_colour_gradient(low="red",high="blue")
# 
# 
# ###
# 
# 
# to_plot_shape <- go_order[go_order$term_name %in% toplot$term_name,]
# to_plot_shape$top <- rep('n', nrow(to_plot_shape))
# for(i in 1:nrow(toplot)){
#   to_plot_shape[to_plot_shape$query==toplot$query[i] & to_plot_shape$term_name == toplot$term_name[i], ]$top <- 'top'
# }
# 
# 
# ggplot(to_plot_shape,aes(x=query,y=term_name,color=logp, shape=top)) +geom_point(size=10) +scale_colour_gradient(low="red",high="blue") + theme_light()
# 
# 






