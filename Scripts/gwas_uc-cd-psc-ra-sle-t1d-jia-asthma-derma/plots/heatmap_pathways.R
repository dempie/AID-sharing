library(data.table)
library(ChIPpeakAnno)
library(stringr)
library(gprofiler2)
library(GenomicRanges)
library(EnsDb.Hsapiens.v75)
library(biomaRt)
library(ComplexHeatmap)
library(gtools)
library(RColorBrewer)
library(formattable)
library(ggplot2)

mart <- useDataset("hsapiens_gene_ensembl", useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl"))
loci.table <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/loci_table_names_nicola_genes.csv', data.table = F)

kegg <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/kegg_pathway_analysis.RDS')
go <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/go_pathway_analysis.RDS')

#-----heatmap kegg results------------------------------------------------------

he <- kegg$result
genes <- unique(strsplit(paste0(he$intersection, collapse = ','), split = ',')[[1]])
genes_symbol <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values=  genes,mart=mart )
genes_symbol <- genes_symbol[!duplicated(genes_symbol$ensembl_gene_id), ]


t_t <- data.frame(genes) 
for(i in 1:length(he$term_name)){
  t_t[, paste0(paste(strsplit(he$term_name[i], split=' ')[[1]], collapse = '_'),'_' , he$query[i])] <- as.numeric(genes %in% strsplit(he[he$term_name==he$term_name[i] & he$query==he$query[i], 'intersection'], split=',' )[[1]])
  
}

t_t
#add gene symbol and remove ENSEMBLE
rownames(t_t) <- genes_symbol$hgnc_symbol[match(t_t$genes, genes_symbol$ensembl_gene_id)]   
t_t$genes <- NULL
t_t

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

ht_opt$ROW_ANNO_PADDING = unit(0.2, "cm")
#plot heatmap_gene_and_pathways.pdf             
pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure_pathways/heatmap_gene_and_pathways_only_path_genes.pdf', height = 8, width = 16)
Heatmap(tt, column_title = "Genes and pathways KEGG genes",
        rect_gp = gpar(col = "white", lwd = 0.25, type='2'), 
        column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        left_annotation = rowAnnotation( '-log10(Padj)' = anno_barplot(   axis_param = list(direction = "reverse"),-log10(as.numeric(tt[,'p_value'])), width = unit(2, "cm"))),
        row_split = tt[,'trait'],
        column_split = sep,
        column_gap = unit(2, "mm"),
        border=T,
        row_gap = unit(5, "mm"), 
        row_names_gp = gpar(col =   RColorBrewer::brewer.pal(5, 'Paired')[c(1,3,5)], fontsize=8),
        column_labels = c(colnames(tt)[1:c(nrow(t_t)-2)], 'Factor', '') ,
        column_names_gp = gpar(fontsize=7),
        column_names_side = 'bottom',
        show_heatmap_legend = F,
        row_title = "Pathways",
        na_col = 'red', 
        col = colori, 
        heatmap_width = unit(28, 'cm'), heatmap_height = unit(12, 'cm')
)
dev.off()


#---- plot the heatmpat with all genes -----------------------------------------

he <- kegg$result
genes <- unique(loci.table[order(loci.table$trait),]$closest_gene)
genes_symbol <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values=  genes,mart=mart )
genes_symbol <- genes_symbol[!duplicated(genes_symbol$ensembl_gene_id), ]


t_t <- data.frame(genes) 
for(i in 1:length(he$term_name)){
  t_t[, paste0(paste(strsplit(he$term_name[i], split=' ')[[1]], collapse = '_'),'_' , he$query[i])] <- as.numeric(genes %in% strsplit(he[he$term_name==he$term_name[i] & he$query==he$query[i], 'intersection'], split=',' )[[1]])
  
}

genes_symbol[!str_detect(genes_symbol$hgnc_symbol,'' ),]$hgnc_symbol <- genes_symbol[!str_detect(genes_symbol$hgnc_symbol,'' ),]$ensembl_gene_id
#add gene symbol and remove ENSEMBLE
rownames(t_t) <- genes_symbol$hgnc_symbol[match(t_t$genes, genes_symbol$ensembl_gene_id)]   
t_t$genes <- NULL
t_t

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

ht_opt$ROW_ANNO_PADDING = unit(0.2, "cm")
#plot heatmap_gene_and_pathways.pdf         


#plot heatmap_gene_and_pathways_with_all_genes.pdf

for(i in 1:nrow(tt)){
    rownames(tt)[i] <- str_replace(string = rownames(tt)[i],pattern = '_', replacement = ' ' )
    rownames(tt)[i] <- str_remove(rownames(tt)[i], pattern = '_f1')
    rownames(tt)[i] <- str_remove(rownames(tt)[i], pattern = '_f2')
    rownames(tt)[i] <- str_remove(rownames(tt)[i], pattern = '_f3')
}

pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure_pathways/heatmap_gene_and_pathways_all_genes_full_matrix.pdf', height = 8, width = 16)
Heatmap(tt, column_title = "Genes and pathways all genes KEGG",
        rect_gp = gpar(col = "white", lwd = 0.25, type='2'), 
        column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        left_annotation = rowAnnotation( '-log10(Padj)' = anno_barplot(   axis_param = list(direction = "reverse"),-log10(as.numeric(tt[,'p_value'])), width = unit(2, "cm"))),
        row_split = tt[,'trait'],
        column_split = sep,
        column_gap = unit(2, "mm"),
        border=T,
        row_gap = unit(5, "mm"), 
        row_names_gp = gpar(col =   RColorBrewer::brewer.pal(5, 'Paired')[c(1,3,5)], fontsize=10),
        column_labels = c(colnames(tt)[1:c(nrow(t_t)-2)], 'Factor', '') ,
        column_names_gp = gpar(fontsize=0),
        column_names_side = 'bottom',
        show_heatmap_legend = F,
        row_title = "Pathways",
        na_col = 'red', 
        col = colori, 
        heatmap_width = unit(30, 'cm'), heatmap_height = unit(14, 'cm')
)
dev.off()





#-------------------heatmpap and go terms------------------------------------

he <- go$result
genes <- unique(strsplit(paste0(he$intersection, collapse = ','), split = ',')[[1]])
genes_symbol <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values=  genes,mart=mart )
genes_symbol <- genes_symbol[!duplicated(genes_symbol$ensembl_gene_id), ]

t_t <- data.frame(genes) 
for(i in 1:length(he$term_name)){
  t_t[, paste0(paste(strsplit(he$term_name[i], split=' ')[[1]], collapse = '_'),'_' , he$query[i])] <- as.numeric(genes %in% strsplit(he[he$term_name==he$term_name[i] & he$query==he$query[i], 'intersection'], split=',' )[[1]])
  
}

t_t
#add gene symbol and remove ENSEMBLE
rownames(t_t) <- genes_symbol$hgnc_symbol[match(t_t$genes, genes_symbol$ensembl_gene_id)]   
t_t$genes <- NULL
t_t

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



#----plot the GO terms heatmap -------------------------------------------------

pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure_pathways/heatmap_GO_gene_and_pathways_only_path_genes.pdf', height = 20, width = 20)
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



#-------------------table of pathways-------------------------------------------

#kegg table
kg <- kegg$result[, c('query','p_value', 'term_name', 'intersection')] 

#use gene symbols
he <- kegg$result
genes <- unique(strsplit(paste0(he$intersection, collapse = ','), split = ',')[[1]])
genes_symbol <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values=  genes,mart=mart )
genes_symbol <- genes_symbol[!duplicated(genes_symbol$ensembl_gene_id), ]

kg$hgnc <- rep('NA', nrow(kg))
for(i in 1:nrow(kg)){
  kg$hgnc[i] <- paste0(genes_symbol[ genes_symbol$ensembl_gene_id %in% strsplit( kg$intersection[i], split = ',')[[1]] ,]$hgnc_symbol, collapse = ',  ')
}

colnames(kg) <- c('Trait', 'P-value (adjusted)', 'Pathway', 'ENSEMBL', 'Gene symbols')
kg$Trait <- toupper(kg$Trait)
plott <- kg[,c(3,1,2, 5)]
plott$`P-value (adjusted)` <- signif(plott$`P-value (adjusted)`, 2)



formattable(plott,align =c("l","c","c" ,"r") ,
            list( 'Pathway' = formatter( "span", style = style(color = "grey",font.weight = "bold")))
)

fwrite(plott, file = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure_pathways/kegg_results.csv',sep = ',', col.names = T, quote = F, row.names = F)




#go table
gg <-  go$result[, c('query','p_value', 'term_name', 'intersection')] 
gg <- head(gg[order(gg$p_value), ],20)
genes <- unique(strsplit(paste0(gg$intersection, collapse = ','), split = ',')[[1]])
genes_symbol <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values=  genes,mart=mart )
genes_symbol <- genes_symbol[!duplicated(genes_symbol$ensembl_gene_id), ]
gg$hgnc <- rep('NA', nrow(gg))

for(i in 1:nrow(gg)){
  gg$hgnc[i] <- paste0(genes_symbol[ genes_symbol$ensembl_gene_id %in% strsplit( gg$intersection[i], split = ',')[[1]] ,]$hgnc_symbol, collapse = ',  ')
}

colnames(gg) <- c('Trait', 'P-value (adjusted)', 'GO term', 'ENSEMBL', 'Gene symbols')
gg$Trait <- toupper(gg$Trait)
plott <- gg[,c(3,1,2, 5)]
plott$`P-value (adjusted)` <- signif(plott$`P-value (adjusted)`, 2)

rownames(plott) <- NULL
formattable(plott,align =c("l","c","c" ,"r") ,
            list( 'GO term' = formatter( "span", style = style(color = "grey",font.weight = "bold")))
)


fwrite(plott, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure_pathways/go_results.csv',sep = ',', col.names = T, quote = F, row.names = F)



################################################################################

#----------- KEGG semaphor------------------------------------------------------

he <- kegg$result
genes <- unique(strsplit(paste0(he$intersection, collapse = ','), split = ',')[[1]])
genes_symbol <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values=  genes,mart=mart )

#ecxclude the duplicated ones
genes_symbol <- genes_symbol[!duplicated(genes_symbol$ensembl_gene_id), ]


t_t <- data.frame(genes) 
for(i in 1:length(he$term_name)){
  t_t[, paste0(paste(strsplit(he$term_name[i], split=' ')[[1]], collapse = '_'),'_' , he$query[i])] <- as.numeric(genes %in% strsplit(he[he$term_name==he$term_name[i] & he$query==he$query[i], 'intersection'], split=',' )[[1]])
  
}

#add gene symbol and remove ENSEMBLE
rownames(t_t) <- genes_symbol$hgnc_symbol[match(t_t$genes, genes_symbol$ensembl_gene_id)]   
t_t$genes <- NULL

#--------------
kg <- kegg$result[,c('query', 'p_value', 'term_name')]

mtp <- matrix(nrow = length(unique(kg$term_name)), ncol = 3)
colnames(mtp) <- c('f1', 'f2', 'f3')
rownames(mtp) <- unique(kg$term_name)

        for(i in c('f1', 'f2', 'f3')){
          kg[ kg$query==i,]$term_name
          mtp[, i][base::match(  kg[ kg$query==i, ]$term_name, names(mtp[, i]))] <- -log10(kg[ kg$query==i, ]$p_value) #importat to set the radius of the circle
        }
        
        
        for(i in c('f1', 'f2', 'f3')){
          mtp[,i][is.na(mtp[,i])] <- 0
        }
        
        colnames(mtp) <- toupper(colnames(mtp))
        
        #order mtp with the shared pathways on top and the genes to be clustereb by factors
        otp <- mtp[c(1,3,5,8,9,4,2,6,7,10,11,12,13,14,15),]

       

#add the genes

gtp <- matrix(nrow = length(unique(kg$term_name)), ncol = nrow(genes_symbol))
rownames(gtp) <- unique(kg$term_name)
colnames(gtp) <- genes_symbol$hgnc_symbol
kres <- kegg$result



g_s <- list()
u <- vector()
for(k in 1:length(unique(kres$term_name)) ){

      tt <- unique(kres$term_name)[k]
      
      #if there is an intersection between the genes for the same pathway 
      if( length(intersect( str_split(kres[kres$term_name==tt,][1, ]$intersection, pattern = ',')[[1]],  str_split(kres[kres$term_name==tt,][2 ,]$intersection, pattern= ',')[[1]] ))>0 ){
        
                  for(i in 1:nrow(kres[kres$term_name==tt,])){
                    g_s[[i]] <-  genes_symbol[genes_symbol$ensembl_gene_id %in% str_split( kres[kres$term_name==tt,][i,]$intersection, pattern = ','  )[[1]] ,]$hgnc_symbol
                    gtp[tt, g_s[[i]]] <- kres[kres$term_name==tt,][i,]$query
                    u[i] <-  kres[kres$term_name==tt,][i,]$query
                  }        
        a <- intersect(g_s[[1]], g_s[[2]]) 
       gtp[tt ,a] <- paste0(u, collapse = '-')
       
       
      #if there is no intersection proceed as follows
      } else{
                
                for(i in 1:nrow(kres[kres$term_name==tt,])){
                        ges <-  genes_symbol[genes_symbol$ensembl_gene_id %in% str_split( kres[kres$term_name==tt,][i,]$intersection, pattern = ','  )[[1]] ,]$hgnc_symbol
                        gtp[tt ,ges] <- kres[kres$term_name==tt,][i,]$query
                    }
      }
      
}

ggtp <-  gtp[rownames(otp),rownames(t_t)[1:48]]

#---------------------------------------------------
colori<- ComplexHeatmap:::default_col(gtp[,rownames(t_t)[1:48]])
colori[c('f1', 'f2', 'f3', 'f1-f3')] <-  rep('white', 4)

a<-Heatmap(ggtp, 
           #column_title = "Genes and pathways KEGG",
           rect_gp = gpar(col = brewer.pal(5, 'Greys')[2], lwd = 0.25), 
           column_title_gp = gpar(fontsize = 20, fontface = "bold"),
           col=colori, 
           border=T,
           row_names_gp = gpar( fontsize=0),
           column_names_gp = gpar(fontsize=25),
           column_names_side = 'bottom',
           show_heatmap_legend = F,
           row_title = "Pathways",
           na_col = 'white', 
           height = unit(6.5, 'cm'), width = 28, 
           top_annotation =  columnAnnotation(legend = anno_empty(border = F, width =  unit(28, "cm"), height = unit(2, 'cm'))), #make space for the annotation
           
           
           cell_fun = function(j, i, x, y, width, height, fill) {
             if(!is.na(ggtp[,rownames(t_t)[1:48]][i, j] )){
               
               if(ggtp[,rownames(t_t)[1:48]][i, j]== 'f1' ){
                 grid.points(x = x, y = y, pch = 22, size = unit(18, "mm"), gp = gpar(col = brewer.pal(5, 'Paired')[2], fill=brewer.pal(5, 'Paired')[1] ))
               }
               
               
               if(ggtp[,rownames(t_t)[1:48]][i, j]== 'f2' ){
                 grid.points(x = x, y = y, pch = 22, size = unit(18, "mm"), gp = gpar(col = brewer.pal(5, 'Paired')[4], fill=brewer.pal(5, 'Paired')[3]))
               }
               
               if(ggtp[,rownames(t_t)[1:48]][i, j]== 'f3' ){
                 grid.points(x = x, y = y, pch = 22, size = unit(18, "mm"), gp = gpar(col = brewer.pal(6, 'Paired')[6], fill=brewer.pal(5, 'Paired')[5]))
               }
               
               if(ggtp[,rownames(t_t)[1:48]][i, j]== 'f1-f3' ){
                 grid.points(x = x, y = y, pch=22, size = unit(18, "mm"), gp = gpar(col = 'grey', fill='black'))
               }
             }
             
           }
)


b <- Heatmap(otp,  
             rect_gp = gpar(col = brewer.pal(5, 'Greys')[2], lwd = 0.25), 
             cluster_rows = F,
             col = c('white', 'white'), 
             border=T,
             show_column_dend = F,
             show_row_dend = F,
             row_names_gp = gpar(fontsize=3),
             row_names_side = 'left',
             width = unit(6.5, 'cm'),
             height = unit(28, 'cm'),
             show_heatmap_legend = F,
             column_names_gp = gpar(fontsize=25),
             column_names_side = 'bottom',
             column_names_rot = 360,
             top_annotation =  columnAnnotation(legend = anno_empty(border = F,
                                                  width =  unit(28, "cm"), height = unit(2, 'cm'))), #create an empty annotation for adding the legend elements 

             cell_fun = function(j, i, x, y, width, height, fill) {
               if(otp[i, j] > 0){
                 grid.circle(x = x, y = y, r = unit(otp[i,j], 'mm'), gp = gpar(fill = brewer.pal(7,'BrBG')[3], col = 'black')) #the unit of the -log10(pvalues) is mm as in the legend
                 #grid.text(signif((10^100)^(-otp[i,j]),2), x=x,  y = y, gp=gpar(fontsize=10)) #in case you want to see the pvalues in the centre of the circle
                 
               }
             }
      )

 
 

#save as pdf
pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure_pathways/semaphor_plot.pdf', height = 36, width = 36)
#plot the heatmaps
b+a
#for adding the legend on the pvalues scales, create the same heatmap with fixed values and copy paste them in illustrator 
for(k in 1:6) {
  i <- c(0.05,0.30,0.50, 0.70, 0.80,0.90)[k] #this are the positions 
  
  decorate_annotation('legend', {
    
    grid.circle(x=i, r= unit(-log10(c(10^-10,10^-8,10^-6,10^-4,10^-2,0.05)), 'mm')[k],gp = gpar(fill = brewer.pal(7,'BrBG')[3], col = 'black')) #this are the pvalues that are shown in the legend, the dimension of the -log10(p) is in mm
    grid.text(x=i, round(-log10(c(10^-10,10^-8,10^-6,10^-4,10^-2,0.05)), 3)[k]) #show the p
  })
}

dev.off()



#----- go circle plot-----------------------------------------------------------

go <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/go_pathway_analysis.RDS')


top_10 <- go$result[order(go$result$p_value, decreasing = F),][1:20,c('query', 'p_value', 'term_name')]

top_10$p_log <- -log10(top_10$p_value)

ggplot(top_10, aes(x = term_name, y=p_log, fill=query)) + geom_col(position = position_dodge2(width = 0.9, preserve = "single"), col='grey') + 
  scale_fill_manual(values=brewer.pal(5,'Paired')[c(1,3,5)]) +
  ylab('-log10(P-adjusted)') + xlab('GO term')+ ggtitle('Top 20 GO terms')+
  coord_flip()  +theme_classic() 


#search the top 10 unique pathways and searcj the pvalues for all the factors for the top 10 pathways
top_10_pathways <- unique(go$result[order(go$result$p_value, decreasing = F),][,'term_name'])[(1:10)]
kg <- go$result[go$result$term_name %in% top_10_pathways,c('query', 'p_value', 'term_name')]

mtp <- matrix(nrow = length(unique(kg$term_name)), ncol = 3)
colnames(mtp) <- c('f1', 'f2', 'f3')
rownames(mtp) <- unique(kg$term_name)

for(i in c('f1', 'f2', 'f3')){
  ko <- kg[ kg$query==i,]$term_name
  mtp[kg[ kg$query==i, ]$term_name , i] <- -log10(kg[ kg$query==i, ]$p_value) #importat to set the radius of the circle
}


for(i in c('f1', 'f2', 'f3')){
  mtp[,i][is.na(mtp[,i])] <- 0
}

colnames(mtp) <- toupper(colnames(mtp))

pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure_pathways/go_path_top10.pdf', height = 20, width = 10)
Heatmap(mtp,  
        rect_gp = gpar(col = brewer.pal(5, 'Greys')[2], lwd = 0.25), 
        column_title = 'Top 10 GO terms',
        cluster_rows = F,
        col = c('white', 'white'), 
        border=T,
        show_column_dend = F,
        show_row_dend = F,
        cluster_columns = F,
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
    
    grid.circle(x=i, r= unit(-log10(c(10^-15,10^-12,10^-9,10^-6,10^-3)), 'mm')[k],gp = gpar(fill = brewer.pal(7,'BrBG')[4], col = 'black')) #this are the pvalues that are shown in the legend, the dimension of the -log10(p) is in mm
    grid.text(x=i, round(-log10(c(10^-15,10^-12,10^-9,10^-6,10^-3)), 3)[k]) #show the p
  })
}

dev.off()


#------ react -----------

mart <- useDataset("hsapiens_gene_ensembl", useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl"))
loci.table <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/loci_table_names_nicola_genes.csv', data.table = F)

kegg <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/kegg_pathway_analysis.RDS')
react <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/react_pathway_analysis.RDS')



he <- react$result
genes <- unique(strsplit(paste0(he$intersection, collapse = ','), split = ',')[[1]])
genes_symbol <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values=  genes,mart=mart )

genes
t_t <- data.frame(genes) 
for(i in 1:length(he$term_name)){
  t_t[, paste0(paste(strsplit(he$term_name[i], split=' ')[[1]], collapse = '_'),'_' , he$query[i])] <- as.numeric(genes %in% strsplit(he[he$term_name==he$term_name[i] & he$query==he$query[i], 'intersection'], split=',' )[[1]])
  
}

t_t
#add gene symbol and remove ENSEMBLE
rownames(t_t) <- genes_symbol$hgnc_symbol[match(t_t$genes, genes_symbol$ensembl_gene_id)]   
t_t$genes <- NULL
t_t

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

ht_opt$ROW_ANNO_PADDING = unit(0.2, "cm")
#plot heatmap_gene_and_pathways.pdf             

Heatmap(tt, column_title = "Genes and pathways REACT genes",
        rect_gp = gpar(col = "white", lwd = 0.25, type='2'), 
        column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        left_annotation = rowAnnotation( '-log10(Padj)' = anno_barplot(   axis_param = list(direction = "reverse"),-log10(as.numeric(tt[,'p_value'])), width = unit(2, "cm"))),
        row_split = tt[,'trait'],
        column_split = sep,
        column_gap = unit(2, "mm"),
        border=T,
        row_gap = unit(5, "mm"), 
        row_names_gp = gpar(col =   RColorBrewer::brewer.pal(5, 'Paired')[c(1,3,5)], fontsize=8),
        column_labels = c(colnames(tt)[1:c(nrow(t_t)-2)], 'Factor', '') ,
        column_names_gp = gpar(fontsize=7),
        column_names_side = 'bottom',
        show_heatmap_legend = F,
        row_title = "Pathways",
        na_col = 'red', 
        col = colori, 
        heatmap_width = unit(28, 'cm'), heatmap_height = unit(12, 'cm')
)



