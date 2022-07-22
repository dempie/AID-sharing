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


#load the dataset 

loci.table <- read.table('loci_definitions/final_locus_table.tsv', header = T)
#add the loci names
loci.table$final.locus=paste0(loci.table$Chr,":",loci.table$start,"-",loci.table$end,"_",loci.table$sub_locus)
fwrite(loci.table,'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/loci_table_names_nicola.csv', sep = ',')

#create a table of only factors
loci.table <- loci.table[loci.table$trait %in% c('f1', 'f2', 'f3'), ]
fwrite(loci.table,'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/loci_table_names_nicola_only_factors.csv', sep = ',')

#-----closest gene--------------------------------------------------------------

ref_genome <- EnsDb.Hsapiens.v75
ref_genes <- genes(ref_genome)
ref_genes <- ref_genes[ref_genes@elementMetadata$gene_biotype=='protein_coding',]
reg_genes_tss <- resize(ref_genes,width = 1 ,fix = 'start')

f_lead <- list()
G_list <- list()

mart <- useDataset("hsapiens_gene_ensembl", useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")) #select the database to convert and maake sure it is build 37 
loci.table$closest_gene <- rep('-', nrow(loci.table))

for(i in 1:3){
  tt <- c('f1', 'f2', 'f3')[i]
  f_lead[[tt]] <- GRanges(seqnames  =loci.table[loci.table$trait==tt, c('Chr')],   IRanges(names =loci.table[loci.table$trait==tt, c('SNP')] , start = loci.table[loci.table$trait==tt, 'bp']))
  elementMetadata(f_lead[[tt]])[['ensembl_gene_id']] <-  reg_genes_tss[nearest(f_lead[[tt]], reg_genes_tss),]@ranges@NAMES
  
  #add a column with the ensemble gene id into the factor loci table
  loci.table[loci.table$trait==tt, ][match( f_lead[[tt]]@ranges@NAMES, loci.table[loci.table$trait==tt, ]$SNP),]$closest_gene  <- unlist(elementMetadata(f_lead[[tt]])[['ensembl_gene_id']])
}


#------ make a list of genes for all the factors with ensemble_gened and gene symbol ----------
f_list_genes <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values= loci.table$closest_gene ,mart= mart)
loci.table$closest_gene_name <- rep(NA, nrow(loci.table))


for(i in 1:nrow(loci.table)){
  loci.table[i, 'closest_gene_name'] <- unique(f_list_genes[f_list_genes$ensembl_gene_id==loci.table[i, 'closest_gene'], ]$hgnc_symbol)
}
#add the ensemble gene id to the genes that do not have a hgnc symbol 
loci.table[loci.table$closest_gene_name=='', ]$closest_gene_name <- loci.table[loci.table$closest_gene_name=='', ]$closest_gene


#save the output
fwrite(loci.table, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/loci_table_names_nicola_genes.csv', sep=',')



#------ load the dataset and remove hla regione -------------------------------------------------------
loci.table <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/loci_table_names_nicola_genes.csv', data.table = F)

#remove hla locus
no_hla <- loci.table[! (loci.table$Chr==6 & data.table::between(loci.table$bp, 29000000, 34000000)),]
fwrite(no_hla, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/loci_table_names_nicola_genes_no_HLA.csv', sep = ',', col.names = T)


#----- make an upset plot-------------------------------------------------------
q <- make_comb_mat(list(f1=loci.table[loci.table$trait=='f1',]$final.locus, f2=loci.table[loci.table$trait=='f2',]$final.locus, f3=loci.table[loci.table$trait=='f3',]$final.locus))
pdf(width = 10, height = 5, file = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/upset_plot_final_loci_nicola.pdf')
UpSet(q, set_order = c("f1", "f2", "f3"), 
      comb_order = c(5,6,7,2,3,4,1),
      comb_col = c(brewer.pal(12, 'Paired')[c(12,12,12,12,1,3,5)]),
      top_annotation = upset_top_annotation(q, add_numbers = TRUE, height = unit(6, "cm")),
      right_annotation = upset_right_annotation(q, add_numbers = TRUE, width = unit(5,'cm') ),
      row_title = "Factor", 
      column_title = "Intersection of all loci")
dev.off()


#-------pathway analysis--------------------------------------------------------

#kegg analysis
kegg <- gost(query = list(f1=unique(loci.table[loci.table$trait=='f1',]$closest_gene), f2=unique(loci.table[loci.table$trait=='f2',]$closest_gene), f3=unique(loci.table[loci.table$trait=='f3',]$closest_gene)),
                           sources = c('KEGG'), significant = T, evcodes = T)
saveRDS(kegg, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/kegg_pathway_analysis.RDS')

#go analysis 
go <- gost(query = list(f1=unique(loci.table[loci.table$trait=='f1',]$closest_gene), f2=unique(loci.table[loci.table$trait=='f2',]$closest_gene), f3=unique(loci.table[loci.table$trait=='f3',]$closest_gene)),
           sources = c('GO'), significant = T, evcodes = T)
saveRDS(go, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/go_pathway_analysis.RDS')

react <- gost(query = list(f1=unique(loci.table[loci.table$trait=='f1',]$closest_gene), f2=unique(loci.table[loci.table$trait=='f2',]$closest_gene), f3=unique(loci.table[loci.table$trait=='f3',]$closest_gene)),
           sources = c('REAC'), significant = T, evcodes = T)
saveRDS(react, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/react_pathway_analysis.RDS')


#------ make a csv of genes for all the factors with ensemble_gened and gene symbol ----------

#---------------------list of genes for plotting the pathways-------------------
kegg <-  readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/kegg_pathway_analysis.RDS')
keggs <- kegg$result[,c('term_name', 'intersection', 'query', 'term_id')]
out <- list()
ent <- list()
for(i in 1:length(unique( keggs$term_name[duplicated(keggs$term_name)]))){
  pp<- keggs[duplicated(keggs$term_name),]$term_name[i]
  kid <- keggs[duplicated(keggs$term_name),]$term_id[i]
  for(k in 1:nrow( kegg$result[kegg$result$term_name==pp,])){
    out[[pp]][[kegg$result[kegg$result$term_name==pp,][k, ]$query]] <- kegg$result[kegg$result$term_name==pp,][k, ]$intersection
    ent[[paste0(pp, kid,collapse = ':') ]][[paste0(kegg$result[kegg$result$term_name==pp,][k, ]$query, kid, collapse = ':')]]<- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values=  kegg$result[kegg$result$term_name==pp,][k, ]$intersection,mart= mart)
  }
  
}


saveRDS(ent        , 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/kegg_genes_in_shared_pathway_analysis.RDS')


for(i in 1:6){
  tt<- names(ent)[i]
  
  for( k in 1:2){
    ff <- names(ent[[tt]])[k]
    ent[[tt]][[ff]]$trait <- rep(ff, nrow(ent[[tt]][[ff]]))
    ent[[tt]][[ff]]$pathway <- rep(tt, nrow(ent[[tt]][[ff]]))
  }
  
}


tocall <- list()

for(i in 1:6){
  tt <- names(ent)[i]
  
  for( k in 1:2){
    ff <- names(ent[[tt]])[k]
    
    tocall[[paste0(tt, ff, collapse = '-')]] <- ent[[tt]][[ff]]
  }
  
}

export <- do.call(rbind,tocall)

write.table(export, row.names = T, quote = F, sep=',', file='outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/genes_in_shared_pathways_table.csv')


#-------------------table of pathways---------------------------------

#kegg table
kg <- kegg$result[, c('query','p_value', 'term_name', 'intersection')] 

#use gene symbols
he <- kegg$result
genes <- unique(strsplit(paste0(he$intersection, collapse = ','), split = ',')[[1]])
genes_symbol <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values=  genes,mart=mart )

kg$hgnc <- rep('NA', nrow(kg))
for(i in 1:nrow(kg)){
kg$hgnc[i] <- paste0(genes_symbol[ genes_symbol$ensembl_gene_id %in% strsplit( kg$intersection[i], split = ',')[[1]] ,]$hgnc_symbol, collapse = ',  ')
}

colnames(kg) <- c('Trait', 'P-value (adjusted)', 'Pathway', 'ENSEMBL', 'Gene symbols')
kg$Trait <- toupper(kg$Trait)
plott <- kg[,c(3,1,2, 5)]



formattable(plott,align =c("l","c","c" ,"r") ,
            list( 'Pathway' = formatter( "span", style = style(color = "grey",font.weight = "bold")))
)


#go table

gg <-  go$result[, c('query','p_value', 'term_name', 'intersection')] 
gg <- head(gg[order(gg$p_value), ],20)
genes <- unique(strsplit(paste0(gg$intersection, collapse = ','), split = ',')[[1]])
genes_symbol <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values=  genes,mart=mart )
gg$hgnc <- rep('NA', nrow(gg))

for(i in 1:nrow(gg)){
  gg$hgnc[i] <- paste0(genes_symbol[ genes_symbol$ensembl_gene_id %in% strsplit( gg$intersection[i], split = ',')[[1]] ,]$hgnc_symbol, collapse = ',  ')
}

colnames(gg) <- c('Trait', 'P-value (adjusted)', 'GO term', 'ENSEMBL', 'Gene symbols')
gg$Trait <- toupper(gg$Trait)
plott <- gg[,c(3,1,2, 5)]
formattable(plott,align =c("l","c","c" ,"r") ,
            list( 'GO term' = formatter( "span", style = style(color = "grey",font.weight = "bold")))
)


#----plot the heatmap of the pathways-------------------------------------------

he <- kegg$result
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

ht_opt$ROW_ANNO_PADDING = unit(0.2, "cm")
#plot heatmap_gene_and_pathways.pdf             
pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/heatmap_gene_and_pathways_only_path_genes.pdf', height = 8, width = 16)
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


t_t <- data.frame(genes) 
for(i in 1:length(he$term_name)){
  t_t[, paste0(paste(strsplit(he$term_name[i], split=' ')[[1]], collapse = '_'),'_' , he$query[i])] <- as.numeric(genes %in% strsplit(he[he$term_name==he$term_name[i] & he$query==he$query[i], 'intersection'], split=',' )[[1]])
  
}

genes_symbol[!str_detect(genes_symbol$hgnc_symbol,'' ),]$hgnc_symbol <- genes_symbol[!str_detect(genes_symbol$hgnc_symbol,'' ),]$ensembl_gene_id
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

ht_opt$ROW_ANNO_PADDING = unit(0.2, "cm")
#plot heatmap_gene_and_pathways.pdf         


#plot heatmap_gene_and_pathways_with_all_genes.pdf             
pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/heatmap_gene_and_pathways_all_genes_full_matrix.pdf', height = 8, width = 16)
Heatmap(tt, column_title = "Genes and pathways all genes KEGG",
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
        column_names_gp = gpar(fontsize=2),
        column_names_side = 'bottom',
        show_heatmap_legend = F,
        row_title = "Pathways",
        na_col = 'red', 
        col = colori, 
        heatmap_width = unit(30, 'cm'), heatmap_height = unit(10, 'cm')
)
dev.off()



#-------------------heatmpap and go terms------------------------------------



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



#----plot the GO terms heatmap -------------------------------------------------

pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/heatmap_GO_gene_and_pathways_only_path_genes.pdf', height = 20, width = 20)
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


head(go$result[order(go$result$p_value), c('query', 'p_value', 'term_name')],10)

go_f <- list()
for(i in 1:3){
  tt <- c('f1', 'f2', 'f3')[i]
go_f[[tt]] <- go$result[ go$result$query==tt,]
}


#---- vem diagram go terms------------------------------------------------------
library(VennDiagram)
venn.diagram(
  x =list(f1=go_f$f1$term_name, f2=go_f$f2$term_name, f3=go_f$f3$term_name),
  category.names = c("F1" , "F2 " , "F3"),
  filename = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/GO_venn_diagramm.png',
  output=F,
  
  # Output features
  imagetype="png" ,
  height = 1200 , 
  width = 1200 , 
  resolution = 300,
  main= 'GO terms shared among Factors',

  # Circles
  lwd = 2,
  lty = 'blank',
  fill =brewer.pal(5, 'Paired')[c(1,3,5)] ,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)



#-------p value distribution of shared verus unique ----------------------------

Reduce(intersect, list(go_f$f1$term_name, f2=go_f$f2$term_name, f3=go_f$f3$term_name))
Reduce(setdiff, list(go_f$f1$term_name, f2=go_f$f2$term_name, f3=go_f$f3$term_name))
a<-go$result[go$result$term_name %in% Reduce(intersect, list(go_f$f1$term_name, f2=go_f$f2$term_name, f3=go_f$f3$term_name)),]$p_value 
b <- go$result[go$result$term_name %in% Reduce(setdiff, list(go_f$f1$term_name, f2=go_f$f2$term_name, f3=go_f$f3$term_name)),]$p_value 

  f1_u <- go$result[go$result$term_name %in% Reduce(setdiff, list(go_f$f1$term_name, f2=go_f$f2$term_name, f3=go_f$f3$term_name)),]$p_value 
  f2_u <- go$result[go$result$term_name %in% Reduce(setdiff, list(go_f$f2$term_name, f2=go_f$f1$term_name, f3=go_f$f3$term_name)),]$p_value 
  f3_u <- go$result[go$result$term_name %in% Reduce(setdiff, list(go_f$f3$term_name, f2=go_f$f2$term_name, f1=go_f$f1$term_name)),]$p_value 
  

  
  
  a1 <- data.frame(a=a, info=rep('shared', length(a)))
  b1 <- data.frame(a=c(f1_u, f2_u, f3_u), info=rep('unique', length(c(f1_u, f2_u, f3_u))))
  ee <- rbind(a1, b1)


library(ggplot2)

ggplot(ee) + geom_density(aes( x=a, group=info,col=info)) +xlab('p_value') 
ggplot(ee, aes( y=-log10(a), x=info,col=info)) + geom_boxplot() +  geom_jitter() +scale_y_log10() +ylab('-log10(pvalue)') 


#----- go barplot---------------------------------------------------------------

go <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/go_pathway_analysis.RDS')


top_10 <- go$result[order(go$result$p_value, decreasing = F),][1:20,c('query', 'p_value', 'term_name')]

top_10$p_log <- -log10(top_10$p_value)

ggplot(top_10, aes(x = term_name, y=p_log, fill=query)) + geom_col(position = position_dodge2(width = 0.9, preserve = "single"), col='grey') + 
  scale_fill_manual(values=brewer.pal(5,'Paired')[c(1,3,5)]) +
  ylab('-log10(P-adjusted)') + xlab('GO term')+ ggtitle('Top 20 GO terms')+
  coord_flip()  +theme_classic() 




#---circle map

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

#order mtp with the shared pathways on top and the genes to be clustered by factors
#otp <- mtp[c(1,3,5,8,9,4,2,6,7,10,11,12,13,14,15),]

mtp

colnames(mtp) <- toupper(colnames(mtp))

Heatmap(mtp,  
        rect_gp = gpar(col = brewer.pal(5, 'Greys')[2], lwd = 0.25), 
        column_title = 'Top 10 GO terms',
        cluster_rows = F,
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










