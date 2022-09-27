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

loci.table <- fread('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/09_loci_definitions_Nicola/loci_definitions/final_locus_table.tsv', data.table = F)

loci.table <- loci.table[loci.table$trait %in% c('f1', 'f2', 'f3'),]

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
fwrite(loci.table, 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/10_genes_and_pathways/loci_table_names_nicola_genes.csv', sep=',')





#-------pathway analysis--------------------------------------------------------

#kegg analysis
kegg <- gost(query = list(f1=unique(loci.table[loci.table$trait=='f1',]$closest_gene), f2=unique(loci.table[loci.table$trait=='f2',]$closest_gene), f3=unique(loci.table[loci.table$trait=='f3',]$closest_gene)),
             sources = c('KEGG'), significant = T, evcodes = T)
saveRDS(kegg, 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/10_genes_and_pathways/kegg_pathway_analysis.RDS')

#go analysis 
go <- gost(query = list(f1=unique(loci.table[loci.table$trait=='f1',]$closest_gene), f2=unique(loci.table[loci.table$trait=='f2',]$closest_gene), f3=unique(loci.table[loci.table$trait=='f3',]$closest_gene)),
           sources = c('GO'), significant = T, evcodes = T)
saveRDS(go, 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/10_genes_and_pathways/go_pathway_analysis.RDS')

#react
react <- gost(query = list(f1=unique(loci.table[loci.table$trait=='f1',]$closest_gene), f2=unique(loci.table[loci.table$trait=='f2',]$closest_gene), f3=unique(loci.table[loci.table$trait=='f3',]$closest_gene)),
              sources = c('REAC'), significant = T, evcodes = T)
saveRDS(react, 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/10_genes_and_pathways/react_pathway_analysis.RDS')



#--------------------plotting the heatmaps--------------------------------------

#----------------------------------------------------------------------------
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
pdf('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/10_genes_and_pathways/heatmap_gene_and_pathways_only_path_genes.pdf', height = 8, width = 16)
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


#---------------------------- Semaphor plot-------------------------------------




















