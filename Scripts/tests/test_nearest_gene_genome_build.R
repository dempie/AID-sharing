library(data.table)
library(ChIPpeakAnno)
library(stringr)
library(gprofiler2)
library(GenomicRanges)
library(EnsDb.Hsapiens.v75)
library(biomaRt)
library(ComplexHeatmap)
library(clusterProfiler)
library(enrichplot)
library(GOSemSim)
library(dplyr)
library(gtools)

#--------- first find the nearest genes to each lead SNP------------------------


#I used ensmble names, on build 37
factor_loci <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/moloc_factors/factor_loci_moloc_info.txt', data.table = F) 
# 
# ref_genome <-TxDb.Hsapiens.UCSC.hg19.knownGene
ref_genome <- EnsDb.Hsapiens.v75
ref_genes <- genes(ref_genome)
ref_genes <- ref_genes[ref_genes@elementMetadata$gene_biotype=='protein_coding',]
ref_genes@elementMetadata

ranges(ref_genes)
ranges(reg_genes_tss)
start(ref_genes)
reg_genes_tss <- resize(ref_genes,width = 1 ,fix = 'start')


# ref_genes <- genes(ref_genome, single.strand.genes.only=FALSE )
# ref_genes<- c(unlist(ref_genes))
f_lead <- list()
G_list <- list()

mart <- useDataset("hsapiens_gene_ensembl", useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")) #select the database to convert and maake sure it is build 37 
factor_loci$closest_gene <- rep('-', nrow(factor_loci))

for(i in 1:3){
  tt <- c('f1', 'f2', 'f3')[i]
  f_lead[[tt]] <- GRanges(seqnames  =factor_loci[factor_loci$trait==tt, c('chr')],   IRanges(names =factor_loci[factor_loci$trait==tt, c('SNP')] , start = factor_loci[factor_loci$trait==tt, 'BP']))
  elementMetadata(f_lead[[tt]])[['ensembl_gene_id']] <-  ref_genes[nearest(f_lead[[tt]], reg_genes_tss, select=c('arbitrary')),]@ranges@NAMES
  
  #add a column with the ensemble gene id into the factor loci table
  factor_loci[factor_loci$trait==tt, ][match(factor_loci[factor_loci$trait==tt, ]$SNP, f_lead[[tt]]@ranges@NAMES),]$closest_gene  <- unlist(elementMetadata(f_lead[[tt]])[['ensembl_gene_id']])
}


fwrite(factor_loci, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/factor_loci_moloc_closest_gene_info.txt') 
factor_loci <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/factor_loci_moloc_closest_gene_info.txt', data.table = F) 

#make a list of genes for all the factors with ensemble_gened and gene symbol
f_list_genes <- list()
for(i in 1:3){
  tt <- c('f1', 'f2', 'f3')[i]
  f_list_genes[[tt]] <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values= factor_loci[factor_loci$trait==tt, ]$closest_gene ,mart= mart)
}





factor_loci

facto_notss <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/factor_loci_moloc_closest_gene_info.txt', data.table = F) 


facto_notss$closest_gene %in% factor_loci$closest_gene



gost_tss <- gost(query = list(f1=f_list_genes$f1$ensembl_gene_id, f2=f_list_genes$f2$ensembl_gene_id, f3=f_list_genes$f3$ensembl_gene_id),
                           sources = c('KEGG'), significant = T, evcodes = T)
gost_no_tss <-  gost(query = list(f1=facto_notss[facto_notss$trait=='f1', 'closest_gene'], f2=facto_notss[facto_notss$trait=='f2', 'closest_gene'] , f3=facto_notss[facto_notss$trait=='f3', 'closest_gene']),
                     sources = c('KEGG'), significant = T, evcodes = T)

gost_no_tss$result[,colnames(gost_no_tss$result)[1:12]]
gost_tss$result[,colnames(gost_tss$result)[1:12]]

#ààààà


he <-gost_tss$result
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
library(rcartocolor)
display_carto_all(colorblind_friendly = TRUE)

colori<- ComplexHeatmap:::default_col(x = tt)
colori[c('f1', 'f2', 'f3')] <-  rcartocolor::carto_pal(12, 'Safe')[c(4,5,10)]#colorblind safe
colori[c('0','1')] <- c('#a6cee3','#1f78b4')
c('#a6cee3','#1f78b4','#b2df8a')
#column split
sep <- c(rep('a', 37), rep('b', 2))
names(sep)<- colnames(tt)

Heatmap(tt, column_title = "Genes and pathways all genes",
        rect_gp = gpar(col = "white", lwd = 0.25), 
        column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        left_annotation = rowAnnotation( '-log10(Padj)' = anno_barplot(   axis_param = list(direction = "reverse"),-log10(as.numeric(tt[,'p_value'])), width = unit(2, "cm"))),
        row_split = tt[,'trait'],
        column_split = sep,
        column_gap = unit(3, "mm"),
        row_gap = unit(5, "mm"), 
        row_names_gp = gpar(col =  rcartocolor::carto_pal(12, 'Safe')[c(4,5,10)], fontsize=10),
        column_labels = c(colnames(tt)[1:37], 'Factor', '') ,
        column_names_side = 'bottom',
        show_heatmap_legend = F,
        row_title = "Pathways",
        na_col = 'black', 
        col = colori, 
        heatmap_width = unit(28, 'cm'), heatmap_height = unit(15, 'cm')
)




factor_loci[69,]












