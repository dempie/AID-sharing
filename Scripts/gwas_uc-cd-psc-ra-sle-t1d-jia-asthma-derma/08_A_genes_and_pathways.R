
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
library(gtools)
library(RColorBrewer)
library(circlize)
#--------- first find the nearest genes to each lead SNP------------------------


#I used ensmble names, on build 37
factor_loci <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/hyprcoloc_only_factors/factor_loci.txt', data.table = F) 
factor_loci <- factor_loci[factor_loci$trait %in% c('f1', 'f2', 'f3'),]

ref_genome <- EnsDb.Hsapiens.v75
ref_genes <- genes(ref_genome)
ref_genes <- ref_genes[ref_genes@elementMetadata$gene_biotype=='protein_coding',]
reg_genes_tss <- resize(ref_genes,width = 1 ,fix = 'start')

f_lead <- list()
G_list <- list()

mart <- useDataset("hsapiens_gene_ensembl", useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")) #select the database to convert and maake sure it is build 37 
factor_loci$closest_gene <- rep('-', nrow(factor_loci))

for(i in 1:3){
  tt <- c('f1', 'f2', 'f3')[i]
  f_lead[[tt]] <- GRanges(seqnames  =factor_loci[factor_loci$trait==tt, c('chr')],   IRanges(names =factor_loci[factor_loci$trait==tt, c('SNP')] , start = factor_loci[factor_loci$trait==tt, 'BP']))
  elementMetadata(f_lead[[tt]])[['ensembl_gene_id']] <-  reg_genes_tss[nearest(f_lead[[tt]], reg_genes_tss),]@ranges@NAMES
  
  #add a column with the ensemble gene id into the factor loci table
  factor_loci[factor_loci$trait==tt, ][match(factor_loci[factor_loci$trait==tt, ]$SNP, f_lead[[tt]]@ranges@NAMES),]$closest_gene  <- unlist(elementMetadata(f_lead[[tt]])[['ensembl_gene_id']])
}


fwrite(factor_loci, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/factor_loci_moloc_closest_gene_info.txt') 


#make a list of genes for all the factors with ensemble_gened and gene symbol
f_list_genes <- list()
for(i in 1:3){
  tt <- c('f1', 'f2', 'f3')[i]
  f_list_genes[[tt]] <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values= factor_loci[factor_loci$trait==tt, ]$closest_gene ,mart= mart)
}




#plot the upset plot of the genes
q <- make_comb_mat(list(f1=f_list_genes$f1$ensembl_gene_id, f2=f_list_genes$f2$ensembl_gene_id, f3=f_list_genes$f3$ensembl_gene_id))
UpSet(q)
pdf(width = 10, height = 5, file = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/nearest_genes_upset_plot_factors.pdf')
UpSet(q, set_order = c("f1", "f2", "f3"), 
      comb_order = order(comb_size(q), decreasing = T),
      comb_col = c(rcartocolor::carto_pal(12, 'Safe')[c(1,2,3)][comb_degree(q)]),
      top_annotation = upset_top_annotation(q, add_numbers = TRUE, height = unit(6, "cm")),
      right_annotation = upset_right_annotation(q, add_numbers = TRUE, width = unit(5,'cm') ),
      row_title = "Factor", 
      column_title = "Intersection of all genes ")
dev.off()




#--- gprofiler KEGG ------------------------------------------------------------
factor_loci<- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/factor_loci_moloc_closest_gene_info.txt', data.table=F)
gost_test_region_2 <- gost(query = list(f1=factor_loci[factor_loci$trait=='f1',]$closest_gene, f2=factor_loci[factor_loci$trait=='f2',]$closest_gene, f3=factor_loci[factor_loci$trait=='f3',]$closest_gene),
                          sources = c('KEGG'), significant = T, evcodes = T)
saveRDS(gost_test_region_2, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/gprofiler_kegg_only.RDS')



#-----gprofiler GO-------------------------------------------------------------
factor_loci<- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/factor_loci_moloc_closest_gene_info.txt', data.table=F)
gost_test_GO <- gost(query = list(f1=factor_loci[factor_loci$trait=='f1',]$closest_gene, f2=factor_loci[factor_loci$trait=='f2',]$closest_gene, f3=factor_loci[factor_loci$trait=='f3',]$closest_gene),
                     sources = c('GO'), significant = T, evcodes = T)
saveRDS(gost_test_GO, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/goprofiler_GO_only.RDS')



#---------------------list of genes for plotting the pathways-------------------
kegg <-  readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/gprofiler_kegg_only.RDS')
keggs <- kegg$result[,c('term_name', 'intersection', 'query')]
out <- list()
for(i in 1:3){
      pp<- keggs[duplicated(keggs$term_name),]$term_name[i]
            for(k in 1:nrow( kegg$result[kegg$result$term_name==pp,])){
          out[[pp]][[kegg$result[kegg$result$term_name==pp,][k, ]$query]] <- kegg$result[kegg$result$term_name==pp,][k, ]$intersection
      }
          
}

saveRDS(out, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/genes_in_shared_pathways.RDS')

#---- heatmap -----------------------------------------------------------------

gost_all <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/gprofiler_kegg_only.RDS')
he <- gost_all$result
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
names(sep) <- colnames(tt)
                     

#plot heatmap_gene_and_pathways.pdf             
pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/heatmap_gene_and_pathways_all_genes.pdf', height = 8, width = 16)
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
dev.off()

saveRDS(tt, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/heatmap_gene_and_pathways_all_genes_table.RDS')



#---------- heatmap to show all genes -----------------------------------------


he_all_genes <- gost_all$result
genes <- unique(do.call(rbind, f_list_genes)$ensembl_gene_id)
genes_symbol <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values=  genes,mart=mart )
#add ENSID to the ones that do not have gene id
genes_symbol[which(genes_symbol$hgnc_symbol==''),]$hgnc_symbol <- genes_symbol[which(genes_symbol$hgnc_symbol==''),]$ensembl_gene_id 

genes
t_t <- data.frame(genes) 
for(i in 1:length(he$term_name)){
  t_t[, paste0(paste(strsplit(he_all_genes$term_name[i], split=' ')[[1]], collapse = '_'),'_' , he$query[i])] <- as.numeric(genes %in% strsplit(he_all_genes[he_all_genes$term_name==he_all_genes$term_name[i] & he_all_genes$query==he_all_genes$query[i], 'intersection'], split=',' )[[1]])
  
}

t_t
#add gene symbol and remove ENSEMBLE
rownames(t_t) <- genes_symbol$hgnc_symbol[match(t_t$genes, genes_symbol$ensembl_gene_id)]   
t_t$genes <- NULL
t_t

t_t['trait',] <- c(he_all_genes$query) 
t_t['p_value',] <- c(he_all_genes$p_value)
tt <- t(t_t)


#color selection
library(rcartocolor)
display_carto_all(colorblind_friendly = TRUE)

colori<- ComplexHeatmap:::default_col(x = tt)
colori[c('f1', 'f2', 'f3')] <-  rcartocolor::carto_pal(12, 'Safe')[c(4,5,10)]#colorblind safe
colori[c('0','1')] <- c('#a6cee3','#1f78b4')
c('#a6cee3','#1f78b4','#b2df8a')
#column split
sep <- c(rep('a', 179), rep('b', 2))
names(sep)<- colnames(tt)


#plot heatmap_gene_and_pathways_with_all_genes.pdf             
pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/heatmap_gene_and_pathways_all_genes_full_matrix.pdf', height = 8, width = 16)
Heatmap(tt, column_title = "Genes and pathways all genes",
        rect_gp = gpar(col = "white", lwd = 0.25), 
        column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        left_annotation = rowAnnotation( '-log10(Padj)' = anno_barplot(   axis_param = list(direction = "reverse"),-log10(as.numeric(tt[,'p_value'])), width = unit(2, "cm"))),
        row_split = tt[,'trait'],
        column_split = sep,
        column_gap = unit(3, "mm"),
        row_gap = unit(5, "mm"), 
        row_names_gp = gpar(col =  rcartocolor::carto_pal(12, 'Safe')[c(4,5,10)], fontsize=10),
        column_labels = c(colnames(tt)[1:179], 'Factor', '') ,
       column_names_gp = gpar(fontsize=3),
        column_names_side = 'bottom',
        show_heatmap_legend = F,
        row_title = "Pathways",
        na_col = 'black', 
        col = colori, 
        heatmap_width = unit(30, 'cm'), heatmap_height = unit(15, 'cm')
)
dev.off()



#----clusterprofiler go terms -------------------------------------------------

# ggo <- enrichGO(gene     = as.character(f3_genes$entrezgene_id),
#                OrgDb    = org.Hs.eg.db,
#                ont      = "BP")
# 
# ggo <- pairwise_termsim(ggo)
# ggo@result
# ggo <- simplify(ggo, cutoff=0.7, by="p.adjust", select_fun=min)
# emapplot(ggo)
# 
# 
# kk <- list()
# plot_kk <- list()
# for(i in 1:3){
#   tt <- c('f1','f2', 'f3')[i]
# kk[[tt]]<- clusterProfiler::enrichKEGG(f_list_genes[[tt]]$entrezgene_id)
# plot_kk[[tt]]<- cnetplot(kk[[tt]], showCategory = 5, circular=T, categorySize="pvalue") + ggtitle(paste0(tt,' KEGG pathways and genes'))
# }
# 
# dotplot(ggo, 
#         showCategory = 12000, 
#         title = "Enriched Pathways",
#         font.size = 8)
# 
# ?dotplot
# 
# plot_kk
# gost_test_region_2$result
# 
# #dotplot 
# dotplot(
#   ggo,
#   x = "Cluster",
#   color = "p.adjust",
#   showCategory = 5,
#   split = NULL,
#   font.size = 12,
#   title = "",
#   size = NULL,
#   label_format = 30,
# )
# 
# #dotplot KEGG
# cc <- compareCluster(list(f1=f_list_genes$f1$entrezgene_id, f2=f_list_genes$f2$entrezgene_id , f3=f_list_genes$f3$entrezgene_id), fun = "enrichKEGG")
# pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/cluster_profiler/dotplot_kegg_cluster_profiler_all_factors.pdf', height = 14, width = 14)
# dotplot(
#   cc,
#   x = "Cluster",
#   color = "p.adjust",
#   showCategory = 50,
#   font.size = 5,
#   title = "",
# #  by = "geneRatio",
#   size = NULL,
#   includeAll = F,
#   label_format = 5
# )
# dev.off()
# 
# #cnetplot KEGG
# pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/cluster_profiler/network_kegg_cluster_profiler_all_factors.pdf', height = 14, width = 14)
# cnetplot(cc, showCategory = 10, colorEdge=T, circular=F, node_label='all')
# dev.off()




#----try using semantic similarity measurement----------------------------------

# gprofiler_all<- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/g')
# 
# go_f <- list()
# for(i in 1:3){
#       tt <- c('f1', 'f2', 'f3')[i]
#       go_f[[tt]] <- c(gprofiler_all$result[gprofiler_all$result$query==tt,]$term_id)
#       go_f[[tt]] <- go_f[[tt]][str_detect(go_f[[tt]], 'GO:')]
# }
# 
# 
# 
# hsGO <- godata('org.Hs.eg.db', ont="BP")
# mgoSim(go_f$f1, go_f$f2, semData=hsGO, measure="Wang", combine='BMA')
# goSim(as.list(go_f$f1), semData=hsGO, measure="Wang")
# semsim <- mclusterSim(list(f1=f_list_genes$f1$entrezgene_id, f2=f_list_genes$f2$entrezgene_id , f3=f_list_genes$f3$entrezgene_id), semData=hsGO, measure="Wang", combine="BMA")
# corrplot(semsim, is.corr = T, addCoef.col = 'black', type = 'upper') #the information content shows that the genes are highly correlated in their semantics, althiugh they are different 
# 




#-----unique genes replication -------------------------------------------------

#replicate the same analysis by taking out the shared genes among all the factors
a <- combinations(n=3, r=2, c('f1', 'f2', 'f3'), repeats.allowed = F) #find the possible combinations
shared_genes <- list()

for(i in 1:3){
    shared_genes[[paste0(a[i,1], '-',a[i,2])]] <- intersect(f_list_genes[[a[i,1]]]$ensembl_gene_id, f_list_genes[[a[i,2]]]$ensembl_gene_id) 
}


sh_genes <- unique(do.call(c, shared_genes))

f_unique_genes <- list()

for(i in 1:3){
    #create a list of genes that are unique and not shared
    tt <- c('f1', 'f2', 'f3')[i]
    f_unique_genes[[tt]] <- f_list_genes[[tt]][which(!f_list_genes[[tt]]$ensembl_gene_id %in% sh_genes),]
}

f_unique_genes
gost_unique <- gost(  query = list(f1=f_unique_genes$f1$ensembl_gene_id, f2=f_unique_genes$f2$ensembl_gene_id, f3=f_unique_genes$f3$ensembl_gene_id),
                           sources = c( 'KEGG'), significant = T, evcodes = T)
gost_all <- gost(  query = list(f1=f_list_genes$f1$ensembl_gene_id, f2=f_list_genes$f2$ensembl_gene_id, f3=f_list_genes$f3$ensembl_gene_id),
                   sources = c( 'KEGG'), significant = T, evcodes = T) 
  
  
gost_all$result
gost_unique$result



#the heatmap 
he_u <- gost_unique$result[, c('query', 'term_name', 'p_value', "intersection") ] 
genes <- unique(strsplit(paste0(he_u$intersection, collapse = ','), split = ',')[[1]])
genes_symbol <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values=  genes,mart= mart)

genes
t_t <- data.frame(genes) 
for(i in 1:length(he_u$term_name)){
  t_t[, paste0(paste(strsplit(he_u$term_name[i], split=' ')[[1]], collapse = '_'),'_' , he_u$query[i])] <- as.numeric(genes %in% strsplit(he_u[he_u$term_name==he_u$term_name[i] & he_u$query==he_u$query[i], 'intersection'], split=',' )[[1]])
  
}

t_t
#add gene symbol and remove ENSEMBLE
rownames(t_t) <- genes_symbol$hgnc_symbol[match(t_t$genes, genes_symbol$ensembl_gene_id)]   
t_t$genes <- NULL
t_t

t_t['trait',] <- c(he_u$query) 
t_t['p_value',] <- c(he_u$p_value)
tt <- t(t_t)


#color selection
library(rcartocolor)
display_carto_all(colorblind_friendly = TRUE)

colori<- Com
colori[c('f1', 'f2', 'f3')] <-  rcartocolor::carto_pal(12, 'Safe')[c(4,5,10)]#colorblind safe
colori[c('0','1')] <- c('#a6cee3','#1f78b4')
c('#a6cee3','#1f78b4','#b2df8a')
#column split
sep <- c(rep('a', 35), rep('b', 2))
names(sep)<- colnames(tt)
dim(tt)

#plot heatmap_gene_and_pathways.pdf             
pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/heatmap_gene_and_pathways_factor_unique_genes.pdf', height = 8, width = 16)
Heatmap(tt, column_title = "Genes and pathways Factor Unique Genes",
        rect_gp = gpar(col = "white", lwd = 0.25), 
        column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        left_annotation = rowAnnotation( '-log10(p_adjusted)' = anno_barplot(   axis_param = list(direction = "reverse"),-log10(as.numeric(tt[,'p_value'])), width = unit(2, "cm"))),
        row_split = tt[,'trait'],
        column_split = sep,
        column_gap = unit(3, "mm"),
        row_gap = unit(5, "mm"), 
        row_names_gp = gpar(col =  rcartocolor::carto_pal(12, 'Safe')[c(4,5,10)], fontsize=10),
        column_labels = c(colnames(tt)[1:35], 'Factor', '') ,
        column_names_side = 'bottom',
        show_heatmap_legend = F,
        row_title = "Pathways",
        na_col = 'black', 
        col = colori, 
        heatmap_width = unit(28, 'cm'), heatmap_height = unit(15, 'cm')
)
dev.off()





