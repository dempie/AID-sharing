library(formattable)
library(biomaRt)
library(data.table)


mart <- useDataset("hsapiens_gene_ensembl", useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")) #select the database to convert and maake sure it is build 37 

kegg <-  readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/kegg_pathway_analysis.RDS')

#-------------------table of pathways---------------------------------

#--------- kegg table complete -------------------------------------------------
kg <- kegg$result[, c('query','p_value', 'term_name', 'intersection')] 

#use gene symbols
he <- kegg$result
genes <- unique(strsplit(paste0(he$intersection, collapse = ','), split = ',')[[1]])
genes_symbol <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values=  genes,mart=mart )

#add the gene symbols to each row (i.e. each pathway)
kg$hgnc <- rep('NA', nrow(kg))
for(i in 1:nrow(kg)){
  kg$hgnc[i] <- paste0(genes_symbol[ genes_symbol$ensembl_gene_id %in% strsplit( kg$intersection[i], split = ',')[[1]] ,]$hgnc_symbol, collapse = ',  ')
}

#format for exporting with the rigth names
colnames(kg) <- c('Trait', 'P-value (adjusted)', 'Pathway', 'ENSEMBL', 'Gene symbols')
kg$Trait <- toupper(kg$Trait)
kegg.table <- kg[,c(3,1,2, 5)]


formattable(kegg.table,align =c("l","c","c" ,"r") ,
            list( 'Pathway' = formatter( "span", style = style(color = "grey",font.weight = "bold")))
)

#save the csv
write.table(kegg.table, sep = ',', file = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/supplementary/kegg.table.csv')



#--------- go table complete ---------------------------------------------------

#go table complete, (all the results)
go <-  readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/go_pathway_analysis.RDS')
gg <-  go$result[, c('query','p_value', 'term_name', 'intersection')] 
#gg <- head(gg[order(gg$p_value), ],20)
genes <- unique(strsplit(paste0(gg$intersection, collapse = ','), split = ',')[[1]])
genes_symbol <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values=  genes,mart=mart )
gg$hgnc <- rep('NA', nrow(gg))

for(i in 1:nrow(gg)){
  gg$hgnc[i] <- paste0(genes_symbol[ genes_symbol$ensembl_gene_id %in% strsplit( gg$intersection[i], split = ',')[[1]] ,]$hgnc_symbol, collapse = ',  ')
}
#add the names
colnames(gg) <- c('Trait', 'P-value (adjusted)', 'GO term', 'ENSEMBL', 'Gene symbols')
gg$Trait <- toupper(gg$Trait)
go.table <- gg[,c(3,1,2, 5)]

#save the table ascsv
write.table(go.table, sep = ',', file = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/supplementary/go.table.csv')

formattable(plott,align =c("l","c","c" ,"r") ,
            list( 'GO term' = formatter( "span", style = style(color = "grey",font.weight = "bold")))
)


#--------- go table top 20  ---------------------------------------------

#add the names
colnames(gg) <- c('Trait', 'P-value (adjusted)', 'GO term', 'ENSEMBL', 'Gene symbols')
gg$Trait <- toupper(gg$Trait)
go.table_top_20 <- gg[1:20,c(3,1,2, 5)]

#save the table ascsv
write.table(go.table_top_20, sep = ',', file = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/supplementary/go.table_top20.csv')


formattable(go.table_top_20,align =c("l","c","c" ,"r") ,
            list( 'GO term' = formatter( "span", style = style(color = "grey",font.weight = "bold")))
)




#--------- reactome table top --------------------------------------------------
 react <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/react_pathway_analysis.RDS')


react <-  react$result[, c('query','p_value', 'term_name', 'intersection')] 
genes <- unique(strsplit(paste0(react$intersection, collapse = ','), split = ',')[[1]])
genes_symbol <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values=  genes,mart=mart )
react$hgnc <- rep('NA', nrow(react))

for(i in 1:nrow(react)){
  react$hgnc[i] <- paste0(genes_symbol[ genes_symbol$ensembl_gene_id %in% strsplit( react$intersection[i], split = ',')[[1]] ,]$hgnc_symbol, collapse = ',  ')
}
#add the names
colnames(react) <- c('Trait', 'P-value (adjusted)', 'Reactome pathway', 'ENSEMBL', 'Gene symbols')
react$Trait <- toupper(react$Trait)
react.table <- react[,c(3,1,2, 5)]

#save the table ascsv
write.table(react.table, sep = ',', file = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/supplementary/react.table.csv')

formattable(react.table,align =c("l","c","c" ,"r") ,
            list( 'Reactome' = formatter( "span", style = style(color = "grey",font.weight = "bold")))
)












