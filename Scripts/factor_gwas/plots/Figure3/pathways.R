library(data.table)
library(stringr)
library(EnsDb.Hsapiens.v75)
library(biomaRt)
library(ComplexHeatmap)
library(RColorBrewer)
library(dplyr)


kegg <- readRDS('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/10_genes_and_pathways/kegg_pathway_analysis.RDS')
mart <- useDataset("hsapiens_gene_ensembl", useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl"))

#----------- KEGG semaphor------------------------------------------------------

he <- kegg$result
genes <- unique(strsplit(paste0(he$intersection, collapse = ','), split = ',')[[1]])
genes_symbol <- getBM(filters= "ensembl_gene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol'),values=  genes,mart=mart )

#ecxclude the duplicated ones
genes_symbol <- genes_symbol[!duplicated(genes_symbol$ensembl_gene_id), ]


#create a table with both the pathways and the genes for each factor
t_t <- data.frame(genes) 
for(i in 1:length(he$term_name)){
  t_t[, paste0(paste(strsplit(he$term_name[i], split=' ')[[1]], collapse = '_'),'_' , he$query[i])] <- as.numeric(genes %in% strsplit(he[he$term_name==he$term_name[i] & he$query==he$query[i], 'intersection'], split=',' )[[1]])
  
}

#add gene symbol as row nameand remove ENSEMBLE
rownames(t_t) <- genes_symbol$hgnc_symbol[match(t_t$genes, genes_symbol$ensembl_gene_id)]   
t_t$genes <- NULL

#--------------------------------------
kg <- kegg$result[,c('query', 'p_value', 'term_name')]

#create mtp, a matrix containing the unique pathway name and the columns containing the factors
mtp <- matrix(nrow = length(unique(kg$term_name)), ncol = 3)
colnames(mtp) <- c('f1', 'f2', 'f3')
rownames(mtp) <- unique(kg$term_name)

#populate mtp with the z scores
for(i in c('f1', 'f2', 'f3')){
      kg[ kg$query==i,]$term_name
      mtp[, i][base::match(  kg[ kg$query==i, ]$term_name, names(mtp[, i]))] <- -log10(kg[ kg$query==i, ]$p_value) #importat to set the radius of the circle
}

# add zeros to the empty rows
for(i in c('f1', 'f2', 'f3')){
  mtp[,i][is.na(mtp[,i])] <- 0
}

#toupper case
colnames(mtp) <- toupper(colnames(mtp))

#order mtp with the shared pathways on top and the genes to be clustereb by factors
otp <- mtp[c(1,5,4,2,11,3,6,7,8,9,10,12,13),]



#create a table and add the information about the genes, shared or not
gtp <- matrix(nrow = length(unique(kg$term_name)), ncol = nrow(genes_symbol))
rownames(gtp) <- unique(kg$term_name)
colnames(gtp) <- genes_symbol$hgnc_symbol
kres <- kegg$result

#allocate vectors and lists
g_s <- list()
u <- vector()


#populate the matrix
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

ggtp <-  gtp[rownames(otp),rownames(t_t)]
#look at the combinations
names(table(ggtp))



#---------------------------------------------------
colori <- ComplexHeatmap:::default_col(gtp[,rownames(t_t)])
colori[c('f1', 'f2', 'f3', 'f1-f3', 'f2-f3')] <-  rep('white', 5)

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
             if(!is.na(ggtp[,rownames(t_t)][i, j] )){
               
               if(ggtp[,rownames(t_t)][i, j]== 'f1' ){
                 grid.points(x = x, y = y, pch = 22, size = unit(18, "mm"), gp = gpar(col = brewer.pal(6, 'Paired')[1], fill=brewer.pal(5, 'Paired')[2] ))
               }
               
               
               if(ggtp[,rownames(t_t)][i, j]== 'f2' ){
                 grid.points(x = x, y = y, pch = 22, size = unit(18, "mm"), gp = gpar(col = brewer.pal(5, 'Paired')[3], fill=brewer.pal(5, 'Paired')[4]))
               }
               
               if(ggtp[,rownames(t_t)][i, j]== 'f3' ){
                 grid.points(x = x, y = y, pch = 22, size = unit(18, "mm"), gp = gpar(col = brewer.pal(6, 'Paired')[5], fill=brewer.pal(6, 'Paired')[6]))
               }
               
               if(ggtp[,rownames(t_t)][i, j]== 'f1-f3' ){
                 grid.points(x = x, y = y, pch=22, size = unit(18, "mm"), gp = gpar(col = 'grey', fill='black'))
               }
               
               if(ggtp[,rownames(t_t)][i, j]== 'f2-f3' ){
                 grid.points(x = x, y = y, pch=22, size = unit(18, "mm"), gp = gpar(col = 'grey', fill='yellow'))
               }
             }
             
           }
)


b <- Heatmap(otp,  
             rect_gp = gpar(col = brewer.pal(5, 'Greys')[2], lwd = 0.25), 
             cluster_rows = F,
             cluster_columns = F,
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
pdf('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure3/semaphor_plot.pdf', height = 36, width = 36)
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
