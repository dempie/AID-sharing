library(pathview)
library(RColorBrewer)
library(stringr)


ent <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways_conditional/kegg_genes_in_shared_pathway_analysis.RDS')
setwd('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure_pathways/kegg_color_pathview/')
colfactor <- RColorBrewer::brewer.pal(5, 'Paired')[c(1,3,5)]
names(colfactor) <- c('f1', 'f2', 'f3')
for(i in 1:length(ent)){
            fs_inv <- list()
            tt <- ent[i]
            
            path <- strsplit(names(tt), split = ':')[[1]][2] 
            
            for(k in 1:length(tt[[1]])){
                fs_inv[[k]] <- names(tt[[1]][k])
            }
            
              a <- rep(1, nrow(tt[[1]][[fs_inv[[1]]]]))
              names(a) <- tt[[1]][[fs_inv[[1]]]]$entrezgene_id
              
              b <- rep(-1, nrow(tt[[1]][[fs_inv[[2]]]]))
              names(b) <- tt[[1]][[fs_inv[[2]]]]$entrezgene_id
              
          
          if(length(intersect(names(a),names(b)))!=0){
                    q <- rep(0, length(intersect(names(a),names(b))))
                    names(q) <- intersect(names(a),names(b))
                    a <- a[!names(a) %in% names(q)]
                    b <- b[!names(b) %in% names(q)]
                    
                    pv.out <- pathview(gene.data =  c(a,b,q) , pathway.id = path,
                                       species = "hsa", sample.layer=F, out.suffix=paste0(path,'_', paste(unlist(fs_inv), collapse = '-'), collapse = '_' ), 
                                       low = list(gene = colfactor[fs_inv[[2]]], cpd = "blue"), 
                                       mid = list(gene = "yellow", cpd = "blue"), 
                                       high = list(gene = colfactor[fs_inv[[1]]], cpd = "yellow"))
          } else {
      
                    pv.out <- pathview(gene.data =  c(a,b) , pathway.id = path ,
                                   species = "hsa", sample.layer=F, out.suffix=paste0(path,'_', paste(unlist(fs_inv), collapse = '-'), collapse = '_' ),
                                   low = list(gene = colfactor[fs_inv[[2]]], cpd = "blue"), 
                                   mid = list(gene = "yellow", cpd = "blue"), 
                                   high = list(gene = colfactor[fs_inv[[1]]], cpd = "yellow"))
              }
          
}



