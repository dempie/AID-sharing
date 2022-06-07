#in this script plotting and analysis of colocalization will be done

library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(tidyr)
library(ChIPpeakAnno)
library(stringr)
library(RColorBrewer)
library(moloc)
library(gprofiler2)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene )
library(biomaRt)


#function to create a genomic ranges object useful to find overlapps between loci 

create_range <- function(res, chr='CHR', start='start', end='end' ) {
  #check if there are overlapping ranges
  cnd<- length(GRanges( seqnames =  res$chr ,IRanges(names = rownames(res) ,start = as.numeric(res$start), end = as.numeric(res$end)))@ranges) != (length(reduce(GRanges( seqnames =  res$chr ,IRanges(names = rownames(res) ,start = as.numeric(res$start), end = as.numeric(res$end))))))
  if(cnd==T) {warning('There are overlaps inside this thing') }

  return(GRanges( seqnames =  res$chr ,IRanges(names = rownames(res) ,start = as.numeric(res$start), end = as.numeric(res$end))) )
 
}

#----------analysis of the factor gwas------------------------------------------

coloc_factors <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/hyprcoloc_only_factors/colocalize_factors.RDS')
loci_factors <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/hyprcoloc_only_factors/factor_loci.txt', data.table=F)
f1_range <- create_range(loci_factors[loci_factors$trait=='f1',])
f2_range <- create_range(loci_factors[loci_factors$trait=='f2',])
f3_range <- create_range(loci_factors[loci_factors$trait=='f3',])


#find the loci that are physically overlapping among all the three factors 
ovl_f <- findOverlapsOfPeaks(f1_range, f2_range, f3_range , connectedPeaks = 'keepAll')

#a loop to create as many list as the number of traits that have at least one unique locus
#it takes the output from findoverlapps function and the names of the traits as in findoverlapps!!
lists_forupset <- function(ovl_f, traits){ 
    #allocate lists
    require(ComplexHeatmap)
    un <- list()
    final <- list(list())
    sh <- list()
    #create the lists and put it the eleemtnts that are unique
    for( i in 1:length(traits)){
          tt<-traits[i]
          filt <- str_starts(names(ovl_f$uniquePeaks@ranges),tt, negate = FALSE)
          un[[tt]] <- names(ovl_f$uniquePeaks@ranges)[filt]
    }
    
    #put into the lists the element that are shared among the lists
    for(k in 1:length(ovl_f$mergedPeaks$peakNames)){
    
    chek_in <- traits %in% substring(ovl_f$mergedPeaks$peakNames[[k]], first = 1, last = 2) #check if there is f1, f2 or f3
    
            for(u in 1:length(chek_in)){
              tt<-traits[u]
              if(chek_in[u]==T){
                un[[tt]][[length(un[[tt]])+1]]<- paste0(ovl_f$mergedPeaks$peakNames[[k]], collapse = '-')
              }
              
            }
    }
    
    #-generate a matrix showing which locus is present in each disease
    output <- list()
    for(k in 1:length(un)){
          tt <- c(names(un))[k]
          to_get <- strsplit(unlist(strsplit(un[[tt]], split = '-')), split='__')
          
          for(i in c(1:length(to_get))){
            output[[tt]][[i]]<-to_get[[i]][2]
            
          }
          output[[tt]]<- unlist(output[[tt]])
      
    }
    output <- list_to_matrix(output)
    output <- output[order(as.numeric(rownames(output))),]
    
    #generate the output
    return(list(loci=output, traits_overlap=un))
}



up_n <- lists_forupset(ovl_f = ovl_f, traits = c('f1', 'f2', 'f3')) 


cm <- make_comb_mat(up_n$traits_overlap, mode = 'distinct')
pdf(width = 10, height = 5, file = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/upset_plot_factors.pdf')
UpSet(cm, set_order = c("f1", "f2", "f3"), comb_order = order(comb_size(cm), decreasing = T),
      comb_col = c((brewer.pal(6, 'Set3'))[4:6][comb_degree(cm)]),
       top_annotation = upset_top_annotation(cm, add_numbers = TRUE, height = unit(6, "cm")),
      right_annotation = upset_right_annotation(cm, add_numbers = TRUE, width = unit(5,'cm') )
       )
dev.off()


#----- function to run moloc-----------------------------------------------------
#https://github.com/clagiamba/moloc

molocalize <- function(list_of_paths, trait_names, loci, N_hat){
      #packages 
      # require(data.table)
      # require(dplyr)
      # require(GenomicRanges)
      # require(moloc)
      
      SNPs <- list()
      list_of_files <- list()
      
      for( i in c(1:length(trait_names))){
                #load the sumstats and put them into a list
                list_of_files[[i]] <- fread(list_of_paths[[i]], data.table = F)
                SNPs[[i]]<- list_of_files[[i]]$SNP
              }
              
      names(list_of_files) <- unlist(trait_names)
              #start by searching the commong SNPs between all the gwas and generate a referenc set 
      shared_SNPs <- Reduce(intersect, SNPs) 
      cat(paste0('The number of shared SNPs is ', length(shared_SNPs )))
              
      reference_file <- fread('SNP/reference.1000G.maf.0.005.txt.gz', data.table = F) #load the reference file
      ref_set <- reference_file[reference_file$SNP   %in%  shared_SNPs, ] #create a reference set of the positions of the shared SNPs
      loc_index <- sort(unique(loci$pan_locus))
              
      betas <- list()
      SE <- list()
      res_coloc <- list()
      report <- data.frame(row.names = paste0('locus_', loc_index )  ) #allocate space to produce a report of the number of SNPs per locus
      lo <- list() 
      
      for(i in c(1:length(unique(loci$pan_locus)))){
                
                    #select the active locus, the beginning and the end of the locus
                    start_loc <- min(loci[loci$pan_locus==loc_index[i], 'start'] )
                    end_loc <- max(loci[loci$pan_locus==loc_index[i], 'end'] )
                    chr_loc <- unique(loci[loci$pan_locus==loc_index[i], 'chr'])
                    
                    #the active SNP for locus i
                    SNP_active <- ref_set[c( between(ref_set$BP,  start_loc ,  end_loc)  & ref_set$CHR== chr_loc ), ]$SNP
                    
                    #if there are less than 10 SNPs raise the window of SNPs to take of +-100 kb 
                    if(length(SNP_active)<10) {
                      start_loc_2 <-  as.numeric(ifelse((as.numeric(start_loc) - 100000)<0, 0, as.numeric(start_loc) - 100000 )) #not below zero 
                      end_loc_2 <- as.numeric(as.numeric(end_loc) + 100000 )
                      #select the active SNPs
                      SNP_active <- ref_set[c( between(ref_set$BP,  start_loc_2 ,  end_loc_2)  & ref_set$CHR== chr_loc ), ]$SNP
                    } else {
                      #select the active SNPs
                      SNP_active <- ref_set[c( between(ref_set$BP,  start_loc ,  end_loc)  & ref_set$CHR== chr_loc ), ]$SNP
                    }
                    
                    #allocate the space
                    the_locus <- data.frame(row.names = c(1:length(SNP_active)))
                    to_moloc <- list()
                    #for each GWAS take out the BETAs and the SE
                    for(k in c(1:length(unique(unlist( loci[loci$pan_locus==i,]$trait))))){
                              #take only the ones that are significant fot that loci
                              tt<- unique(unlist((loci[loci$pan_locus==i,]$trait)))[k]   #important the unique(), as  in case of multiple locus from the same trait, the last iteration will overwright the previous, no importance as it will be identical 
                              tr_SNP_active <- list_of_files[[tt]][ list_of_files[[tt]]$SNP %in% SNP_active ,]
                              tr_SNP_active <- tr_SNP_active[!duplicated(tr_SNP_active ), ] #remove duplicated SNPs 
                              #exctract the betas and se
                              the_locus[, 'SNP'] <- SNP_active
                              the_locus[, 'BETA'] <- select(tr_SNP_active, 'BETA')
                              the_locus[, 'SE'] <- select(tr_SNP_active, 'SE')
                              the_locus[, 'MAF'] <- ref_set[ref_set$SNP %in% tr_SNP_active$SNP, ]$MAF
                              the_locus[, 'N'] <- rep(N_hat[[tt]], length(SNP_active))
                              to_moloc[[tt]] <- the_locus #important in case of multiple locus from the same trait, the last iteration will overwright the previous, no importance as it will be identical 
                            
                            }
                    
                    #run coloc only if the number of traits for that locus is more than one (otherwise you will have an error) and if there are at least 10 SNPs
                    if(length(to_moloc)>1 & length(SNP_active)>10){
                      res_coloc[[paste0('locus_', i)]] <- moloc_test(to_moloc)
                      
                    } 
                    
                    lo[[paste0('locus_', i)]] <- to_moloc
                    report[i,'SNP_active'] <-  length(SNP_active)
                    report[i,'traits'] <- paste0(loci[loci$pan_locus==i,]$trait, collapse = ',')
                    print(i)
                  }
      
      
      output <- list(lo=lo, res_coloc=res_coloc, report=report) 
      return(output)
}

#--------------- run moloc on the factors --------------------------------------
##Calculate Effective Sample Size for Factors
#restrict to MAF of 40% and 10%
N_hat_F <- list()
for(k in 1:3){
f <- fread(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/factor', k ,'_gwas_final_withQindex.txt'), data.table = F)
subsetf<-subset(f, f$MAF <= .4 & f$MAF >= .1)
subsetf <- subsetf[!is.na(subsetf$SE),]

N_hat_F[[paste0('f', k)]]<-mean(1/((2*subsetf$MAF*(1-subsetf$MAF))*subsetf$SE^2))
}


#gather te info about the paths and the gwas names
sstats_names <- list.files('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/')
sstats_names  <- sstats_names[c(-8, -9, -12)]
the_paths<- as.list(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/', sstats_names))

#loop for loading all the files, create a vector of names and find the SNPs that are present in all the GWAS
gwas_names <- list()


for( i in c(1:length(sstats_names ))){
  #createa lit of trait names
  gwas_names[i] <- strsplit(sstats_names , '_')[[i]][[1]]
}
#-------------------------------------------------------------

#run molocalize
factor_loci <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/hyprcoloc_only_factors/factor_loci.txt', data.table = F)
factor_coloc <- molocalize(list_of_paths = the_paths[c(4,5,6)], trait_names = gwas_names[c(4,5,6)],loci = factor_loci, N_hat = N_hat_F)
saveRDS(factor_coloc, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/moloc_factors/factor_moloc_output.RDS')
factor_coloc <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/moloc_factors/factor_moloc_output.RDS')



factor_coloc$report #look into locus 119 and 95 as one single locus from one factor is overlapping with two from the other factor 

locus_119_a <- loci_factors[loci_factors$pan_locus==119, ][c(1,3),]
locus_119_a$pan_locus <- 1
moloc_119 <-  molocalize(list_of_paths = the_paths[c(4,5,6)], trait_names = gwas_names[c(4,5,6)],loci = locus_119_a , N_hat = N_hat_F) #posterior 9.973760e-01


locus_119_b <- loci_factors[loci_factors$pan_locus==119, ][c(2,3),]
locus_119_b$pan_locus <- 1
moloc_119_b <-  molocalize(list_of_paths = the_paths[c(4,5,6)], trait_names = gwas_names[c(4,5,6)],loci = locus_119_b , N_hat = N_hat_F) # 9.973760e-01
#the problem is that f2 range contains both f1 and f2, so no way of distinguishing it 
#-------------------------------------------------------------------------------

#script for adding the results of coloc to the loci table 
#locmoloc column indicates if the locus is co localizing or not, if it is indicates the traits involved in the co localization

factor_loci$configuration <- rep('-', nrow(factor_loci))
factor_loci$posterior_configuration <- rep('-', nrow(factor_loci))
factor_loci$locmoloc <- rep(0, nrow(factor_loci))

for(k in 1:length(unique(factor_loci$pan_locus))){
  if(nrow(factor_loci[factor_loci$pan_locus==k, ])==1){
    
    factor_loci[factor_loci$pan_locus==k, ]$locmoloc <- factor_loci[factor_loci$pan_locus==k, ]$trait
    factor_loci[factor_loci$pan_locus==k, ]$configuration <- 'no_over'
    
  } else{
    act <- factor_coloc$res_coloc[[paste0("locus_", k)]]$priors_lkl_ppa
    m<- which.max(c(act$PPA))
    factor_loci[factor_loci$pan_locus==k, ]$configuration <- paste0(rownames(act)[m])
    factor_loci[factor_loci$pan_locus==k, ]$posterior_configuration <- act$PPA[m]
    
          if( paste0(rownames(act)[m]) %in% c( "abc"  ,  "ab" )){
                  factor_loci[factor_loci$pan_locus==k, ]$locmoloc<- paste0(unique(factor_loci[factor_loci$pan_locus==k, ]$trait),collapse = '-' ) #if they colocalize
                  
          } else if (paste0(rownames(act)[m]) %in% c( "a,b"  ,  "a,b,c" ) ){
                  for(r in 1:nrow(factor_loci[factor_loci$pan_locus==k, ])){
                          factor_loci[factor_loci$pan_locus==k, ][r,]$locmoloc<- paste0(unique(factor_loci[factor_loci$pan_locus==k, ][r,]$trait),collapse = '-' )
                  }
            
          } else if (paste0(rownames(act)[m]) ==  "ab,c" ) {
                  factor_loci[factor_loci$pan_locus==k & factor_loci$trait %in% c('f1', 'f2'),  ]$locmoloc <-   paste0(unique(factor_loci[factor_loci$pan_locus==k & factor_loci$trait %in% c('f1', 'f2'),  ]$trait),collapse = '-' ) 
                  factor_loci[factor_loci$pan_locus==k & factor_loci$trait %in% c('f3'),  ]$locmoloc <-   paste0(unique(factor_loci[factor_loci$pan_locus==k & factor_loci$trait %in% c('f3'),  ]$trait),collapse = '-' ) 
                
          } else if (paste0(rownames(act)[m]) ==  "a,bc" ) {
                  factor_loci[factor_loci$pan_locus==k & factor_loci$trait %in% c('f2', 'f3'),  ]$locmoloc <-   paste0(unique(factor_loci[factor_loci$pan_locus==k & factor_loci$trait %in% c('f2', 'f3'),  ]$trait),collapse = '-' )
                  factor_loci[factor_loci$pan_locus==k & factor_loci$trait %in% c('f1'),  ]$locmoloc <-   paste0(unique(factor_loci[factor_loci$pan_locus==k & factor_loci$trait %in% c('f1'),  ]$trait),collapse = '-' )
          }
  }
  
}


#save the output 
fwrite(factor_loci[order(factor_loci$pan_locus),], 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/moloc_factors/factor_loci_moloc_info.txt')
factor_loci <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/moloc_factors/factor_loci_moloc_info.txt', data.table = F) 




#plot an upset plot of the loci after colocalization---------------------------- 
traits <- c('f1', 'f2', 'f3')
ups <- list()
for(q in 1:length(traits)){
  tt <- traits[q]
  ups[[tt]] <- paste0(factor_loci[str_detect(factor_loci$locmoloc, tt) & factor_loci$trait==tt , ]$pan_locus_name, '_', factor_loci[str_detect(factor_loci$locmoloc, tt) & factor_loci$trait==tt , ]$locmoloc)
  
}

a <- make_comb_mat(ups, mode = 'distinct')
pdf(width = 10, height = 5, file = 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/colocalization_upset_plot_factors.pdf')
UpSet(a, set_order = c("f1", "f2", "f3"), comb_order = order(comb_size(a), decreasing = T),
      comb_col = c((brewer.pal(6, 'Set3'))[4:6][comb_degree(a)]),
      top_annotation = upset_top_annotation(a, add_numbers = TRUE, height = unit(6, "cm")),
      right_annotation = upset_right_annotation(a, add_numbers = TRUE, width = unit(5,'cm') )
)
dev.off()


#------------------ gprofiler --------------------------------------------------


#--------- first find the nearest genes to each lead SNP------------------------
#I used ensmble names, on build 37
factor_loci <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/moloc_factors/factor_loci_moloc_info.txt', data.table = F) 

ref_genome <-TxDb.Hsapiens.UCSC.hg19.knownGene 
ref_genes <- genes(ref_genome)

f_lead <- list()
G_list <- list()

mart <- useDataset("hsapiens_gene_ensembl", useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")) #select the database to convert and maake sure it is build 37 
factor_loci$closest_gene <- rep('-', nrow(factor_loci))

for(i in 1:3){
        tt <- c('f1', 'f2', 'f3')[i]
        f_lead[[tt]] <- GRanges(seqnames  =paste0('chr',factor_loci[factor_loci$trait==tt, c('chr')]),   IRanges(names =factor_loci[factor_loci$trait==tt, c('SNP')] , start = factor_loci[factor_loci$trait==tt, 'BP']))
        elementMetadata(f_lead[[tt]])[['entrezgene_id']] <-  ref_genes[nearest(f_lead[[tt]], ref_genes, ignore.strand=T),]@elementMetadata@listData[[1]]
        #check number of nearest genes equal to number of SNPs
        ifelse(isTRUE(length(f_lead[[tt]]@seqnames) == length(elementMetadata(f_lead[[tt]])[['entrezgene_id']] )), print(paste0('All good with nearest genes numbers in ', tt)), warning('not equal number or ranges and nearest genes in ', tt) )
      
        G_list[[tt]] <- getBM(filters= "entrezgene_id", attributes= c("entrezgene_id","ensembl_gene_id", 'hgnc_symbol' ),values=  elementMetadata(f_lead[[tt]])[['entrezgene_id']],mart= mart)
        elementMetadata(f_lead[[tt]])[['ensembl_gene_id']] <- G_list[[tt]][base::match(f_lead[[tt]]@elementMetadata[[1]],  G_list[[tt]]$entrezgene_id),][,2]
        
        #add a column with the ensemble gene id into the factor loci table
        factor_loci[factor_loci$trait==tt, ][match(factor_loci[factor_loci$trait==tt, ]$SNP, f_lead[[tt]]@ranges@NAMES),]$closest_gene  <- unlist(elementMetadata(f_lead[[tt]])[['ensembl_gene_id']])
}


fwrite(factor_loci, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/moloc_factors/factor_loci_moloc_closest_gene_info.txt') 



#gprofiler----------------------------------------------------------------------
gost_test_region <- gost( organism = , query = list(f1=unique(G_list$f1$ensembl_gene_id), f2=unique(G_list$f2$ensembl_gene_id), f3=unique(G_list$f3$ensembl_gene_id)),
                          sources = c("GO:BP", "REAC", 'KEGG'), significant = T)
gostplot(gost_test_region, capped = F)

#-------------------------------------------------------------------------------
#write a matrix to be uploaded intro stringr to build the matrix
fwrite(matrix(factor_loci$closest_gene, ncol = 1, nrow = length(factor_loci$closest_gene)), 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/test.string.csv')

#read the file generated by string with protein ids and make it readable by biomaRt
string <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/string_protein_annotations.csv', data.table = F)
string <- string[,c(1,2)]
string[,2] <- str_sub(string$identifier, start = 6)
string

#find ensembl_gene_id to the proteins from string
s <- getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id","ensembl_gene_id"),values=  string[,2] ,mart= mart)
string$ensembl_gene_id <- s[base::match(string$identifier, s$ensembl_peptide_id),]$ensembl_gene_id

#chech that the enemble ids are the ones present in origin in my dataset
mean(string$ensembl_gene_id %in% factor_loci$closest_gene) # 1 all of them


string$trait <- rep(NA, nrow(string))
ttr <- list()
factor_loci_noNA<- factor_loci[!is.na(factor_loci$closest_gene),]
string_no_NA <- string[!is.na(string$ensembl_gene_id),]
string$ensembl_gene_id[is.na(string$ensembl_gene_id)] <- 0

for(i in 1:length(string_no_NA$ensembl_gene_id)) {
  a<- string_no_NA$ensembl_gene_id[i]
   string_no_NA$trait[i]<- paste0(factor_loci_noNA[factor_loci_noNA$closest_gene==a, ]$trait, collapse = '-')
  
}
string_no_NA

string$trait <- string_no_NA[match(string$ensembl_gene_id, string_no_NA$ensembl_gene_id),]$trait

fwrite(string, 'outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/string_mapping_columns.csv')
