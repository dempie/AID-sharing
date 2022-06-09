
library(data.table)
library(qqman)
library(dplyr)
library(ggpubr)
library(Gviz)
library(EnsDb.Hsapiens.v75)
library(biomaRt)

#a function to plot specific loci, provide the files and the position of the locus
plot_locus <- function(list_of_files, trait_names, start, end, chr, y_a=0, y_b=20, locus_name='-'){
      #packages 
    # require(data.table)
    #require(qqman)
     names(list_of_files) <- unlist(trait_names)
    #select the active locus, the beginning and the end of the locus
    start_loc <- start
    end_loc <- end
    chr_loc <- chr
    #allocate the space
    to_plot <- list()
    #for each GWAS take out the BETAs and the SE
    for(k in 1:length(list_of_files)){

          tt<- trait_names[[k]]   
          tr_SNP_active <- list_of_files[[tt]]
          tr_SNP_active <- tr_SNP_active[c(tr_SNP_active$CHR==chr & between(tr_SNP_active$BP, left = start, right = end)),]
          tr_SNP_active <- tr_SNP_active[!duplicated(tr_SNP_active ), ] #remove duplicated SNPs 
          #exctract pvalue
          to_plot[[tt]] <- data.frame('SNP'=tr_SNP_active$SNP, 'P'=as.numeric(tr_SNP_active$P), 'CHR'=as.numeric(tr_SNP_active$CHR), 'BP'=as.numeric(tr_SNP_active$BP) )
    }
    #check the shared SNPs among the files
    tr_SNP <- list()
    for(k in 1:length(list_of_files)){
          tt<- trait_names[[k]]  
          tr_SNP[[tt]] <- to_plot[[tt]]$SNP
    }
    shared_snp <- Reduce(intersect, tr_SNP)
    
    #take only the shared snp to plot
    to_p <- list()
    for(k in 1:length(list_of_files)){
          tt<- trait_names[[k]] 
          to_p[[tt]]<-to_plot[[tt]][to_plot[[tt]]$SNP %in% shared_snp,]   
    }
    
    par(mfrow=c(length(list_of_files),1))
    par(mar=c(5,5,3,3))
    for(n in 1:length(to_plot)){
                tp <-trait_names[[n]]
                manhattan(to_p[[tp]], chr="CHR", bp="BP", snp="SNP", p="P" ,ylim=c(y_a,y_b),xlim=c(start, end), main=paste0(tp,'  chr', chr,'_', start,'_', end, '  locus name ', locus_name ), suggestiveline =F , 
                          col = c("brown2"), highlight = c(to_p[[tp]][which.min(to_p[[tp]]$P), ]$SNP) )  
                grid(nx = NULL, ny = NULL,
                     lty = 2,      # Grid line type
                     col = "gray", # Grid line color
                     lwd = 1)
                
    }
   
}


#----- plot the locus ----------------------------------------------------------

#gather te info about the paths and the gwas names
sstats_names <- list.files('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/')
sstats_names  <- sstats_names[c(-8, -9, -12)]
the_paths <- as.list(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/', sstats_names))
the_paths <- the_paths[c(4,5,6)]
#loop for loading all the files, create a vector of names and find the SNPs that are present in all the GWAS
gwas_names <- list()
for( i in c(1:length(sstats_names ))){
  #createa lit of trait names
  gwas_names[i] <- strsplit(sstats_names , '_')[[i]][[1]]
}


list_of_files <- list()

for( i in c(1:3)){
  #load the sumstats and put them into a list
  list_of_files[[i]] <- fread(the_paths[[i]], data.table = F)
}




plot_locus(list_of_files = list_of_files, trait_names =gwas_names[c(4,5,6)], start =11020255, end = 11400000, chr = 16, y_b= 25, locus_name = 205)


library(Gviz)
library(biomaRt)

ref_genes <- genes(EnsDb.Hsapiens.v75)
ref_genes <- ref_genes[ref_genes@elementMetadata$gene_biotype=='protein_coding',]


a <- GRanges(seqnames = 16 ,IRanges( start =11020255, end = 11400000))


#function for plotting
granges2df <- function(x) {
  df <- as(x, "data.frame")
  df <- df[,c("seqnames","start","end","strand","group_name", 'exon_id')]
  colnames(df)[1] <- "chromosome"
  colnames(df)[5] <- "transcript"
  df
}

txdf <- select(EnsDb.Hsapiens.v75,
               keys=keys(EnsDb.Hsapiens.v75, "GENEID"),
               columns=c("GENEID","TXID", 'SYMBOL'),
               keytype="GENEID")

ebt <- exonsBy(EnsDb.Hsapiens.v75, by="tx")

d <- list()
for(i in 1:length(a@elementMetadata$gene_id)){

idx <- txdf$GENEID ==  a@elementMetadata$gene_id[i]
txs <- txdf$TXID[idx]
#all the xons for these transcripts
ebt2 <- ebt[txs]
df <- granges2df(ebt2)
df$gene <- a@elementMetadata$gene_id[i]
df$symbol <-  txdf[match(df$gene, txdf$GENEID), ]$SYMBOL 
d[[i]] <- df
}


d <- do.call(rbind, d) 
grt <- GeneRegionTrack(d)


















d <- list()
for(i in 1:length(a@elementMetadata$gene_id)){
idx <- txdf$GENEID == a@elementMetadata$gene_id[i]
exons <- txdf$TXID[idx]
exons <- ebt[exons]

df <- granges2df(exons)
df$gene <- a@elementMetadata$gene_id[i]
d[[i]] <- df
}
d <- do.call(rbind, d) 
grt <- GeneRegionTrack(df)
gax <- GenomeAxisTrack()
plotTracks(list(gax,grt))

df







#a function to plot specific loci, provide the files and the position of the locus
plot_locus <- function(list_of_files, trait_names, start, end, chr, y_a=0, y_b=20, locus_name='-'){
  #packages 
  # require(data.table)
  #require(qqman)
  names(list_of_files) <- unlist(trait_names)
  #select the active locus, the beginning and the end of the locus
  start_loc <- start
  end_loc <- end
  chr_loc <- chr
  #allocate the space
  to_plot <- list()
  #for each GWAS take out the BETAs and the SE
  for(k in 1:length(list_of_files)){
    
    tt<- trait_names[[k]]   
    tr_SNP_active <- list_of_files[[tt]]
    tr_SNP_active <- tr_SNP_active[c(tr_SNP_active$CHR==chr & between(tr_SNP_active$BP, left = start, right = end)),]
    tr_SNP_active <- tr_SNP_active[!duplicated(tr_SNP_active ), ] #remove duplicated SNPs 
    #exctract pvalue
    to_plot[[tt]] <- data.frame('SNP'=tr_SNP_active$SNP, 'P'=as.numeric(tr_SNP_active$P), 'CHR'=as.numeric(tr_SNP_active$CHR), 'BP'=as.numeric(tr_SNP_active$BP) )
  }
  #check the shared SNPs among the files
  tr_SNP <- list()
  for(k in 1:length(list_of_files)){
    tt<- trait_names[[k]]  
    tr_SNP[[tt]] <- to_plot[[tt]]$SNP
  }
  shared_snp <- Reduce(intersect, tr_SNP)
  
  #take only the shared snp to plot
  list_of_granges <- list()
  
  for(k in 1:length(list_of_files)){
    tt <- trait_names[[k]] 
    to_p <- to_plot[[tt]][to_plot[[tt]]$SNP %in% shared_snp,]   
    gr_tplot <- GRanges(seqnames = to_p$CHR, IRanges(start = to_p$BP, end = to_p$BP, names = to_p$SNP) )
    elementMetadata(gr_tplot)[['p']] <- -log(to_p[match(gr_tplot@ranges@NAMES, to_p$SNP), ]$P)
    list_of_granges[[k]] <- gr_tplot
  }
  
  dTrack <- list()
  for(n in 1:length(to_plot)){
    tt <- trait_names[[n]]
    dTrack[[n]] <- DataTrack(list_of_granges[[n]], name = tt, start = start, end = end, chromosome = chr, genome = 'hg19', col='red', ylim=c(y_a, y_b)
                             ,grid=T, frame=T, baseline=0, col.baseline= "darkgreen", lty.baseline= 1, lwd.baseline=2, col.axis= "darkgreen",
                             lty.grid= 2, lwd.grid=0.3, col.title="darkgreen" 
                             )
    
  }
  
  #retrieve the info on the genes
  ref_genes <- genes(EnsDb.Hsapiens.v75)
  ref_genes <- ref_genes[ref_genes@elementMetadata$gene_biotype=='protein_coding',]
  
  
  #function for plotting
  granges2df <- function(x) {
    df <- as(x, "data.frame")
    df <- df[,c("seqnames","start","end","strand","group_name", 'exon_id')]
    colnames(df)[1] <- "chromosome"
    colnames(df)[5] <- "transcript"
    df
  }
  
  txdf <- select(EnsDb.Hsapiens.v75,
                 keys=keys(EnsDb.Hsapiens.v75, "GENEID"),
                 columns=c("GENEID","TXID", 'SYMBOL'),
                 keytype="GENEID")
  
  ebt <- exonsBy(EnsDb.Hsapiens.v75, by="tx")
  
  d <- list()
  for(i in 1:length(a@elementMetadata$gene_id)){
    
    idx <- txdf$GENEID ==  a@elementMetadata$gene_id[i]
    txs <- txdf$TXID[idx]
    #all the xons for these transcripts
    ebt2 <- ebt[txs]
    df <- granges2df(ebt2)
    df$gene <- a@elementMetadata$gene_id[i]
    df$symbol <-  txdf[match(df$gene, txdf$GENEID), ]$SYMBOL 
    d[[i]] <- df
  }
  
  d <- do.call(rbind, d) 
  grt <- GeneRegionTrack(d)
  
  
  #plot it
  #ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = chr)
  axisTrack <- GenomeAxisTrack()
        plotTracks(c(dTrack, axisTrack, grt) , from = start, to = end, transcriptAnnotation="symbol", collapseTranscripts = 'longest', title.width=0.5, margin=40
                  )
        
}


plot_locus(list_of_files = list_of_files, trait_names =gwas_names[c(4,5,6)], start =10850255, end = 11800000, chr = 16, y_b= 35, locus_name = 205)







