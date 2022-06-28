
library(data.table)
library(GenomicRanges)
library(biomaRt)
library(stringr)
library(qqman)
library(RColorBrewer)
library(Gviz)
library(EnsDb.Hsapiens.v75)
library(ggpubr)




factor_loci <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/factor_loci_moloc_closest_gene_info.txt', data.table = F) 
colo <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/moloc_factors/factor_loci_moloc_info.txt') 

#create the Granges object for plotting with karyotipe R
f_ranges <- list()
f_list <- list()
for(i in 1:3){
        tt <- c('f1','f2', 'f3')[i]
        f_list[[tt]] <- fread(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/',tt,'_munged_build37.txt'), data.table = F)
        f_ranges[[tt]] <- GRanges(seqnames = paste0('chr', f_list[[tt]]$CHR), IRanges(start= f_list[[tt]]$BP, end=f_list[[tt]]$BP, names = f_list[[tt]]$SNP ), pval=f_list[[tt]]$P) 
}
  
#add the gene names to the factor loci table
mart <- useDataset("hsapiens_gene_ensembl", useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")) #select the database to convert and maake sure it is build 37 
genes <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", 'hgnc_symbol'),values= factor_loci$closest_gene ,mart= mart)
genes$hgnc_symbol[!str_detect(genes$hgnc_symbol, '')] <- genes$ensembl_gene_id[!str_detect(genes$hgnc_symbol, '')]

factor_loci$gene_name <- rep(NA, nrow(factor_loci)) 
for(i in 1:nrow(factor_loci)){
      factor_loci[i ,]$gene_name <- genes[ genes$ensembl_gene_id==factor_loci[i ,]$closest_gene ,]$hgnc_symbol
} 


#add the q values to the factor loci list
factor_loci[,'Q_index_pvalue'] <- rep(0, nrow(factor_loci))
gwas_list <- list()
for(i in 1:3){
      #load the gwas 
      tt <- c('f1', 'f2', 'f3')[i]
      gwas_list[[tt]] <- fread(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/factor',i,'_gwas_final_withQindex.txt'), data.table = F)
      #add the qpval   
      for(k in 1:nrow(factor_loci[factor_loci$trait==tt, ])){
              to_take <- factor_loci[factor_loci$trait==tt,]$SNP[k]
              factor_loci[factor_loci$trait==tt & factor_loci$SNP==to_take,]$Q_index_pvalue <- gwas_list[[tt]][gwas_list[[tt]]$SNP==to_take, ]$Q_chisq_pval
      }
  
}

#----- plot the locus plot and the genes ---------------------------------------

plot_locus_genes <- function(list_of_files, trait_names, start, end, chr, y_a=0, y_b=20, locus_name='-', colore){
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
    tr_SNP_active <- tr_SNP_active[c(tr_SNP_active$CHR==chr & data.table::between(tr_SNP_active$BP, lower = start, upper = end)),]
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
    elementMetadata(gr_tplot)[['p']] <- -log10(to_p[match(gr_tplot@ranges@NAMES, to_p$SNP), ]$P)
    list_of_granges[[k]] <- gr_tplot
  }
  
  dTrack <- list()
  for(n in 1:length(to_plot)){
    tt <- trait_names[[n]]
    dTrack[[n]] <- DataTrack(list_of_granges[[n]], name = tt, start = start, end = end, chromosome = chr, genome = 'hg19', ylim=c(y_a, y_b),
                           grid=F, frame=T, baseline= -log10(5e-8), col.baseline= brewer.pal(9, 'Greys')[9], lty.baseline= 1, lwd.baseline=0.2, col.axis= "black",
                              col.title="black" , background.title = c(brewer.pal(3,'Blues')[1],brewer.pal(3, 'Greens')[1], brewer.pal(3,'Reds')[1])[n],
                              pch=21, col='black', fill=colore[n]
    )
    
  }
  
  #retrieve the info on the genes
  ref_genes <- genes(EnsDb.Hsapiens.v75)
  ref_genes <- ref_genes[ref_genes@elementMetadata$gene_biotype=='protein_coding',]
  ind <- findOverlaps(GRanges(seqnames = chr ,IRanges(start = start, end = end)), ref_genes, type = 'any' )
  a <- ref_genes[ind@to,]
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
  grt <- GeneRegionTrack(d, showTitle=F, fill='black', background.title='white', size=0.5)
  
  
  #plot it
  #ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = chr)
  axisTrack <- GenomeAxisTrack()
  plotTracks(c(dTrack, axisTrack, grt) ,from = start, to = end, transcriptAnnotation="symbol", collapseTranscripts = 'longest', title.width=1,
             type=c('p'), lwd.grid=0.05
  )
  
}

#------ plots faster manathan plots of -----------------------------------------

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
    tr_SNP_active <- tr_SNP_active[c(tr_SNP_active$CHR==chr & data.table::between(tr_SNP_active$BP, lower = start, upper = end)),]
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
    manhattan(to_p[[tp]], chr="CHR", bp="BP", snp="SNP", p="P" ,xlim=c(start, end), main=paste0(tp,'  chr', chr,'_', start,'_', end, '  locus name ', locus_name ), suggestiveline =F , 
              col = c("brown2"), highlight = c(to_p[[tp]][which.min(to_p[[tp]]$P), ]$SNP) )  
    grid(nx = NULL, ny = NULL,
         lty = 2,      # Grid line type
         col = "gray", # Grid line color
         lwd = 1)
    
  }
  
}

#-------------------------------------------------------------------------------
#plot examples


vai <- function(x){
  plot_locus(list_of_files = f_list, trait_names = list('f1', 'f2', 'f3'), start = min(factor_loci[factor_loci$gene_name==x,]$start )   ,end = max(factor_loci[factor_loci$gene_name==x,]$end ), chr = factor_loci[factor_loci$gene_name==x,]$chr,  locus_name = '') 
}

plot_locus(list_of_files = f_list, trait_names = list('f1', 'f2', 'f3'), start =204450734   ,end = 204804955 , chr = 2,  locus_name = '') #CTLA4
plot_locus(list_of_files = f_list, trait_names = list('f1', 'f2', 'f3'), start =191163706     ,end = 192504445 , chr = 2,  locus_name = '') #STAT1 ok 
plot_locus(list_of_files = f_list, trait_names = list('f1', 'f2', 'f3'), start = min(factor_loci[factor_loci$gene_name=='BACH2',]$start )   ,end = max(factor_loci[factor_loci$gene_name=='BACH2',]$end ), chr = factor_loci[factor_loci$gene_name=='BACH2',]$chr,  locus_name = '') 
vai('IFNG')
vai('CCR7')
vai('CCR4')
vai('TNFSF4') #ok, significant for f3 and f1
vai('CXCR5')
vai('IL2RA')
vai('JAK2') #only f1 nice
vai('TRAF6') #only f3 significant  
vai('BATF3')
vai('PRKCQ')
vai('ETS1')
vai('IL2RB')
vai('STAT1') #only F2  
vai('CTLA4') #only F2 but close in F1
vai('TNFSF15') 
vai('RORC')

plot_locus_genes( list_of_files = f_list, trait_names = list('f1', 'f2', 'f3'), start =191163706     ,end = 192504445 , chr = 2 ,colore = brewer.pal(5, 'Paired')[c(1,3,5)] ) #stat1 and stat4 gene significant for f2
plot_locus_genes( list_of_files = f_list, trait_names = list('f1', 'f2', 'f3'), start = 172000000    ,end =  174000000 , chr = 1 ,colore = brewer.pal(5, 'Paired')[c(1,3,5)] ) ##ok, significant for f3 and f1 TNFSF4
plot_locus_genes( list_of_files = f_list, trait_names = list('f1', 'f2', 'f3'), start = 4900000    ,end =  5400000 , chr = 9 ,colore = brewer.pal(5, 'Paired')[c(1,3,5)] ) ##ok, only f1 is JAK2
plot_locus_genes( list_of_files = f_list, trait_names = list('f1', 'f2', 'f3'), start = 36000000    ,end =  37000000 , chr = 11 ,colore = brewer.pal(5, 'Paired')[c(1,3,5)] ) ##ok, only f3 is TRAF6
plot_locus_genes( list_of_files = f_list, trait_names = list('f1', 'f2', 'f3'), start = 37000000    ,end =  38000000 , chr = 22 ,colore = brewer.pal(5, 'Paired')[c(1,3,5)] ) #IL2RB significant for f2 anf f3  
plot_locus_genes( list_of_files = f_list, trait_names = list('f1', 'f2', 'f3'), start = 204300000    ,end =  205300000 , chr = 2 ,colore = brewer.pal(5, 'Paired')[c(1,3,5)],  y_a = 0, y_b = 25  ) #only F2 but close in F1 CTLA$
plot_locus_genes( list_of_files = f_list, trait_names = list('f1', 'f2', 'f3'), start = 117000000    ,end =  118000000 , chr = 9 ,colore = brewer.pal(5, 'Paired')[c(1,3,5)] ) #only F1 TNFSF15
plot_locus_genes( list_of_files = f_list, trait_names = list('f1', 'f2', 'f3'), start = 150000000    ,end =  152000000 , chr = 1 ,colore = brewer.pal(5, 'Paired')[c(1,3,5)] ) #rorc



#--------- the selected ones ---------------------------------------------------

#stat1 and stat4 gene significant for f2
pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure1/locus_zoom/locus_chr2_191163706_192004445_stat1_f2.pdf', width = 14, height = 10)
plot_locus_genes( list_of_files = f_list, trait_names = list('f1', 'f2', 'f3'), start =191200000     ,end = 192200000 , chr = 2 ,colore = brewer.pal(5, 'Paired')[c(1,3,5)], y_a = 0, y_b = 20 ) #stat1 and stat4 gene significant for f2
dev.off()

pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure1/locus_zoom/locus_chr9_4963838_5270603_jak2_f1.pdf', width = 14, height = 10)
plot_locus_genes( list_of_files = f_list, trait_names = list('f1', 'f2', 'f3'), start = 4850000    ,end =  5400000 , chr = 9 ,colore = brewer.pal(5, 'Paired')[c(1,3,5)] , y_a = 0, y_b = 20) ##ok, only f1 is JAK2
dev.off()

pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure1/locus_zoom/locus_chr11_36336263_36530644_traf6_f3.pdf', width = 14, height = 10)
plot_locus_genes( list_of_files = f_list, trait_names = list('f1', 'f2', 'f3'), start = 36200000, end = 36650000 , chr = 11 ,colore = brewer.pal(5, 'Paired')[c(1,3,5)] ,  y_a = 0, y_b = 15) ##ok, only f3 is TRAF6
dev.off()

pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure1/locus_zoom/locus_chr2_204450734_204804955_ctla4_f2.pdf', width = 14, height = 10)
plot_locus_genes( list_of_files = f_list, trait_names = list('f1', 'f2', 'f3'), start = 204450000    ,end =  204950000 , chr = 2 ,colore = brewer.pal(5, 'Paired')[c(1,3,5)],  y_a = 0, y_b = 25  ) #only F2 but close in F1 CTLA$
dev.off()


pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure1/locus_zoom/locus_chr1_67392496_68018455_c1orf141_il23r_f1.pdf', width = 14, height = 10)
plot_locus_genes( list_of_files = f_list, trait_names = list('f1', 'f2', 'f3'), start = 67200000   ,end = 68100000 , chr = 1 ,colore = brewer.pal(5, 'Paired')[c(1,3,5)],  y_a = 0, y_b = 35  ) #only F1, c1orf141_il23r
dev.off()


pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure1/locus_zoom/locus_chr17_47222663_47677191_ZNF652_f3.pdf', width = 14, height = 10)
plot_locus_genes( list_of_files = f_list, trait_names = list('f1', 'f2', 'f3'), start = 47000000 ,end = 47900000 , chr = 17 ,colore = brewer.pal(5, 'Paired')[c(1,3,5)],  y_a = 0, y_b = 25 ) #only F3 ZNF652
dev.off()

pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure1/locus_zoom/locus_chr3_48541016_50290016_apeh_f1.pdf', width = 14, height = 10)
plot_locus_genes( list_of_files = f_list, trait_names = list('f1', 'f2', 'f3'), start = 48400000,end =  50300000 , chr = 3 ,colore = brewer.pal(5, 'Paired')[c(1,3,5)],  y_a = 0, y_b = 30 ) #only F1, apeh, huge locus
dev.off()


pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure1/locus_zoom/locus_chr8_81092889_81325152_tpd52_f3.pdf', width = 14, height = 10)
plot_locus_genes( list_of_files = f_list, trait_names = list('f1', 'f2', 'f3'), start = 80500000,end =  81400000 , chr = 8 ,colore = brewer.pal(5, 'Paired')[c(1,3,5)],  y_a = 0, y_b = 25 ) #only F3
dev.off()


#q index significant

pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure1/locus_zoom/locus_chr1_1_113759855_114651076_bcl2l15_ptpn22_f2_Q_index_significant.pdf', width = 14, height = 10)
plot_locus_genes( list_of_files = f_list, trait_names = list('f1', 'f2', 'f3'), start = 113600000    ,end =114900000 , chr = 1 ,colore = brewer.pal(5, 'Paired')[c(1,3,5)],  y_a = 0, y_b = 35  ) #only F2 ptpn22 and bcl2l15, Q index significant
dev.off()


