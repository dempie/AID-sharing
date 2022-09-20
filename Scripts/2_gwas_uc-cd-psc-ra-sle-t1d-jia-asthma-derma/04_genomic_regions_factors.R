#in this scirpt I will perform downstrean analysis of the GWAS results.
library(data.table)
library(GenomicRanges)
library(dplyr)
library(qqman)
library(ChIPpeakAnno)
library(RColorBrewer)
library(ComplexHeatmap)

#----locus.breaker function by Nicola, it identifies loci-----------------------
locus.breaker=function(res,p.sig=5e-8, p.limit=1e-5,hole.size=250000
                       ,p.label="p",chr.label="chr",pos.label="pos"){
  
  res = res[order(as.numeric(res[, chr.label]), as.numeric(res[,pos.label])),]
  
  res=res[which(res[,p.label]<p.limit),]
  trait.res=c()
  for(j in 1:22){
    
    res.chr=res[which(res[,chr.label]==j),]
    if(nrow(res.chr)>1){
      holes=res.chr[,pos.label][-1]-res.chr[,pos.label][-length(res.chr[,pos.label])] 
      gaps=which(holes>hole.size)
      if(length(gaps)>0){
        for(k in 1:(length(gaps)+1)){
          
          if(k==1){
            res.loc=res.chr[1:(gaps[k]),]  
          }else if(k==(length(gaps)+1)){
            res.loc=res.chr[(gaps[k-1]+1):nrow(res.chr),]  
          }else{
            res.loc=res.chr[(gaps[k-1]+1):(gaps[k]),]
          }
          if(min(res.loc[,p.label])<p.sig){
            
            start.pos=min(res.loc[,pos.label],na.rm=T)
            end.pos=max(res.loc[,pos.label],na.rm=T)
            chr=j
            best.snp=res.loc[which.min(res.loc[,p.label]),]
            line.res=c(chr,start.pos,end.pos,unlist(best.snp))
            trait.res=rbind(trait.res,line.res)
          }
          
          
        }
      }else{
        res.loc=res.chr
        if(min(res.loc[,p.label])<p.sig)  {
          
          start.pos=min(res.loc[,pos.label],na.rm=T)
          end.pos=max(res.loc[,pos.label],na.rm=T)
          chr=j
          best.snp=res.loc[which.min(res.loc[,p.label]),]
          line.res=c(chr,start.pos,end.pos,unlist(best.snp))
          trait.res=rbind(trait.res,line.res)
        }
        
      }
      
    }else if(nrow(res.chr)==1){
      
      res.loc=res.chr
      if(min(res.loc[,p.label])<p.sig){
        start.pos=min(res.loc[,pos.label],na.rm=T)
        end.pos=max(res.loc[,pos.label],na.rm=T)
        chr=j
        best.snp=res.loc[which.min(res.loc[,p.label]),]
        line.res=c(chr,start.pos,end.pos,unlist(best.snp))
        trait.res=rbind(trait.res,line.res)
      }
      
      
    }
  }
  
  print(trait.res)
  trait.res=as.data.frame(trait.res,stringsAsFactors=FALSE)
  trait.res=trait.res[,-(which(names(trait.res)==chr.label))]
  names(trait.res)[1:3]=c("chr","start","end")
  trait.res
}


#------ run locus_breaker for defining the genomic regions----------------------
sstats <- list()
regions <- list()
for(i in c('f1', 'f2', 'f3')){
  
  sstats[[i]]<- fread(paste0('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/03_estimation/summarystats_',i,'_.txt'), data.table = F)
  #run locus breaker on each of the factors summary stats
  regions[[i]] <- locus.breaker(res=sstats[[i]], p.label = 'P', chr.label = 'CHR', pos.label = 'BP')
  regions[[i]][, 'trait'] <- i 
  
}

regions_factors <- do.call(rbind, regions)
rownames(regions_factors) <-NULL 

regions_factors
fwrite(regions_factors, 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/04_genomic_regions_factors/regions_factors.txt', row.names = F, col.names = T, sep = '\t')


#------- search for ovelapping genomic regions between factors -----------------
f_ranges <- list()
for(i in c('f1', 'f2', 'f3')){

f_ranges[[i]] <- GRanges(seqnames = regions_factors[regions_factors$trait==i,]$chr, IRanges(as.numeric(regions_factors[regions_factors$trait==i,]$start),as.numeric(regions_factors[regions_factors$trait==i,]$end)))

}

#plot venndiagram
pdf(width = 14, height = 14, file = 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/04_genomic_regions_factors/venn_f1_f2_f3_loci.pdf')
venn_f1_f2_f3_loci<- makeVennDiagram(Peaks=f_ranges,connectedPeaks = 'keepAll',
                                     NameOfPeaks=c("F1", "F2", "F3"),
                                     imagetype="png" ,
                                     height = 480 , 
                                     width = 480 , 
                                     resolution = 300,
                                     compression = "lzw",
                                     lwd = 1,
                                     cex = 1.5,
                                     cat.cex = 2,
                                     cat.default.pos = "outer",
                                     cat.pos = c(-27, 27, 135),
                                     cat.dist = c(0.055, 0.055, 0.085),
                                     cat.fontfamily = "sans",
                                     rotation = 1
)

dev.off()

#----- plot UpSet plot ---------------------------------------------------------

#function for obtaining a combination matrix for plotting and UpSet plot from Granges names
lists_forupset <- function(peaks, traits){ 
  #allocate lists
  require(ComplexHeatmap)
  require(ChIPpeakAnno)
  require(stringr)
  
  ovl_f <- findOverlapsOfPeaks(peaks)
  
  un <- list()
  final <- list(list())
  sh <- list()
  #create the lists and put it the elements that are unique for each trait
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
        un[[tt]][[length(un[[tt]])+1]] <- paste0(ovl_f$mergedPeaks$peakNames[[k]], collapse = '-')
      }
      
    }
  }
  
  #generate the output
 
  cm <- make_comb_mat(un, mode = 'distinct')
  return(cm)
  
}


#plot

cm <- lists_forupset(f_ranges, traits = c('f1', 'f2', 'f3')) 

pdf(width = 14, height = 14, file = 'outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/04_genomic_regions_factors/upset_f1_f2_f3_loci.pdf')
UpSet(cm, set_order = c("f1", "f2", "f3"), comb_order = order(comb_size(cm), decreasing = F),
      comb_col = c(brewer.pal(8, 'Set2')[c(7,6,6,6)],brewer.pal(12, 'Paired')[c(2,4,6)]),
      top_annotation = upset_top_annotation(cm, add_numbers = TRUE, height = unit(6, "cm")),
      left_annotation = upset_left_annotation(cm, add_numbers = TRUE, width = unit(3,'cm'), gp = gpar(fill = brewer.pal(5, 'Greys')[1])  ),
      row_title = "Factor", 
      column_title = " Physical intersection of factor loci"
)
dev.off()


#---------------- unique new loci ----------------------------------------------







