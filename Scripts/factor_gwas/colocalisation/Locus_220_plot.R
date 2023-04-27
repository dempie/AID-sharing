source("scripts/Multi_coloc_funtions.R")

loci.table=fread("/project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/hyprcoloc_all_against_all/loci_all_traits.txt")
loci.table=loci.table[order(loci.table$pan_locus),]
for(i in loci.table$pan_locus){
  
  idx=which(loci.table$pan_locus==i)
  loci.table$start[idx]=min(loci.table$start[idx])
  loci.table$end[idx]=max(loci.table$end[idx])
  
}



#mappa=fread("/processing_data/shared_datasets/ukbiobank/genotypes#/LD_reference/ld_reference_bfiles#/ukbb_all_30000_random_unrelated_white_british.bim")
#
#freq=fread("/processing_data/shared_datasets/ukbiobank/genotypes#/LD_reference/ld_reference_bfiles#/ukbb_all_30000_random_unrelated_white_british.frqx")
#mappa$freq=(freq$`C(HOM A1)`*2+freq$`C(HET)`)/(rowSums(freq[,5:7])*2)
#names(mappa)=c("CHR","SNP","","BP","A1","A2","FREQ")
#mappa=mappa[,-3]
#mappa$MAF[mappa$MAF>0.5]=1-mappa$MAF[mappa$MAF>0.5]
#mappa=mappa[,c("SNP","CHR","BP","MAF","A1","A2")]
#write.table(mappa,file="data/UKBB_30k_map.tsv",row.names=F,quote=F,sep="\t")


mappa=fread("data/UKBB_30k_map.tsv")



colocalization.table.all=c()
colocalization.table.H4=c()
bfile="/processing_data/shared_datasets/ukbiobank/genotypes/LD_reference/ld_reference_bfiles/ukbb_all_30000_random_unrelated_white_british"

## GWAS info
file.list=system("ls /project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/*.txt",intern=T)
labels=gsub("/project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/","",gsub("_munged_build37.txt","",file.list))

proportions=readRDS("/project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/04_sumstat_outputs/case_controls_fraction_nic.RDS")

reference.map=mappa
file.list=file.list
labels=labels
t.type=c("cc","cc","cc","quant","quant","quant","cc","cc","cc","cc","cc","cc")
prev.var=c(0.2696713,0.434383,0.02902916,1,1,1,0.3593954,0.2388718,0.3269585,0.5736819,0.0377603,0.3679372)

file.table=data.frame(file=file.list,trait=labels,type=t.type,prev.var=prev.var)
col.order=c("trait","Chr",  "start", "end", "SNP" ,  "bp","refA","othA","freq","b","se","p","bJ","bJ_se","pJ","LD_r" ,"n" , "pan.locus" , "sub_locus" )
final.locus.table=c()

locus=220
#for( locus in unique(loci.table$pan_locus)){
  
  loci.table.tmp=loci.table[loci.table$pan_locus==locus,]
  locus.info=loci.table.tmp

  #Define genomic region
  start=min(loci.table.tmp$start)-100000
  end=max(loci.table.tmp$end)+100000
  chr=loci.table.tmp$chr[1]
  n.table=c()
  ##### load and harmonise the SNP data
  
  datasets=list()
  mappa.loc=mappa[which(mappa$CHR==chr & mappa$BP>=start & mappa$BP<=end),]
  for(i in unique(loci.table.tmp$trait)){
    
    print(dim(mappa.loc))
    dataset=fread(file.table$file[file.table$trait==i])
    if(!("MAF" %in% colnames(dataset)) | all(is.na((dataset$MAF)))){
      if("FRQ" %in% colnames(dataset)){
        
        dataset$MAF=dataset$FRQ
        dataset$MAF[dataset$MAF>0.5]=1-dataset$MAF[dataset$MAF>0.5]
        
      }else{
        
        dataset$MAF=mappa$MAF[match(dataset$SNP,mappa$SNP)]
      }
    }
    
    
    idx=which(dataset$MAF>0.1 & dataset$MAF<0.4)
    if(!(any(names(dataset)%in%"N"))){
      N_hat<-median(1/((2*dataset$MAF[idx]*(1-dataset$MAF[idx]))*dataset$SE[idx]^2),na.rm = T)
      dataset$N=ceiling(N_hat)
    }
    mappa.loc=reference.map[which(reference.map$CHR==chr & reference.map$BP>=start & reference.map$BP<=end),]
    dataset=dataset[match(mappa.loc$SNP,dataset$SNP),]
    flip=dataset[,c("SNP","CHR","BP","A2","A1","BETA")]
    names(flip)=c("rsid","chr","pos","a0","a1","beta")
    names(mappa.loc)=c("rsid","chr","pos","maf","a1","a0")
    flip.t=snp_match(sumstats = flip,info_snp = mappa.loc,join_by_pos=FALSE,strand_flip = FALSE,match.min.prop=0)
    dataset=dataset[match(flip.t$rsid,dataset$SNP),]
    dataset$A1=flip.t$a1
    dataset$A2=flip.t$a0
    dataset$BETA=flip.t$beta
    dataset$varbeta=dataset$SE^2
    dataset=dataset[,c("SNP","CHR","BP","A1","A2","BETA","varbeta","P","MAF","N")]
    dataset$type=file.table$type[file.table$trait==i]
    if(file.table$type[file.table$trait==i]=="cc"){
      
      dataset$s=file.table$prev.var[file.table$trait==i]
      names(dataset)=c("snp","chr","pos","a1","a0","beta","varbeta","pvalues","MAF","N","type","s")
      
    }else if(file.table$type[file.table$trait==i]=="quant"){
      
      dataset$sdY=file.table$prev.var[file.table$trait==i]
      names(dataset)=c("snp","chr","pos","a1","a0","beta","varbeta","pvalues","MAF","N","type","sdY")
      
      
    }else{
      
      stop("Type has to be either 'cc' or 'quant'")
      
    }
    
    datasets[[i]]=dataset
    
    
  }
  
  datasets$f1$trait="F1"
  datasets$f1$log.pvalues=-log10(datasets$f1$pvalues)

  datasets$f2$trait="F2"
  datasets$f2$log.pvalues=-log10(datasets$f2$pvalues)

  datasets$f3$trait="F3"
  datasets$f3$log.pvalues=-log10(datasets$f3$pvalues)

  factor.plot=rbind(datasets$f1,datasets$f2,datasets$f3)
  p1=ggplot(factor.plot,aes(x=pos,y=log.pvalues,fill=trait,color=trait))+
    geom_point(shape=21,alpha=1,size=7)+scale_color_manual(values = c("#114466" ,"#1a5217" ,"#8a1112"))+
    theme_minimal()+facet_wrap(trait~.,nrow = 3, scales = "free")+scale_fill_manual(values = c("#1F78B4" ,"#33A02C" ,"#E31A1C"))
  
  
   
  
  
  
  
  ### Conditional datasets
  conditional.datasets=list()
  max.loci=1
  for(i in 1:length(datasets)){
    
    tmp=cojo.ht(D=datasets[[i]],p.tresh = 1e-4,bfile = bfile)
    conditional.datasets[[i]]=(tmp)
    names(conditional.datasets)[i]=names(datasets)[i]
    max.loci=max(max.loci,nrow(tmp$ind.snps))
    
  }
  
  conditional.datasets$f1$ind.snps[which(conditional.datasets$f1$ind.snps$pJ<1e-6 | conditional.datasets$f1$ind.snps$p<5e-8),]
  
  
  conditional.datasets$f1$results$rs12922863$trait="F1"
  conditional.datasets$f1$results$rs12922863$sub_trait="F1-rs12922863"
  
  conditional.datasets$f1$results$rs416603$trait="F1"
  conditional.datasets$f1$results$rs416603$sub_trait="F1-rs416603"
 
  conditional.datasets$f1$results$rs13335254$trait="F1"
  conditional.datasets$f1$results$rs13335254$sub_trait="F1-rs13335254"
  
  conditional.datasets$f2$ind.snps[which(conditional.datasets$f2$ind.snps$pJ<1e-6 | conditional.datasets$f2$ind.snps$p<5e-8),]
  
  
  
  conditional.datasets$f2$results$rs34306440$trait="F2"
  conditional.datasets$f2$results$rs34306440$sub_trait="F2-rs34306440"
  
  conditional.datasets$f2$results$rs12928968$trait="F2"
  conditional.datasets$f2$results$rs12928968$sub_trait="F2-rs12928968"
  
  conditional.datasets$f3$ind.snps[which(conditional.datasets$f3$ind.snps$pJ<1e-6 | conditional.datasets$f3$ind.snps$p<5e-8),]

  conditional.datasets$f3$results$rs2110605$trait="F3"
  conditional.datasets$f3$results$rs2110605$sub_trait="F3-rs2110605"
  
  conditional.datasets$f3$results$rs12925474$trait="F3"
  conditional.datasets$f3$results$rs12925474$sub_trait="F3-rs12925474"
  
  conditional.datasets$f3$results$rs7206912$trait="F3"
  conditional.datasets$f3$results$rs7206912$sub_trait="F3-rs7206912"
  
   
  
  
  splitup.results=rbind(conditional.datasets$f1$results$rs12922863
                        ,conditional.datasets$f1$results$rs416603
                        ,conditional.datasets$f1$results$rs13335254
                        ,conditional.datasets$f2$results$rs34306440
                        ,conditional.datasets$f2$results$rs12928968
                        ,conditional.datasets$f3$results$rs2110605
                        ,conditional.datasets$f3$results$rs12925474
                        ,conditional.datasets$f3$results$rs7206912)

  p2=ggplot(splitup.results,aes(x=bp,y=-log10(as.numeric(pC)),fill=trait,color=trait))+
    geom_point(shape=21,alpha=1,size=5)+geom_hline(yintercept = 0)+scale_color_manual(values = c("#114466" ,"#1a5217" ,"#8a1112"))+
    theme_minimal()+facet_wrap(sub_trait~.,ncol = 1, scales = "free")+scale_fill_manual(values = c("#1F78B4" ,"#33A02C" ,"#E31A1C"))+coord_cartesian(xlim=c(start, end))

  library(patchwork)
  p1+p2
    
  
  
  
  
  #plot_locus_genes <- function(start, end, chr){
    require(Gviz)
    require(EnsDb.Hsapiens.v75)
    require(GenomicRanges)
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
    

  
  gene.starts=by(d$start,d$symbo,FUN = min)
  gene.end=by(d$end,d$symbo,FUN = max)
  g.table=as.data.frame(cbind(gene.starts,gene.end ))
  g.table=g.table[order(g.table$gene.starts),]
  g.table$y=1
  for(i in 2:nrow(g.table)){
    
    if(g.table$gene.starts[i]<g.table$gene.end[i-1] & g.table$y[i-1]==1){
      g.table$y[i]=-1
    }
    
  }
  g.table$col=c(1,2,1,2,1,2,1,2,1,2,1,2,1,2,1)
  g.table["TNP2",]$y=-1
  g.table["PRM3",]$y=-1
  g.table["PRM2",]$y=-1
  g.table["PRM1",]$y=-1
  lab.table=rowMeans(g.table[,1:2])
  d$y=g.table$y[match(d$symbol,row.names(g.table))]
  library(ggrepel)
  p3=ggplot()+geom_segment(data=g.table,aes(x =gene.starts , xend = gene.end, y = y, yend = y),size=2, 
                                         color = "blue", alpha = 1)+coord_cartesian(xlim=c(start, end),ylim=c(-1.8,1.8))+
    geom_rect(data=d,aes(xmin =start , xmax= end, ymin = y+0.5, ymax = y-0.5), 
                 fill = "black", alpha = 0.4, color = "black")+theme_void()+
    geom_text_repel(aes(y=g.table$y+0.6,x=lab.table,label=names(lab.table))) 
  
    
  
  library(patchwork)
  p4=p1/p3+plot_layout(heights  = c(4, 1))
  p5=p2+plot_layout(heights  = c(6, 1))
  
  pdf("Locus_220_decomposition.pdf",height=21,width=30)
  ((p1/p3+plot_layout(heights  = c(9, 1)))|p2)+plot_layout(widths = c(3,1),guides = "collect")
  dev.off()
  
  
  
  g.table=as.data.frame(cbind(gene.starts,gene.end ))
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  