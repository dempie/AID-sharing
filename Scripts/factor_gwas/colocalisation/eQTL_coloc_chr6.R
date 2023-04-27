source("scripts/Multi_coloc_funtions.R")



args = commandArgs(trailingOnly=TRUE)

an.n=as.numeric(args[1])

ct=c("BIN"
,"BMem"
,"CD4ET"
,"CD4NC"
,"CD4SOX4"
,"CD8ET"
,"CD8NC"
,"CD8S100B"
,"DC"
,"MonoC"
,"MonoNC"
,"NK"
,"NKR"
,"Plasma")

an.list=c()
for(i in 6){
  
  an.list=rbind(an.list,cbind(i,ct))
  
}



chr=as.numeric(an.list[an.n,1])
cell.type=an.list[an.n,2]

#eqtl.map=fread("data/onek1k_eqtl_dataset.tsv",skip=1,select = c(3,4,7,8,9,10,11))
#eqtl.map=unique(eqtl.map)
#
#names(eqtl.map)=c("SNP","snpid","CHR","BP","A1","A2","FreqA2")
#
#write.table(eqtl.map,file="data/One1k_snp.info",row.names=F,quote=F,sep="\t")

loci.table=fread("output/loci_definitions/final_locus_table.tsv")
loci.table=as.data.frame(loci.table)
loci.table=loci.table[loci.table$Chr==chr,]

loci.table=loci.table[which(!(loci.table$pan.locus%in%c(101, 104, 105))),]


mappa=fread("data/UKBB_30k_map.tsv")

mappa=mappa[mappa$CHR==chr,]

loci.table=loci.table[order(loci.table$start),]
loci.table$final_locus=paste(loci.table$pan.locus,loci.table$sub_locus,sep="_")

loci.table=loci.table[order(loci.table$final_locus),]

for(i in loci.table$pan_locus){
  
  idx=which(loci.table$pan_locus==i)
  loci.table$start[idx]=min(loci.table$start[idx])
  loci.table$end[idx]=max(loci.table$end[idx])
  
}





colocalization.table.all=c()
colocalization.table.H4=c()

sc.res=fread(paste0("data/onek_scRNA/round1/",cell.type,"_chr",chr,"_correlation_results.tsv"))

eqtl.map=fread("data/One1k_snp.info")
sc.res$SNP=eqtl.map$SNP[match(sc.res$snpid,eqtl.map$snpid)]
sc.res$A1=eqtl.map$A1[match(sc.res$snpid,eqtl.map$snpid)]
sc.res$A2=eqtl.map$A2[match(sc.res$snpid,eqtl.map$snpid)]
sc.res$FREQ=eqtl.map$FreqA2[match(sc.res$snpid,eqtl.map$snpid)]
sc.res$BP=eqtl.map$BP[match(sc.res$snpid,eqtl.map$snpid)]



bfile="/processing_data/shared_datasets/ukbiobank/genotypes/LD_reference/ld_reference_bfiles/ukbb_all_30000_random_unrelated_white_british"
file.list=system("ls /project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/*.txt",intern=T)
labels=gsub("/project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/","",
            gsub("_munged_build37.txt","",file.list))
reference.map=mappa
reference.map=as.data.frame(reference.map)
file.list=file.list
labels=labels
t.type=c("cc","cc","cc","quant","quant","quant","cc","cc","cc","cc","cc","cc")
prev.var=c(0.2696713,0.434383,0.02902916,1,1,1,0.3593954,0.2388718,0.3269585,0.5736819,0.0377603,0.3679372)
file.table=data.frame(file=file.list,trait=labels,type=t.type,prev.var=prev.var)







for(locus in 1:nrow(loci.table)){
  
  
  loci.table.tmp=loci.table[locus,]
  locus.info=loci.table.tmp
  
  
  start=min(loci.table.tmp$start)-100000
  end=max(loci.table.tmp$end)+100000
  n.table=c()
  
  eqtl.results=sc.res[sc.res$BP>=start & sc.res$BP<=end,]
  
  if(min(eqtl.results$p.value)<=1e-6){
  
    datasets=list()
    mappa.loc=reference.map[which(reference.map$CHR==chr & reference.map$BP>=start & reference.map$BP<=end),]
    for(i in unique(loci.table.tmp$trait)){
      
      print(dim(mappa.loc))
      dataset=fread(file.table$file[file.table$trait==i],data.table = FALSE)
      dataset=dataset[which(dataset$SNP%in% mappa.loc$SNP), ]
      
      if(!("MAF" %in% colnames(dataset)) | all(is.na((dataset$MAF)))){
        if("FRQ" %in% colnames(dataset)){
          
          dataset$MAF=dataset$FRQ
          dataset$MAF[dataset$MAF>0.5]=1-dataset$MAF[dataset$MAF>0.5]
          
        }else{
          
          dataset$MAF=reference.map$MAF[match(dataset$SNP,reference.map$SNP)]
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
    
    gene.list=unique(eqtl.results$geneid)
    
    
    for(i in gene.list){   
      
      
      dataset=eqtl.results[eqtl.results$geneid==i,]
      dataset=as.data.frame(dataset)
      if(any(dataset$p.value<1e-6)){
        
        dataset$MAF=dataset$FREQ
        dataset$MAF[dataset$MAF>0.5]=1-dataset$MAF[dataset$MAF>0.5]
        chis=qchisq(dataset$p.value,df=1,lower=F)
        dataset$SE=dataset$estimate/sqrt(chis)
        
        
        N_hat<-median(1/((2*dataset$MAF*(1-dataset$MAF))*dataset$SE^2),na.rm = T)
        dataset$N=ceiling(N_hat)
        mappa.loc=reference.map[which(reference.map$CHR==chr & reference.map$BP>=start & reference.map$BP<=end),]
        
        dataset=dataset[match(mappa.loc$SNP,dataset$SNP),]
        dataset$CHR=chr
        flip=dataset[,c("SNP","CHR","BP","A2","A1","estimate")]
        names(flip)=c("rsid","chr","pos","a0","a1","beta")
        names(mappa.loc)=c("rsid","chr","pos","maf","a1","a0")
        flip.t=snp_match(sumstats = flip,info_snp = mappa.loc,join_by_pos=FALSE,strand_flip = FALSE,match.min.prop=0)
        flip.t=flip.t[match(dataset$SNP,flip.t$rsid),]
        nadataset=dataset[match(flip.t$rsid,dataset$SNP),]
        
        dataset$A1=flip.t$a1
        dataset$A2=flip.t$a0
        dataset$BETA=flip.t$beta
        dataset$varbeta=dataset$SE^2
        dataset=dataset[,c("SNP","CHR","BP","A1","A2","BETA","varbeta","p.value","MAF","N")]
        dataset$type='quant'
        print("Yes")
        dataset$sdY=1
        names(dataset)=c("snp","chr","pos","a1","a0","beta","varbeta","pvalues","MAF","N","type","sdY")
        dataset=na.omit(dataset)
        datasets[[i]]=dataset
      }
      
    }
    
    
    
    if(length(datasets)>1){
      conditional.datasets=list()
      max.loci=1
      for(i in 1:length(datasets)){
        
        tmp=cojo.ht(D=datasets[[i]],p.tresh = 1e-4,bfile = bfile)
        conditional.datasets[[i]]=tmp
        names(conditional.datasets)[i]=names(datasets)[i]
        max.loci=max(max.loci,nrow(tmp$ind.snps))
        
      }
      
      pairwise.list=cbind(names(conditional.datasets)[1],names(conditional.datasets)[-1])
      final.colocs=c()
      all.colocs=list()
      k=1
      for(i in 1:nrow(pairwise.list)){
        
        coloc.res=colo.cojo.ht(conditional.datasets[[pairwise.list[i,1]]],conditional.datasets[[pairwise.list[i,2]]])
        coloc.res$t1=pairwise.list[i,1]
        coloc.res$t2=pairwise.list[i,2]
        final.colocs=rbind(final.colocs,coloc.res)
      }
      final.colocs.H4=final.colocs[round(final.colocs$PP.H4.abf)>0.95,]
      final.colocs.H4=as.data.frame(final.colocs.H4)
      if(nrow(final.colocs.H4)>0){
        final.colocs.H4$locus=locus
        colocalization.table.H4=c(colocalization.table.H4,final.colocs.H4)
      }
      final.colocs$locus=loci.table.tmp$final_locus
      colocalization.table.all=rbind(colocalization.table.all,final.colocs)
    }
  }  
  
  
}
colocalization.table.all$cell_type=cell.type

write.table(colocalization.table.all,file=paste0("output/sc_eqtl/Chr",chr,"_Cell_type_",cell.type,".tsv"),row.names=F,quote=F,sep="\t")









