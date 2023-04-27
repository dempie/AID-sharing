library(data.table)
library(liftOver)
library(rtracklayer)
library(R.utils)
library(dplyr)
library(stringr)

path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)


biomart_file = "data/decode_prot_locmap.tsv"

prot.map=fread(biomart_file)
mappa=fread("data/UKBB_30k_map.tsv")
allProteomics_files_done = Sys.glob(paste0("/processing_data/shared_datasets/plasma_proteome/decode/assocs_filtered/cis_eqtls/", "*.tsv"))
allProteomics_done = str_split_fixed(basename(allProteomics_files_done), "_", 4)



for(i in 1:22){
  chr.data=c()
  prot.map.tmp=prot.map[prot.map$chromosome_name==i,]
  file.idx=which(allProteomics_done[,3]%in%prot.map.tmp$external_gene_name)
  allProteomics_files_done.tmp=allProteomics_files_done[file.idx]
  allProteomics_done.tmp=allProteomics_done[file.idx,]
  for(j in 1:length(allProteomics_files_done.tmp)){
    
    
    pt.data=fread(allProteomics_files_done.tmp[j])
    pt.data$gene=allProteomics_done.tmp[j,3]
    names(pt.data)[1:2]=c("chrom","start")
    pt.data$end=pt.data$start
    grang=as(pt.data, "GRanges")
    cur19 = liftOver(grang, ch)
    cur19=as.data.frame(cur19)
    names(cur19)[3:4]=c("Chrom","Pos")
    chr.data=rbind(chr.data,cur19)
     
  }
  fwrite(chr.data,file=paste0("/processing_data/shared_datasets/plasma_proteome/decode/assocs_filtered/cis_eqtls/Chr",i,"_cis_pqtls.tsv"),row.names=F,quote=F,sep="\t")
  
}

