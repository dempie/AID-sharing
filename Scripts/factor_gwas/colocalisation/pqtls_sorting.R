##### =========================== #####

### Script setup

##### =========================== #####

library(data.table)
library(R.utils)
library(dplyr)
library(stringr)
#library(liftOver)
#library(rtracklayer)
#path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
#ch = import.chain(path)

args = commandArgs(trailingOnly=TRUE)
i=as.numeric(args[1])

print(i)

proteomicsDir = "/processing_data/shared_datasets/plasma_proteome/decode/assocs_filtered/"
biomart_file = "data/decode_prot_locmap.tsv"
outDir = "/processing_data/shared_datasets/plasma_proteome/decode/assocs_filtered/cis_eqtls"
mappa=fread("data/UKBB_30k_map.tsv")


##### =========================== #####

### RAN ONCE TO EXTRACT CIS SIGNALS FOR END-POINTS

##### =========================== #####

 
library(biomaRt)



#
#decode_proteomics = system(paste("ls ",proteomicsDir,"/*.gz",sep=""),intern=T)
#
#filenames = basename(decode_proteomics)
#splitName = str_split_fixed(gsub(".txt.gz", "", filenames), "_", 4)


#
#hsapiens = useMart("ENSEMBL_MART_ENSEMBL",
#				"hsapiens_gene_ensembl")
#
#filters = listFilters(hsapiens)
#
#biomart_frame = getBM(attributes = c("entrezgene_id", "external_gene_name",
#									"chromosome_name", "start_position", "end_position"),
#					filters = "external_gene_name",
#					values = "LECT1",
#						mart = hsapiens)
#
#biomart_frame = biomart_frame[!(grepl("CHR_", biomart_frame$chromosome_name)),]
#
#fwrite(biomart_frame, "data/decode_prot_locmap.tsv",
# 		quote = FALSE, row.names = FALSE, sep = "\t")
biomart = fread(biomart_file, data.table = FALSE)
biomart=biomart[!(biomart$chromosome_name%in%c("Y","X")),]
allProteomics_files = Sys.glob(paste0(proteomicsDir, "*.gz"))
allProteomics_files_done = Sys.glob(paste0("/processing_data/shared_datasets/plasma_proteome/decode/assocs_filtered/cis_eqtls/", "*.tsv"))

allProteomics = str_split_fixed(basename(allProteomics_files), "_", 4)
allProteomics_done = str_split_fixed(basename(allProteomics_files_done), "_", 4)
missing=which(!(paste(allProteomics[,1],allProteomics[,2],allProteomics[,3],sep="_") %in% paste(allProteomics_done[,1],allProteomics_done[,2],allProteomics_done[,3],sep="_") ))

allProteomics_files=allProteomics_files[missing]
allProteomics = str_split_fixed(basename(allProteomics_files), "_", 4)



   gene.pos=biomart[biomart$external_gene_name==allProteomics[i,3],]
   if(nrow(gene.pos)>0){
   startPos=gene.pos$start_position-1000000
   end.pos=gene.pos$end_position+1000000
   chr=gene.pos$chromosome_name
   head=unlist(str_split(system(paste("zcat", allProteomics_files[i]," | head -n1"),intern=T),pattern="\t"))
   line=paste0("zcat ",allProteomics_files[i]," | awk -F \"\t\" '{ if(($1 == \"chr",chr,"\" && $2>",startPos," && $2<",end.pos,")) { print } }'")
   res=fread(line)
   names(res)=head
   out.names=paste0("/processing_data/shared_datasets/plasma_proteome/decode/assocs_filtered/cis_eqtls/",paste(allProteomics[i,1:3],collapse="_"),"_cis_pqtl.tsv")
   fwrite(res,file=out.names,row.names=FALSE,quote=F,sep="\t")
   }else{
     
     print(paste("Gene",allProteomics[i,3],"not found in biomart"))
     
     
   }




