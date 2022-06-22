library(data.table)
library(ChIPpeakAnno)
library(stringr)
library(gprofiler2)
library(GenomicRanges)
library(EnsDb.Hsapiens.v75)
library(biomaRt)
library(ComplexHeatmap)
library(clusterProfiler)
library(enrichplot)
library(GOSemSim)
library(gtools)
library(RColorBrewer)
library(circlize)

#------ plot z score heatmap ---------------------------------------------------

factor_loci <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/08_genes_and_pathways/factor_loci_moloc_closest_gene_info.txt', data.table = F) 

f_list <- list()
f_snps <- list()
for(i in 1:3){
  tt <- c('f1','f2', 'f3')[i]
  f_list[[tt]] <- fread(paste0('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/',tt,'_munged_build37.txt'), data.table = F)
  f_snps[[tt]] <- f_list[[tt]][f_list[[tt]]$SNP %in% factor_loci$SNP, ]
}

all_f_s<- do.call(rbind, f_snps)
all_f_s$z_score <- all_f_s$BETA/all_f_s$SE
to_plot <- data.frame(row.names =all_f_s[all_f_s$LHS=='F1', ]$SNP, f1= all_f_s[all_f_s$LHS=='F1', ]$z_score, f2=all_f_s[all_f_s$LHS=='F2', ]$z_score, f3=all_f_s[all_f_s$LHS=='F3', ]$z_score )

to_plot$type <- rep(NA, nrow(to_plot))
for(i in 1:nrow(to_plot)){
  to_plot[i,]$type <-paste0(unique(factor_loci[factor_loci$pan_locus==factor_loci[factor_loci$SNP==rownames(to_plot[i,]),]$pan_locus, ]$trait), collapse='-')
}



pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/heatmap_complete_absolute_z_scores_all_factors.pdf', height = 16, width = 8)
Heatmap(as.matrix(abs(to_plot[,1:3])),
        col = colorRamp2(c(seq(1,13, length.out=9)), c(RColorBrewer::brewer.pal(9,'Blues'))), 
        column_title = "Heatmap of absolute z scores" , 
        cluster_rows = T,
        cluster_columns = F,
        cluster_column_slices = T,
        cluster_row_slices =  F,
        show_row_dend = F,
        row_split = factor(to_plot[,4], levels = c('f1', 'f2', 'f3', 'f1-f2', 'f1-f3', 'f2-f3', 'f1-f2-f3' )),
        column_gap = unit(3, "mm"),
        row_gap = unit(5, "mm") ,
        border = T,
        row_names_gp = gpar( fontsize=0), 
        heatmap_legend_param =list(title = ""))
dev.off()


RColorBrewer::display.brewer.all(colorblindFriendly = T)



