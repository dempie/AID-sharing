library(data.table)
library(stringr)
library(ComplexHeatmap)
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

factor_loci$pan_locus <- paste0('locus', factor_loci$pan_locus)

all_f_s<- do.call(rbind, f_snps)
all_f_s$z_score <- all_f_s$BETA/all_f_s$SE
to_plot <- data.frame(row.names = unique(factor_loci$pan_locus_name), 
                      F1=rep(NA, length(unique(factor_loci$pan_locus_name))),  
                      F2=rep(NA, length(unique(factor_loci$pan_locus_name))),
                      F3=rep(NA, length(unique(factor_loci$pan_locus_name))), 
                      traits=rep(NA, length(unique(factor_loci$pan_locus_name))),
                      leads=rep(NA, length(unique(factor_loci$pan_locus_name))))



for(i in 1:nrow(to_plot)){
  tt <- rownames(to_plot)[i]
  to_plot[tt, ]$traits <- paste0(factor_loci[factor_loci$pan_locus_name== tt,]$trait, collapse = '-')
  to_plot[tt, ]$leads <- paste0(factor_loci[factor_loci$pan_locus_name== tt, ]$SNP, collapse = '-')
  
}


mean(to_plot[to_plot$traits %in% c('f1', 'f2', 'f3'),]$leads %in% all_f_s$SNP) #all SNP are present in the factor specifics group


#generate a list of factor specific zetas
to_plot_f <- to_plot[to_plot$traits %in% c('f1', 'f2', 'f3'),]
for(i in 1:3){
      tt <- c('f1','f2', 'f3')[i]
     
      for(k in 1:nrow(to_plot_f)){
            snp <- to_plot_f[k,]$leads
            to_plot_f[to_plot_f$leads==snp, toupper(tt)]  <-   all_f_s[all_f_s$LHS==toupper(tt) & all_f_s$SNP==snp,]$z_score
      }
  
}


#----------- z scores only factor specific -------------------------------------
RColorBrewer::display.brewer.all(colorblindFriendly = T)


pdf('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure1/heatmap_complete_absolute_z_scores_specifc_factors.pdf', height = 16, width = 6)
Heatmap(as.matrix(abs(to_plot_f[,1:3])),
        col = colorRamp2(c(seq(1,13, length.out=9)), c(RColorBrewer::brewer.pal(9,'Purples'))), 
        column_title = "Heatmap of absolute z scores" , 
        cluster_rows = T,
        cluster_columns = F,
        cluster_column_slices = F,
        cluster_row_slices =  F,
        show_row_dend = F,
        row_split = factor(toupper(to_plot_f[,'traits']), levels = c('F1', 'F2', 'F3')),
        column_gap = unit(3, "mm"),
        row_gap = unit(5, "mm") ,
        border = T,
        row_names_gp = gpar( fontsize=0), 
       # row_title='Lead SNP',
        heatmap_legend_param =list(title = "Absosule z score"))
dev.off()






















