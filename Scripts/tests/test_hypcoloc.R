snp_23<- rownames(betas$locus_23)
f1 <- list_of_files[[4]]
f2 <- list_of_files[[5]]
f3 <- list_of_files[[6]]

f1_23 <- f1[f1$SNP %in% snp_23, c('SNP', 'P')]
colnames(f1_23) <- c('rsid', 'pval')
f2_23 <- f2[f2$SNP %in% snp_23, c('SNP', 'P')]
colnames(f2_23) <- c('rsid', 'pval')
f3_23 <- f3[f3$SNP %in% snp_23, c('SNP', 'P')]
colnames(f3_23) <- c('rsid', 'pval')

a <- 'f1_23.csv'
b <- 'f3_23.csv'

a <- locuscompare(in_fn1 = 'f1_233.tsv', in_fn2 = 'f3_233.tsv' , marker_col1 =  'rsID', pval_col1 = 'pval',combine = F,title1 = 'F1', title2 = 'F3',  )
b <- locuscompare(in_fn1 = 'f3_233.tsv', in_fn2 = 'f3_233.tsv' , marker_col1 =  'rsID', pval_col1 = 'pval',combine = F,title1 = 'F3', title2 = 'F3' )
k <- locuscompare(in_fn1 = 'f2_233.tsv', in_fn2 = 'f2_233.tsv' , marker_col1 =  'rsID', pval_col1 = 'pval',combine = F,title1 = 'F2', title2 = 'F2' )

colnames(f1_23)

fwrite(f1_23, 'f1_233.tsv',  sep = '\t', col.names = T, row.names = F, quote = F)
fwrite(f3_23, 'f3_233.tsv',  sep = '\t', col.names = T, row.names = F, quote = F)
write.table(f1_23, 'f1_233.tsv', quote = F , row.names = F)
write.table(f3_23, 'f3_233.tsv', quote = F , row.names = F)
write.table(f2_23, 'f2_233.tsv', quote = F , row.names = F)
rownames(f2_23) <- NULL

head(f2_23)




in_fn1 = system.file('extdata','gwas.tsv', package = 'locuscomparer')
in_fn2 = system.file('extdata','eqtl.tsv', package = 'locuscomparer')
locuscompare(in_fn1 = in_fn1, in_fn2 = in_fn2)


layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
a
k
b
a <- a$locuszoom1

k<- k$locuszoom1

b <- b$locuszoom2



library(ggplot2)
library(ggpubr)

figure <- ggarrange(a, k, b,
                    labels = c("F1", "F2", "F3"),
                    ncol = 1, nrow = 3)





install_github("whweve/IntAssoPlot")












pi <- locuscompare(in_fn1 = 'f1_233.tsv', in_fn2 = 'f3_233.tsv' , marker_col1 =  'rsID', pval_col1 = 'pval',title1 = 'F1', title2 = 'F3', combine = T, snp =  )
po <- locuscompare(in_fn1 = 'f1_233.tsv', in_fn2 = 'f2_233.tsv' , marker_col1 =  'rsID', pval_col1 = 'pval',title1 = 'F3', title2 = 'F3' , combine = F)
pu <- locuscompare(in_fn1 = 'f2_233.tsv', in_fn2 = 'f3_233.tsv' , marker_col1 =  'rsID', pval_col1 = 'pval',title1 = 'F2', title2 = 'F2' , combine = F)

pi$locuscompare
po$locuscompare
pu$locuscompare





