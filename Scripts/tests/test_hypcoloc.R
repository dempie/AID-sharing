#test hyperloc 

my_sumstats <- readRDS('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/gwas_final_withQindex.RDS')

f1_gwas <- my_sumstats$factor1

uc <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/uc_delange-2017.txt', data.table = F)
cd <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/crohn_delange-2017.txt', data.table = F)
psc <- fread('outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/01_qc_sumstats/ready_for_munge/psc_ji-2016.txt', data.table = F)


jpeg(file='outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/miamiF1_F2_F3.jpeg', width = 800, height = 400, units='mm', res = 600)
par(mfrow=c(3,1))
par(mar=c(5,5,3,3))
manhattan(F1_plot, chr="CHR", bp="BP", snp="SNP", p="Pval_Estimate" ,ylim=c(0,30), main='F1 =~ crohn + uc  + psc', 
          col = c("darksalmon", "darkseagreen4"))  
par(mar=c(5,5,3,3))
manhattan(F2_plot, chr="CHR", bp="BP", snp="SNP", p="Pval_Estimate" ,ylim=c(0,30), main='F2 =~ jia + sle + ra + t1d', 
          col = c("darksalmon", "darkseagreen4"))
par(mar=c(5,5,3,3))
manhattan(F3_plot, chr="CHR", bp="BP", snp="SNP", p="Pval_Estimate" ,ylim=c(0,30), main='F3 =~ asthma + derma', 
          col = c("darksalmon", "darkseagreen4"))
dev.off()