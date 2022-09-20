
library(corrplot)

#-----
ldsc_ok <- readRDS('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/02_sumstats_function/ldsc_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.RDS')
names_plot = c('Type 1 Diabetes','Crohn\'s disease', 'Ulcerative colitis', 
               'Primary sclerosing cholangitis ', 'Juvenile idiopathic arthritis', 'Systemic lupus erythematosus', 
               'Rheumatoid arthritis',  'Asthma', 'Atopic dermatitis')



rownames(ldsc_ok$S_Stand) <- names_plot
colnames(ldsc_ok$S_Stand) <- names_plot

plot <- ldsc_ok$S_Stand[c('Crohn\'s disease', 'Ulcerative colitis','Primary sclerosing cholangitis ', 'Juvenile idiopathic arthritis', 'Systemic lupus erythematosus','Rheumatoid arthritis','Type 1 Diabetes', 'Asthma' , 'Atopic dermatitis'),c('Crohn\'s disease', 'Ulcerative colitis','Primary sclerosing cholangitis ', 'Juvenile idiopathic arthritis', 'Systemic lupus erythematosus','Rheumatoid arthritis', 'Type 1 Diabetes','Asthma' , 'Atopic dermatitis')]
pdf('outputs/2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/plots/figure1/correlation_matrix_uc-cd-psc-ra-sle-t1d-jia-asthma-derma.pdf', height = 14, width = 14)
corrplot(plot, 
         order = 'original',
         addCoef.col = 'black', 
         method = 'square',
         type = 'upper', 
         is.corr = T,
         tl.col='black',
         outline=F,
         number.cex= 1.5,
         tl.cex=1.5,
         cl.cex=1.5,
         tl.srt= 45,
         addgrid.col='grey', 
         title= 'LDSC genetic correlation matrix of autoimmune diseases', mar=c(0,0,1,0)
)

dev.off()
