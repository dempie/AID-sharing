#  LD score regression on auotoimmunity GWAS
## graphing of V3 result 

#---- Libraries ----------------------------------------------------------------

library(corrplot)
library(qgraph)
library(dplyr)

#---- Load the dataset ---------------------------------------------------------

output2 <- readRDS('Outputs/Version3/LDS_output_final')

#look at the matrix

rownames(output2$S_Stand) <- colnames(output2$S_Stand) 
corrplot(output2$S_Stand, is.corr = F, addCoef.col = 'black', type = 'upper')


#- heritability --
cbind(colnames(output2$S_Stand), (diag(output2$S)) )

#---- remove columns we do not want --------------------------------------------
#remove asthma_1 (demeanis is the smallest one compared to asthma_2 that is han e al)
#remove Alzh_1 that is the smallest (kunkle), leave Alzh_2 that is the bigger one (Wightman)
#remove m1 that is interacting with everything because of the low h2
#remove armfat that is not needed, it's only my internal control fo heritability (around 0.20)
t1<- output2$S_Stand[ !rownames(output2$S_Stand) %in% c('ms_1', 'asthma_2', 'armfat', 'alzheimer_1') , !colnames(output2$S_Stand) %in% c('ms_1', 'asthma_2', 'armfat', 'alzheimer_1') ]




pre <- colnames(t1)

colnames(t1) <-c('Allergies',  "Multiple sclerosis ", 'Alzheimer\'s disease', 'Asthma', 'Celiac disease',
                 'Crohn\'s disease', 'Juvenile idiopathic arthritis', 'Primary biliary cholangitis', 'Primary sclerosing cholangitis',
                 'Rheumatoid arthritis' , 'Systemic lupus erythematosus', 'Autoimmune thyroid disease', 'Ulcerative colitis')

rownames(t1) <- colnames(t1)

after <- colnames(t1)
rbind(pre, after)

blue_red <- colorRampPalette(c('navy', 'white', 'red3')) 
corrplot(t1, order = 'hclust', addCoef.col = 'black', 
         type = 'upper', tl.col = 'black', col = blue_red(200),
         number.cex=0.75, 
         tl.cex= 0.8)


#export as PDF
pdf('Graphs/Correlation_matrix_V3.pdf', height = 14, width = 14)
corrplot(t1, order = 'hclust', addCoef.col = 'black', 
         type = 'upper', tl.col = 'black', col = blue_red(200),
         number.cex=1, 
         tl.cex= 1.2,
         tl.srt = 45)
dev.off()

#without correlation coefficients
pdf('Graphs/Correlation_matrix_V3_nocoeff.pdf', height = 14, width = 14)
corrplot(t1, order = 'hclust', 
         type = 'upper', tl.col = 'black', col = blue_red(200),
         number.cex=1, 
         tl.cex= 1.2,
         tl.srt = 45)
dev.off()





#---- remove autocorrealtion for qraph

t_no_auto <- t1
diag(t_no_auto) <- rep(0, 15)
qgraph(t_no_auto,threshold=0.4,layout="spring")

