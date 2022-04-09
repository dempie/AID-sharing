#  LD score regression on auotoimmunity GWAS
## graphing of V3 result 

#---- Libraries ----------------------------------------------------------------

library(corrplot)
library(qgraph)
library(dplyr)

#---- Load the dataset ---------------------------------------------------------

output2 <- readRDS('Outputs/Version3/LDS_output_final')


#---- remove columns we do not want --------------------------------------------

#remove m1 that is interacting with everything because of the low h2
t1<- output2$S_Stand[-11, -11]

#remove asthma_1 (demeanis is the smallest one compared to asthma_2 that is han e al)

t1<- output2$S_Stand[-11, -11]
#remove Alzh_1 that is the smallest (kunkle), leave Alzh_2 that is the bigger one (Wightman)

t1<- output2$S_Stand[-11, -11]

colnames(output2$S_Stand)
colnames(output2$S_Stand) <-c('Multiple sclerosis',  "Alzheimer\'s disease ", 'Arm fat percentage', ''  )

corrplot(t1, order = 'hclust', addCoef.col = 'black', type = 'upper', tl.col = 'black')


#---- remove autocorrealtion for qraph

t_no_auto <- t1
diag(t_no_auto) <- rep(0, 15)
qgraph(t_no_auto,threshold=0.4,layout="spring")

