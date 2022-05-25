
#in jia I noticed that the beta they reported is not exatcly OR when the log(OR) are compared to the beta. 
#in this script I comapred them and theier z scores. 


library(corrplot)
jia <- fread('Summary_Stats/lopezisac-2020_jia_build37_GCST90010715_buildGRCh37.tsv', data.table = F)


cor(log(jia$all_OR), jia$frequentist_add_beta_1) #0.9719502, really good correlation 


rand <- sample(1:nrow(jia), 100000) #take a random set of them to plot 

#taking the log of odds ratio
a <- lm(jia$frequentist_add_beta_1~log(jia$all_OR))
summary(a)

plot(log(jia$all_OR[rand]), jia$frequentist_add_beta_1[rand], xlim = c(-1, +1),ylim=c(-1.5, +1.5), main =paste0('beta ~ log(OR)', '\n','R-squared= 0.9447','\n', 'regression estimate = 1.111')  ) 
abline(lm(jia$frequentist_add_beta_1~log(jia$all_OR)), col="red")# regression line (y~x)
abline(0, 1, col='blue', lty= 4 )



pdf(width = 12, height = 7, 'jia_beta_or.pdf')
#par(mfrow=c(1,2))
plot(log(jia$all_OR[rand]), jia$frequentist_add_beta_1[rand], xlim = c(-1, +1.5),ylim=c(-1, +1.5), main =paste0('beta ~ log(OR)', '\n','R-squared= 0.9447','\n', 'regression estimate = 1.111'), xlab = 'log(OR)', ylab = 'beta'  ) 
abline(lm(jia$frequentist_add_beta_1~log(jia$all_OR)), col="red")# regression line (y~x)
abline(0, 1, col='blue', lty= 4 )
dev.off()




#----z scores

z_jia_or <- log(jia$all_OR) /((log(jia$all_OR_upper) -  log(jia$all_OR))/qnorm(0.975)) 
z_jia_beta <- jia$frequentist_add_beta_1 / jia$frequentist_add_se_1

a<- lm(z_jia_or~ z_jia_beta)

summary(a)

plot(z_jia_beta[rand] , z_jia_or[rand] , main =paste0('z OR ~ z beta', '\n','R-squared= 0.9676 ','\n', 'regression estimate = 0.9614'), xlab = 'z beta', ylab = 'z OR'  ) 
abline(lm(z_jia_or~ z_jia_beta), col="red")# regression line (y~x)
abline(0, 1, col='blue', lty= 4 )












