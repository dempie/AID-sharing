plot.cojo.ht=function(cojo.ht.obj){
  require(ggplot2)
  library(patchwork)
  if(nrow(cojo.ht.obj$ind.snps)>1){
    
    whole.dataset=c()
    for(i in 1:nrow(cojo.ht.obj$ind.snps)){
      
      tmp=cojo.ht.obj$results[[i]]
      tmp$signal=cojo.ht.obj$ind.snps$SNP[i]
      whole.dataset=rbind(whole.dataset,tmp)
     
    }
    p1=ggplot(cojo.ht.obj$results[[i]],aes(x=bp,y=-log10(p)))+geom_smooth(method = 'loess',alpha=0.6,span = 0.05,color="black")+theme_minimal()
    p2=ggplot(whole.dataset,aes(x=bp,y=-log10(pC),color=signal))+geom_smooth(method = 'loess',alpha=0.6,span = 0.05)+
      theme_minimal()+ggtitle("Conditioned results")
    p3=p1/p2+ plot_layout(heights = c(1, nrow(cojo.ht.obj$ind.snps)+0.2))
    
  }else{
    
    
    p3=ggplot(cojo.ht.obj$results[[1]],aes(x=bp,y=-log10(p)))+geom_point(alpha=0.6)+theme_minimal()
    
    
    
  }
  
  print(p3)
  
  
  
}
