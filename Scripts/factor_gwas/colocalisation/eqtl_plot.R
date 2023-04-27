
library(data.table)
library(lemon)
library(ggplot2)
file.list=system("ls output/sc_eqtl/*MR.tsv",intern=T)

all.res=c()

for(i in file.list){
  
  
  tmp=try(fread(i))
  if(!is(tmp,"try-error")){
  
    
    all.res=rbind(all.res,tmp)
    
    
  }
  
  
}

all.res=all.res[which(all.res$pval<(0.05/120)),]

all.res$direction=ifelse(all.res$b>0,1,-1)
pdf("output/sc_eqtl/eqtl_MR_all.pdf",height=60,width=15)
ggplot(all.res,aes(y=trait,x=paste(gene,locus,sep=" "),shape=as.factor(direction),fill=direction))+
  geom_point(size=3)+scale_shape_manual(values = c(25,24))+
  scale_fill_gradient2()+theme_minimal()+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5),axis.line=element_line(color="black"),legend.position="none")+
  facet_rep_grid(cell.type~.,repeat.tick.labels=TRUE)
dev.off()



all.mr.f=all.res[all.res$trait%in%c("f1","f2","f3"),]
pdf("output/sc_eqtl/eqtl_MR_factors.pdf",height=10,width=10)
ggplot(all.mr.f,aes(y=cell.type,x=paste(gene,locus,sep=" "),shape=factor(direction,levels=c("-1","1")),fill=factor(trait,levels=c("f1","f2","f3"))))+
  geom_point(size=4)+scale_shape_manual(values = c(25,24))+
  scale_fill_manual(values = c("#1F78B4","#33A02C","#E31A1C"))+theme_minimal()+
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5),legend.position="none")+coord_flip()
dev.off()
#shapes 24 25



file.list=system("ls output/sc_eqtl/*coloc.tsv",intern=T)
all.res=c()

for(i in file.list){
  
  
  tmp=fread(i)
  if(ncol(tmp)>1){
    
    
    all.res=rbind(all.res,tmp)
    
    
  }
  
  
}



