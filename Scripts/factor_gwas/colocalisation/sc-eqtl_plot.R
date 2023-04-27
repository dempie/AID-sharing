### collect results
library(data.table)
file.list=system("ls output/sc_eqtl/Chr*",intern = T)

all.colocs=c()
for(i in file.list){
  
  chr=gsub("output/sc_eqtl/Chr","",i)|>strsplit(split="_")|>unlist()
  tmp=fread(i,data.table=F)
  tmp$chr=chr[1]
  if(ncol(tmp)>2){
    colo=tmp[which(tmp$PP.H4.abf>0.9),]
    if(nrow(colo)>0){
    
      all.colocs=rbind(all.colocs,tmp)
    
    }
  }
}

all.colocs$pan_locus=as.numeric(apply(t(all.colocs$locus),2,function(x)unlist(strsplit(x,split="_"))[1]))
loci.table=fread("output/loci_definitions/final_locus_table.tsv")
loci.table$final.locus=paste0(loci.table$Chr,":",loci.table$start,"-",loci.table$end,"_",loci.table$sub_locus)
loci.table=loci.table[order(loci.table$sub_locus),]
loci.table=loci.table[order(loci.table$pan.locus),]
all.colocs$locus=loci.table$final.locus[match(all.colocs$hit1,loci.table$SNP)]

all.colocs=all.colocs[order(all.colocs$pan_locus),]
all.colocs=unique(all.colocs)



write.table(all.colocs,file="output/sc_eqtl/All_colocalization_eqtl.tsv",row.names=F,quote=F,sep="\t")


res=fread("output/sc_eqtl/All_colocalization_eqtl.tsv")
f.t=res[res$PP.H4.abf>0.9 & res$t1 %in%c("f1","f2","f3"),]

order=table(f.t$t1,f.t$t2)
order[order>1]=1
order=(hclust(dist(t(order)),'ward.D2'))
order=order$labels[order$order]
a=table(f.t$t2,f.t$cell_type,f.t$t1)
a=a[order,,]
a[a>1]=1
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}

trans=t_col("red",percent = 100)

f1.c=t_col("#1F78B4",percent = 20)
f2.c=t_col("#33A02C",percent = 20)
f3.c=t_col("#E31A1C",percent = 20)
f1.f2=t(a[,,"f1"])+t(a[,,"f2"])
f1.f2[f1.f2==1]=0
f1.f2[f1.f2==2]=1
f1.f2.f3=t(a[,,"f1"])+t(a[,,"f3"])+t(a[,,"f3"])
f1.f2.f3[f1.f2.f3<=2]=0
f1.f2.f3[f1.f2.f3==3]=1


pdf("Coloc_eqtl.pdf",width=21)
corrplot(t(a[,,"f1"]),col=c("red",trans,f1.c),cl.pos='n',bg=NA,tl.col="white",method = "square")
corrplot(t(a[,,"f2"]),col=c("red",trans,f2.c),cl.pos='n',add = T,bg=NA,tl.col="white",method = "square")
corrplot(t(a[,,"f3"]),col=c("red",trans,f3.c),cl.pos='n',add = T,,bg=NA,tl.col="white",method = "square")
corrplot(f1.f2,col=c("red",trans,"purple"),cl.pos='n',add = T,,bg=NA,tl.col="white",method = "square")
corrplot(f1.f2.f3,col=c("red",trans,"black"),cl.pos='n',add = T,,bg=NA,tl.col="black",method = "square")
dev.off()


table()






