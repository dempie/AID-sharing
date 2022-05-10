#in this scirpt I will perform downstrean analysis of the GWAS results.

#----locus.breaker function by Nicola, it identifies loci-----------------------
locus.breaker=function(res,p.sig=5e-8, p.limit=1e-5,hole.size=250000
                       ,p.label="p",chr.label="chr",pos.label="pos"){
  
  res=res[which(res[,p.label]<p.limit),]
  trait.res=c()
  for(j in 1:22){
    
    res.chr=res[which(res[,chr.label]==j),]
    if(nrow(res.chr)>1){
      holes=res.chr[,pos.label][-1]-res.chr[,pos.label][-length(res.chr[,pos.label])] 
      gaps=which(holes>hole.size)
      if(length(gaps)>0){
        for(k in 1:(length(gaps)+1)){
          
          if(k==1){
            res.loc=res.chr[1:(gaps[k]),]  
          }else if(k==(length(gaps)+1)){
            res.loc=res.chr[(gaps[k-1]+1):nrow(res.chr),]  
          }else{
            res.loc=res.chr[(gaps[k-1]+1):(gaps[k]),]
          }
          if(min(res.loc[,p.label])<p.sig){
            
            start.pos=min(res.loc[,pos.label],na.rm=T)
            end.pos=max(res.loc[,pos.label],na.rm=T)
            chr=j
            best.snp=res.loc[which.min(res.loc[,p.label]),]
            line.res=c(chr,start.pos,end.pos,unlist(best.snp))
            trait.res=rbind(trait.res,line.res)
          }
          
          
        }
      }else{
        res.loc=res.chr
        if(min(res.loc[,p.label])<p.sig)  {
          
          start.pos=min(res.loc[,pos.label],na.rm=T)
          end.pos=max(res.loc[,pos.label],na.rm=T)
          chr=j
          best.snp=res.loc[which.min(res.loc[,p.label]),]
          line.res=c(chr,start.pos,end.pos,unlist(best.snp))
          trait.res=rbind(trait.res,line.res)
        }
        
      }
      
    }else if(nrow(res.chr)==1){
      
      res.loc=res.chr
      if(min(res.loc[,p.label])<p.sig){
        start.pos=min(res.loc[,pos.label],na.rm=T)
        end.pos=max(res.loc[,pos.label],na.rm=T)
        chr=j
        best.snp=res.loc[which.min(res.loc[,p.label]),]
        line.res=c(chr,start.pos,end.pos,unlist(best.snp))
        trait.res=rbind(trait.res,line.res)
      }
      
      
    }
  }
  
  print(trait.res)
  trait.res=as.data.frame(trait.res,stringsAsFactors=FALSE)
  trait.res=trait.res[,-(which(names(trait.res)==chr.label))]
  names(trait.res)[1:3]=c("chr","start","end")
  trait.res
}
#------ run locus breaker function----------------------------------------------

#load the complete summary stats
f1_gwas <- readRDS('outputs/version3/05_GWAS_results/F1_completeGWAS.RDS')
f2_gwas <- readRDS('outputs/version3/05_GWAS_results/F2_completeGWAS.RDS')

nrow(f1_gwas) #4309927
nrow(f2_gwas) #4309927

#remove NAs
f1_gwas_ok <- f1_gwas[which(!is.na(f1_gwas$Pval_Estimate)), ]
f2_gwas_ok <- f2_gwas[which(!is.na(f2_gwas$Pval_Estimate)), ]

nrow(f1_gwas_ok) #4281353
nrow(f2_gwas_ok) #4281353

#run the function

f1_lb_output <- locus.breaker(res= f1_gwas_ok, p.label = <-  'Pval_Estimate', 
                              chr.label = 'CHR', pos.label = 'BP'  )



f2_lb_output <- locus.breaker(res= f2_gwas_ok, p.label = 'Pval_Estimate', 
                              chr.label = 'CHR', pos.label = 'BP'  )


dim(f1_lb_output)
dim(f2_lb_output)


a<- f1_lb_output$SNP

b <- f2_lb_output$SNP






