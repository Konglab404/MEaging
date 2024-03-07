setwd("~/data/Entropy/Result/Summary/")
library(BSgenome.Hsapiens.UCSC.hg19)
annos<-c('hg19_basicgenes')
annotations = build_annotations(genome = 'hg19', annotations = annos)
anno2<-annotations[!is.na(annotations$gene_id),]
pros<-anno2[anno2$type == "hg19_genes_promoters",]
pro_seq<-getSeq(BSgenome.Hsapiens.UCSC.hg19,pros)

promoter<-data.frame(id = pros$id)
#promoter<-promoter[,-6]
promoter$GeneID<-rep("-1",nrow(promoter))
promoter$CpG_dens<-rep(-1,nrow(promoter))
promoter$CpGOE<-rep(-1,nrow(promoter))
promoter$c_frq<- -1
promoter$g_frq<- -1
promoter$cg_frq <- -1
promoter$max_cpgoe<- -1

slide_cpgoe<-function(index){
  xulie<-str_sub(temp,index,index + 499)
  cg_frq<-str_count(xulie,"CG")/500+str_count(xulie,"GC")/500
  c_frq<-str_count(xulie,"C")/500
  g_frq<-str_count(xulie,"G")/500
  cpgoe<-cg_frq / (c_frq*g_frq)
  return(cpgoe)
}


for(i in 1:length(pro_seq)){
  temp<-as.character(pro_seq[i])
  promoter[i,]$GeneID<-pros$gene_id[i]
  numC<-str_count(temp,"C")
  numG<-str_count(temp,"G")
  numCG<-str_count(temp,"CG")+str_count(temp,"GC")
  c_frq<-numC/1000
  g_frq<-numG/1000
  cg_frq<-numCG/1000
  promoter[i,]$c_frq<-c_frq; promoter[i,]$g_frq<-g_frq; promoter[i,]$cg_frq<-cg_frq
  promoter[i,]$CpG_dens<-numCG/20; promoter[i,]$CpGOE<-cg_frq / (c_frq*g_frq)
  #promoter[i,]$max_cpgoe<-sapply(seq(1,2000,5),slide_cpgoe) %>% na.omit() %>% max
}

promoter<-na.omit(promoter)
promoter[promoter$CpGOE>3,]$CpGOE<-3
hist(promoter$CpGOE,breaks = 100)

quantile(promoter$CpGOE,c(0.33,0.66,1))
# 33%      66%     100% 
# 1.304713 1.643554 3.000000 
promoter$type<-"NO"

promoter[promoter$CpGOE>1.643554 ,]$type<-"High CpG Promoter"
promoter[promoter$CpGOE<1.304713,]$type<-"Low CpG Promoter"
promoter[promoter$type == "NO",]$type<-"Intermediate CpG Promoter"

save(promoter,file ="promoter_CpGOE_annotatr_hg19.RData")
