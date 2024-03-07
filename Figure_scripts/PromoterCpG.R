library(BSgenome.Hsapiens.UCSC.hg19)
annotations = build_annotations(genome = 'hg19', annotations = annos)
anno2<-annotations[!is.na(annotations$gene_id),]
pros<-getSeq(BSgenome.Hsapiens.UCSC.hg19,anno2)
pros<-anno2[anno2$type == "hg19_gene_promoters",]
pro_seq<-getSeq(BSgenome.Hsapiens.UCSC.hg19,pros)


for(i in 1:length(pro_seq)){
  temp<-as.character(pro_seq[i])
  promoter[i,]$GeneID<-temp %>% names %>% as.numeric()
  numC<-str_count(temp,"C")
  numG<-str_count(temp,"G")
  numCG<-str_count(temp,"CG")+str_count(temp,"GC")
  c_frq<-numC/1000
  g_frq<-numG/1000
  cg_frq<-numCG/1000
  promoter[i,]$c_frq<-c_frq; promoter[i,]$g_frq<-g_frq; promoter[i,]$cg_frq<-cg_frq
  promoter[i,]$CpG_dens<-numCG/20; promoter[i,]$CpGOE<-cg_frq / (c_frq*g_frq)

}

promoter[promoter$CpGOE>3,]$CpGOE<-3
quantile(promoter$CpGOE,c(0.33,0.66,1))
promoter[promoter$CpGOE>1.643554 ,]$type<-"High CpG Promoter"
promoter$type<-"NO"
promoter[promoter$CpGOE>1.643554 ,]$type<-"High CpG Promoter"
promoter[promoter$CpGOE<1.304713,]$type<-"Low CpG Promoter"
promoter[promoter$type == "NO",]$type<-"Intermediate CpG Promoter"

save(promoter,file ="promoter_CpGOE_annotatr_hg19.RData")
