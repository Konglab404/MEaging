library(BSgenome.Hsapiens.UCSC.hg19)
chrs <- names(Hsapiens)[1:24]
cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
cpgr <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 2))))

for(i in 1:22){
  tmp<-start(cpgr[cpgr@seqnames == paste0("chr",i),])
  tmp2<-tmp+1
  tmp3<-rep(0,length(tmp)*2)
  tmp3[seq(1,length(tmp3)-1,2)]<-tmp
  tmp3[seq(2,length(tmp3),2)]<-tmp2
  write.table(tmp3,file=paste0(i,".gp"),quote=F,
              col.names = F,row.names = F)
}
