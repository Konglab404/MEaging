setwd("/home/wanght/data/Entropy_tmp")
Entropy<-readRDS("Entropy_133samples.aggregate.RDS")
load("meth133_range.RData")
library(stringr)
methRanges<-str_split(Entropy$loc,pattern = ":|\\_",simplify = T)
methRanges<-methRanges[,c(1,2,4)];colnames(methRanges)<-c("seqnames","start","end")
#methRange<-as.data.frame(methRanges)
methRange<-GRanges(seqnames=methRanges[,1],
                   IRanges(
                     as.numeric(methRanges[,2]),
                     as.numeric(methRanges[,3])
                   ),strand="*")

##### TE #####
TE<-readRDS("./TE/RepeatMask_hg19_UCSC_repclass.RDS")
TE<-TE[TE$repClass %in% c("LINE","SINE","DNA","LTR"),]
TvE<-findOverlaps(methRange,TE) %>% as.data.frame() #Query is entropy loci
TvE$class<-TE[TvE$subjectHits,]$repClass

Results<-data.frame(ID = colnames(Entropy)[2:134],
                    SINE=0, LINE=0, DNA=0, LTR=0,
                    Enhancer=0, Insulator=0)
for(i in c("LINE","SINE","DNA","LTR")){
  Results[,i]<-apply(Entropy[TvE[TvE$class == i,]$queryHits,2:134],2,mean,na.rm=T)
}

## Enh/Ins #####

load("chromHMM/chromHMM_gr.RData")
chromHMMgr<-na.omit(chromHMMgr)
EvE<-findOverlaps(methRange,chromHMMgr) %>% as.data.frame()
EvE$class<-chromHMMgr[EvE$subjectHits,]$Region
Results[,"Enhancer"]<-apply(Entropy[EvE[grep(EvE$class,pattern="nhancer"),]$queryHits,2:134],2,mean,na.rm=T)
Results[,"Insulator"]<-apply(Entropy[EvE[grep(EvE$class,pattern="sulator"),]$queryHits,2:134],2,mean,na.rm=T)

save(Results,file="TE_enh_insu.RData")
