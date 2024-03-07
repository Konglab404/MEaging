setwd("~/data/Entropy/Result/")

Entropy_SE<-readRDS("Ent_SE_133_0.5_for3cats_rb.RDS")
sampleinfo<-as.data.frame(colData(Entropy_SE))
colnames(meth133)[2:134]<-sampleinfo$ID

meth133<-readRDS("Meth_133samples.aggregate.RDS")
load("DER_CEN_20230701.RData")

meth133<-meth133[meth133$loc %in% Entropy_SE@rowRanges$loc,]
DER_meth<-meth133[meth133$loc %in% DER_CEN$loc,]
sampleinfo<-read.table("../Sampleinfo133")

CEN<-sampleinfo[sampleinfo$group == "Long-lived",]$ID
YC<-sampleinfo[sampleinfo$group == "Younger Control",]$ID
EC<-sampleinfo[sampleinfo$group == "Elder Control",]$ID

CEN_mean<-rowMeans(meth133[,CEN,with=F],na.rm = T)
YC_mean<-rowMeans(meth133[,YC,with=F],na.rm = T)
EC_mean<-rowMeans(meth133[,EC,with=F],na.rm = T)

# CEN_low_mean<-rowMeans(DER_meth[,CEN,with=F],na.rm = T)
# EC_low_mean<-rowMeans(DER_meth[,EC,with=F],na.rm = T)
# YC_low_mean<-rowMeans(DER_meth[,YC,with=F],na.rm = T)

CEN_low_mean<-colMeans(DER_meth[,CEN,with=F],na.rm = T)
EC_low_mean<-colMeans(DER_meth[,EC,with=F],na.rm = T)
YC_low_mean<-colMeans(DER_meth[,YC,with=F],na.rm = T)

boxplot(YC_low_mean,EC_low_mean,CEN_low_mean)

ggplot(data.frame(meth=c(CEN_mean,YC_mean,EC_mean),
                  type=rep(c("CEN","YC","EC"),each=length(CEN_mean))),
       aes(x=meth,color=type))+
  geom_density()+theme_classic()
