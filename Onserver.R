setwd("data/HNLong/DMEAS/Another_81/Int3CG")
load("tmp2.RData")
library(limma)
library(magrittr)

Total<-cbind(CENsm,SPsm) %>% as.matrix

#limma pipeline #

#Design by Cell compositions
#CEN - F1SP


load("Methvalue_89_3CG.RData")
cc<-read.csv("Cell_composition.csv",row.names=1)
rownames(cc)<-paste0("Entropy",rownames(cc))
cc<-cc[colnames(Methvaluesm),]

load("../../../AllSampleInfo.RData") 
sampleinfo$sub<-paste0("Entropy",sampleinfo$sub)
cc<-cc[sampleinfo$sub,] %>% na.omit
cc$type<-sampleinfo$type

design<-model.matrix(~type+Lymph.+Mid.+Gran.,data=cc)

library(limma)
fit<-lmFit(Methvaluesm,design)
fit2<-eBayes(fit)

result_meth<-topTable(fit2,coef=2,adjust="BH",number=Inf)

