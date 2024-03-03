##### The script template for plotting LOLA enrichment result in Figure 2 and Figure S2 #####
setwd("/home/owht/Data/KIZ/data/Entropy2_rev/")
library(LOLA)
library(ggplot2)
library(stringr)
load("LOLA_DER_cen_low.RData")

Result<-locDER_histone[,c("userSet","description","collection","filename",
                          "oddsRatio","qValue")]
Result
Result$cell<-Result$collection

Result$OR<-"<1"
Result[(Result$oddsRatio > 1) & (Result$oddsRatio<60),]$OR<-"1-60"
#Result[(Result$oddsRatio > 5) & (Result$oddsRatio<60),]$OR<-"5-60"
Result[(Result$oddsRatio>60),]$OR<-">60"
Result$OR<-factor(Result$OR,levels = c("<1","1-60",">60"))

Result$q.value<-"n.s."
Result[(Result$qValue < 0.05),]$q.value<-"<0.05"
Result[(Result$qValue < 1e-10),]$q.value<-"<1e-10"
Result[(Result$qValue < 1e-50),]$q.value<-"<1e-50"
Result[(Result$qValue < 1e-100),]$q.value<-"<1e-100"
Result[(Result$qValue < 1e-300),]$q.value<-"<1e-300"
Result$q.value<-factor(Result$q.value,
                       levels = c("n.s.","<0.05","<1e-10","<1e-50","<1e-100","<1e-300"))
library(ggplot2)
library(stringr)
Result$histone<-str_split(Result$description,pattern = "_",simplify = T)[,2]
Result$cells<-str_split(Result$cell,pattern = "_",simplify = T)[,1]
#levels(Result$filename)<-
ggplot(Result,aes(alpha=q.value,x=cell,y=histone))+
  geom_point(aes(size=oddsRatio,fill=cell),shape=21,color="black")+theme_bw()

#################################
load("/home/owht/Data/KIZ/data/Entropy2_rev/Meth/LOLA_DMC_133.RData")
Result<-chromHMM_lola[,c("userSet","description","collection","filename",
                           "oddsRatio","qValue")]

Result
Result$OR<-"<1"
Result[(Result$oddsRatio > 1) & (Result$oddsRatio<1.5),]$OR<-"1-1.5"
Result[(Result$oddsRatio > 1.5) & (Result$oddsRatio<5),]$OR<-">1.5"
#Result[(Result$oddsRatio > 5) ,]$OR<-">5"
Result$OR<-factor(Result$OR,levels = c("<1","1-1.5",">1.5"))
Result$q.value<-"n.s."
Result[(Result$qValue < 0.05),]$q.value<-"<0.05"
Result[(Result$qValue < 0.01),]$q.value<-"<0.01"
Result[(Result$qValue < 0.001),]$q.value<-"<0.001"
Result[(Result$qValue < 1e-5),]$q.value<-"<1e-5"
Result[(Result$qValue < 1e-10),]$q.value<-"<1e-10"
Result[(Result$qValue < 1e-50),]$q.value<-"<1e-100"
Result$q.value<-factor(Result$q.value,
                       levels = c("n.s.","<0.05","<0.01","<0.001",
                                  "<1e-5","<1e-10","<1e-100"))

library(ggplot2)
Result
Result$filename<-factor(Result$filename,
                        levels = c("1_Active_Promoters.bed",
                                   "2_Weak_Promoters.bed",
                                   "3_Poised_Promoters.bed",
                                   "4_Strong_Enhancer.bed",
                                   "5_Strong_Enhancer.bed",
                                   "6_Weak_Enhancer.bed",
                                   "7_Weak_Enhancer.bed",
                                   "8_Insulator.bed",
                                   "9_Txn_Transition.bed",
                                   "10_Txn_Elongation.bed",
                                   "11_Weak_Txn.bed",
                                   "12_Repressed.bed",
                                   "13_Heterochrom.bed",
                                   "14_Repetitive.bed",
                                   "15_Repetitive.bed"))

Result

ggplot(Result,aes(fill=q.value,alpha=q.value,x=collection,y=filename))+
  geom_point(aes(size=OR))

ggplot(Result[Result$collection == "GM12878",],
       aes(alpha=q.value,x=collection,y=filename))+
  geom_point(aes(size=OR),color="red")+theme_bw()+
  theme(text = element_text(size=14))

