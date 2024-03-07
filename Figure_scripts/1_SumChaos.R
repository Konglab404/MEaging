setwd("~/data/Entropy/")
dir.create("Result")
#library(bsseq)
library(stringr)
library(magrittr)
library(methylKit)
library(data.table)
library(ggplot2)
library(ggsignif)
#library(BSgenome.Hsapiens.UCSC.hg38)
sampleinfo<-read.table("Sampleinfo133")
sampleinfo$ID<-paste0("sample",sampleinfo$ID)
dir.create("./Result/Summary")
# Ent_files<-paste0("Ent89/Entropy",sampleinfo$ID) %>% sort
# Entropy1<-fread(Ent_files[1])
# Entropy1<-Entropy1[!startsWith(Entropy1$segment,"_"),2:8]
# Entropy1$loc = paste0(Entropy1$Chr,":",Entropy1$segment)
# meth_aggre11<-data.table(loc=Entropy1$loc,Entropy1=(Entropy1$Level))
# 
# for(i in 2:133){
#   tmp<-fread(Ent_files[i])
#   tmp<-tmp[!startsWith(tmp$segment,"_"),2:8]
#   tmp$loc = paste0(tmp$Chr,":",tmp$segment)
#   tmp<-data.table(loc=tmp$loc,tmp=(tmp$Level))
#   colnames(tmp)[2]<-str_split(Ent_files[i],pattern = "/",simplify = T)[,2]
#   
#   meth_aggre11<-merge.data.table(meth_aggre11,tmp,by="loc",all.x=T,all.y = T)
#   cat(i);cat(" ")
#   gc()
# }
# saveRDS(meth_aggre11,file="Result/Meth_133samples.aggregate.RDS")

######
# sampleinfo$ID<-paste0("Entropy",sampleinfo$ID)
# IDs<-sampleinfo$ID
# CEN<-sampleinfo[sampleinfo$Age>=95,]$ID
# YC<-sampleinfo[sampleinfo$Age <70,]$ID
# EC<-setdiff(IDs,c(CEN,YC))
# 
# methCEN<-meth_aggre11[,..CEN]
# methYC<-meth_aggre11[,..YC]
# methEC<-meth_aggre11[,..EC]
# 
# na_CEN<-apply(methCEN,1,FUN = function(x){return(sum(!is.na(x)))})
# na_YC<-apply(methYC,1,FUN = function(x){return(sum(!is.na(x)))})
# na_EC<-apply(methEC,1,FUN = function(x){return(sum(!is.na(x)))})
# 
# save(na_CEN,na_EC,na_YC,file="NAnums_133.RData")
# 
# #####Estimate sample number #####
# 
# # sum((na_CEN+na_EC+na_YC)>0.6*133)
# # na_EC[(na_CEN+na_EC+na_YC)>0.5*133] %>% summary
# # na_YC[(na_CEN+na_EC+na_YC)>0.5*133] %>% summary
# 
# sum((na_CEN > 0.5*length(CEN)) & (na_EC+na_YC)>55*0.5)
# na_EC[(na_CEN > 0.5*length(CEN)) & (na_EC+na_YC)>55*0.5] %>% summary #5102093 EC median 7
# 
# methRange<-meth_aggre11[(na_CEN > 0.5*length(CEN)) & (na_EC+na_YC)>55*0.5,]$loc
# ID<-c(CEN,YC,EC)
# meth_half<-meth_aggre11[(na_CEN > 0.5*length(CEN)) & (na_EC+na_YC)>55*0.5,..ID]
# save(meth_half,file="meth133_half.RData")
# rm(methCEN,methYC,methEC,tmp,Entropy1,Ent_files,i,ID,IDs,na_CEN,na_EC,na_YC)
# gc()
# save(methRange,file="Result/meth133_range.RData")
##### Chaos calculation 1: CV #####
#library(limma)
library(ggplot2)
meth133<-readRDS("Result/meth133.RDS")
colnames(meth133)<-sampleinfo$ID
load("./Result/methRange.RData")
IDs<-sampleinfo$ID %>% as.character()
CEN<-sampleinfo[sampleinfo$Age>=95,]$ID %>% as.character()
YC<-sampleinfo[sampleinfo$Age <=70,]$ID %>% as.character()
EC<-setdiff(IDs,c(CEN,YC))

CVs<-apply(meth133,2,FUN = sd,na.rm=T)/apply(meth133,2,FUN = mean,na.rm=T)
saveRDS(CVs,file="./Result/Summary/CVs133.RDS")
boxplot(CVs[YC],CVs[EC],CVs[CEN],col=c("green","yellow","darksalmon"))
ggplot(data.frame(age=sampleinfo$Age,CV=CVs),aes(x=age,y=CV))+geom_point(aes(color=age<=80))+
  scale_color_manual(values = c("TRUE"="blue1", 
                                "FALSE" ="red"))+geom_smooth(aes(color=age<=80),method = "lm")+
  theme_classic()


sampleinfo[sampleinfo$ID %in% YC,]$group <- "Younger Control"
sampleinfo[sampleinfo$ID %in% EC,]$group <- "Elder Control"
sampleinfo[sampleinfo$ID %in% CEN,]$group<-"Long-lived"
sampleinfo$group<-as.factor(sampleinfo$group)
levels(sampleinfo$group)<-c("Younger Control","Elder Control","Long-lived")

library(ggsignif)
ggplot(data.frame(group=sampleinfo$group,CV=CVs),aes(x=group,fill=group,y=CV))+
  geom_boxplot()+geom_point(position = position_jitter())+
  scale_fill_manual(values = c("Younger Control"="skyblue1", "Elder Control" = "gold4", "Long-lived"="red"))+
  geom_signif(comparisons = list(c(1,2),c(2,3)),color="black",test = ks.test)+
  theme_classic()
cor.test(sampleinfo[sampleinfo$Age<=80,]$Age,CVs[sampleinfo$Age<=80]) #p=0.000485 cor=0.458

##### Chaos calculation 2: Split 5/10 #####
EntMat<-function(x){
  testvec<-as.numeric(table(x))
  testvec<-testvec/sum(testvec)
  return(-(sum(testvec*log2(testvec))))
}

meth2<-meth133 %/% 10
apply(meth2,2,EntMat)->Ent10
sampleinfo$Ent10<-Ent10
cor.test(sampleinfo[sampleinfo$Age<95,]$Age,sampleinfo[sampleinfo$Age<95,]$Ent10,method = "pearson")
ggplot(data.frame(group=sampleinfo$group,Ent=sampleinfo$Ent10),aes(x=group,fill=group,y=Ent))+
  geom_boxplot()+geom_point(position = position_jitter())+
  scale_fill_manual(values = c("Younger Control"="skyblue1", "Elder Control" = "gold4", "Long-lived"="red"))+
  geom_signif(comparisons = list(c(1,2),c(2,3)),color="black",test = ks.test)+
  theme_classic()
ggplot(data.frame(age=sampleinfo$Age,Ent=sampleinfo$Ent10),aes(x=age,y=Ent))+geom_point(aes(color=age<=80))+
  scale_color_manual(values = c("TRUE"="blue1", 
                                "FALSE" ="red"))+geom_smooth(aes(color=age<=80),method = "lm")+
  theme_classic()


meth2<-meth133 %/% 5
apply(meth2,2,EntMat)->Ent20
sampleinfo$Ent20<-Ent20

ggplot(sampleinfo,aes(x=Age,y=Ent20))+geom_point()
cor.test(sampleinfo[sampleinfo$Age<95,]$Age,sampleinfo[sampleinfo$Age<95,]$Ent20,method = "pearson")
ggplot(sampleinfo,aes(x=Age,y=Ent20))+geom_point()
cor.test(sampleinfo[sampleinfo$Age<95,]$Age,sampleinfo[sampleinfo$Age<95,]$Ent20,method = "pearson")
ggplot(data.frame(group=sampleinfo$group,Ent=sampleinfo$Ent20),aes(x=group,fill=group,y=Ent))+
  geom_boxplot()+geom_point(position = position_jitter())+
  scale_fill_manual(values = c("Younger Control"="skyblue1", "Elder Control" = "gold4", "Long-lived"="red"))+
  geom_signif(comparisons = list(c(1,2),c(2,3)),color="black",test = t.test)+
  theme_classic()
ggplot(data.frame(age=sampleinfo$Age,Ent=Ent20),aes(x=age,y=Ent))+geom_point(aes(color=age<=80))+
  scale_color_manual(values = c("TRUE"="blue1", 
                                "FALSE" ="red"))+geom_smooth(aes(color=age<=80),method = "lm",formula = y~x)+
  theme_classic()

meth2<-meth133 %/% 20
apply(meth2,2,EntMat)->Ent5
sampleinfo$Ent5<-Ent5

ggplot(sampleinfo,aes(x=Age,y=Ent5))+geom_point()
cor.test(sampleinfo[sampleinfo$Age<95,]$Age,sampleinfo[sampleinfo$Age<95,]$Ent5,method = "pearson")
ggplot(sampleinfo,aes(x=Age,y=Ent5))+geom_point()
ggplot(data.frame(group=sampleinfo$group,Ent=sampleinfo$Ent5),aes(x=group,fill=group,y=Ent))+
  geom_boxplot()+geom_point(position = position_jitter())+
  scale_fill_manual(values = c("Younger Control"="skyblue1", "Elder Control" = "gold4", "Long-lived"="red"))+
  geom_signif(comparisons = list(c(1,2),c(2,3)),color="black",test = ks.test)+
  theme_classic()
ggplot(data.frame(age=sampleinfo$Age,Ent=sampleinfo$Ent5),aes(x=age,y=Ent))+geom_point(aes(color=age<=80))+
  scale_color_manual(values = c("TRUE"="blue1", 
                                "FALSE" ="red"))+geom_smooth(aes(color=age<=80),method = "lm")+
  theme_classic()

save(sampleinfo,file="Result/SplieEntropy.RData")

#sampleinfo[sampleinfo$Age %in% 70:80,]$group<-"Elder Control"
#sampleinfo[sampleinfo$Age %in% 40:69,]$group<-"Younger Control"

ggplot(sampleinfo,aes(x=group,y=Ent10,color=group))+geom_boxplot()

wilcox.test(sampleinfo[CEN,]$Ent10,sampleinfo[EC,]$Ent10)

##### Chaos calculation 3: Perumutation #####
library(statcomp)
PE_total_result<-data.frame(sample=sampleinfo$ID,
                            emb3=0,emb4=0,emb5=0,emb6=0,emb7=0)
for(i in 3:7){
  testOrd<-apply(X=meth133,2,FUN=ordinal_pattern_distribution,ndemb=i)
  testOrd<-testOrd[,sampleinfo$ID]
  PE_total<-apply(testOrd,2,permutation_entropy)
  PE_total_result[,i-1]<-PE_total
  cat(i)
}
ggplot(data.frame(group=sampleinfo$group,Ent=PE_total_result$emb7),aes(x=group,fill=group,y=Ent))+
  geom_boxplot()+geom_point(position = position_jitter())+
  scale_fill_manual(values = c("Younger Control"="skyblue1", "Elder Control" = "gold4", "Long-lived"="red"))+
  geom_signif(comparisons = list(c(1,2),c(2,3)),color="black")+
  theme_classic()


PE_total_result$group<-sampleinfo$group;
PE_total_result$Age<-sampleinfo$Age
PE_total_result$ID<-sampleinfo$ID
save(PE_total_result,file = "PE_result.RData")
library(reshape2)
PE_m<-melt(PE_total_result,measure.vars = 2:6,id.vars = 7:9)
PE_m[PE_m$group=="F1SP",]$group<-"Younger"
ggplot(PE_m,aes(x=variable,y=value,fill=group))+geom_boxplot()+theme_classic()+theme(text=element_text(size=18))


#SLides
PE_slide_result<-matrix(data=0,ncol=133,
                        nrow=nrow(meth_rb) %/% 1000 -4)

colnames(PE_slide_result)<-c(CEN,SP)

perc.meth<-perc.meth[,colnames(PE_slide_result)]

testFun<-function(j,data){
  Start<-(j-1)*1000+1
  End<-Start+4999
  testOrd<-apply(data[Start:End,],2,FUN=ordinal_pattern_distribution,ndemb=4)
  return(apply(testOrd,2,permutation_entropy))
}

x<-pblapply(X=1:nrow(PE_slide_result),FUN=testFun,data=meth_rb,cl=cl)
xt<-do.call("rbind",x)