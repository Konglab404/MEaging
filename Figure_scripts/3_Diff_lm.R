setwd("~/data/Entropy/Result/")
library(stringr)
library(magrittr)
library(data.table)
library(SummarizedExperiment)
library(limma)

sampleinfo<-read.table("../Sampleinfo133")
load("../../DNAmAge/EpiDISH_176.RData")
Entropy_SE<-readRDS("Ent_SE_133_0.5_for3cats.RDS")
#assay(Entropy_SE)[assay(Entropy_SE)<0]<-0
cf<-EpiDISH_est$estF[str_replace(colnames(Entropy_SE),pattern = "Entropy",replacement = "sample_"),]
cf<-as.data.frame(cf)
cf$type<-Entropy_SE$group
cf$Age<-Entropy_SE$Age

blood<-read.table("../../DNAmAge/Blood_covariates.txt",header = T)
rownames(blood)<-paste0("Entropy",blood$ID)
blood<-blood[colnames(Entropy_SE),]
#load("meth133_range.RData")

blood<-cbind(blood,cf)
blood$loc<-colData(Entropy_SE)$loc
blood<-cbind(blood,model.matrix(~0+sampleinfo$batch))
colnames(blood)[15:17]<-c("batchA","batchB","batchC")
##### Older vs Younger #####
EY<-cbind(Entropy_SE[,Entropy_SE$Age<90]$Age,blood[colnames(Entropy_SE[,Entropy_SE$Age<90]),])
rownames(cf)<-colnames(Entropy_SE)
cf_Control<-cbind(cf[rownames(EY),],blood[rownames(EY),])
# Output: p-value for age; p-value for if blood can affect meth; OK flag
cor_each<-function(i){ #perform linear regression for each CpG site and ages
  res<-data.frame(p=0,blood_p=0,age_est=0,OK="N")
  try({fit_simple<-lm(as.numeric(assay(Entropy_SE[i,rownames(EY)]))~cf_Control[,9])
       fit_multiple<-lm(as.numeric(assay(Entropy_SE[i,rownames(EY)]))~cf_Control[,1]+cf_Control[,2]+
                          cf_Control[,3]+cf_Control[,4]+cf_Control[,5]+cf_Control[,6]+cf_Control[,9])
       compare<-anova(fit_multiple,fit_simple,test="F") #Anova to compare the two nested-models
       compare$`Pr(>F)`[2]->res$blood_p                 #Lower p value suggest the blood is more associated with meth
       summary(fit_simple)$coefficients[1,4]->res$p
       summary(fit_simple)$coefficients[1,1]->res$age_est;res$OK<-'OK'
       },silent = T)
  return(res)
}
cor_each(1)

cat("Calculating");cat("\n")
library(pbapply)
library(parallel)
cl = makeForkCluster(6)
pblapply(X = 1:nrow(Entropy_SE),cl = cl,FUN = cor_each)->Aging_DER
stopCluster(cl)
DER_Aging<-do.call("rbind",Aging_DER)
DER_Aging$loc<-Entropy_SE@rowRanges$loc

DER_Aging<-na.omit(DER_Aging)
DER_Aging$q_age<-p.adjust(DER_Aging$p,method = "bonferroni")
DER_Aging$q_blood<-p.adjust(DER_Aging$blood_p,method="bonferroni")
sum(DER_Aging$q_age<0.05 & DER_Aging$blood_p>0.05 & DER_Aging$age_est<0)
sum(DER_Aging$q_age < 5e-8 & DER_Aging$blood_p > 0.05 & DER_Aging$age_est>(0.2/35))
Aging_high_DER<-Entropy_SE[DER_Aging$q_age<0.05 & DER_Aging$blood_p>0.05 &DER_Aging$age_est>0,]
rowData(Aging_high_DER)$q_age<-DER_Aging[DER_Aging$q_age<0.05 & DER_Aging$blood_p>0.05 &DER_Aging$age_est>0,]$q_age
rowData(Aging_high_DER)$coef<-DER_Aging[DER_Aging$q_age<0.05 & DER_Aging$blood_p>0.05 &DER_Aging$age_est>0,]$age_est

#######Setup for predict##############
EC<-colData(Entropy_SE)[Entropy_SE$group == "Elder Control",]$ID
YC<-colData(Entropy_SE)[Entropy_SE$group == "Younger Control",]$ID
CEN<-colData(Entropy_SE)[Entropy_SE$group == "Long-lived",]$ID
Ages<-colData(Entropy_SE)[c(EC,YC),]$Age
CENage<-colData(Entropy_SE)[CEN,]$Age

######################################


ttest_each<-function(i){ #perform linear regression for each CpG site and ages
  res<-data.frame(CvsExp=0,p=-1)
  try({fit_simple<-lm(as.numeric(assay(Aging_high_DER[i,c(EC,YC)]))~Ages)
  predict_entropy<-predict(fit_simple,data.frame(Ages=CENage))
  tmp<-na.omit(data.frame(c=as.numeric(assay(Aging_high_DER[i,CEN])),es=predict_entropy))
  fit_ttest<-t.test(tmp$c,tmp$es,paired = T)
  res<-data.frame(CvsExp=fit_ttest$estimate[1],p=fit_ttest$p.value,
                  row.names =Aging_high_DER[i,]@rowRanges$loc)
  },silent = T)
  return(res)
}

ttest_each(9)

cl = makeForkCluster(8)
pblapply(X = 1:nrow(Aging_high_DER),cl = cl,FUN = ttest_each)->CEN_der
stopCluster(cl)
DER_CEN<-do.call("rbind",CEN_der)

library(qvalue)
DER_CEN$q<-p.adjust(DER_CEN$p,method = "bonferroni")
sum(DER_CEN$CvsExp<(-0) & DER_CEN$q<5e-8 & rowData(Aging_high_DER)$q_age<5e-8) #16118

CEN_low_DER<-Entropy_SE[rowRanges(Entropy_SE)$loc %in% 
                          rownames(DER_CEN)[DER_CEN$CvsExp<(-0) & DER_CEN$q<5e-8 & rowData(Aging_high_DER)$q_age<0.05],]
rowData(CEN_low_DER)<-DER_CEN[DER_CEN$CvsExp<(-0) & DER_CEN$q<5e-8 & rowData(Aging_high_DER)$q_age<5e-8,]

save(CEN_low_DER,DER_CEN,file="DER_CEN_20230723_lmfit.RData")
#save(CEN_low_DER,file="DER_CEN_20230723_lmfit.RData")
hist(CEN_low_DER@rowRanges$CvsExp,xlab = "Residuals in centenarians", main=NULL)

#################
load("DER_CEN_20230723_lmfit.RData")

library(qqman)
qq(DER_CEN$p)
p_value=DER_CEN$p
z = qnorm(p_value/ 2)
lambda = round(median(z^2, na.rm = TRUE) / 0.454, 3)

DER_CEN$chr<-str_split(rownames(DER_CEN),pattern = ":|r",simplify = T)[,2]
DER_CEN$chr<-as.numeric(DER_CEN$chr)
DER_CEN$pos<-str_split(rownames(DER_CEN),pattern = ":|_",simplify = T)[,3]
DER_CEN$pos<-as.numeric(DER_CEN$pos)
DER_CEN$loc<-rownames(DER_CEN)
manhattan(DER_CEN, chr = "chr", bp = "pos", p = "q", snp = "loc",
          col = c("grey", "deepskyblue"), chrlabs = NULL,cex = 0.4,
          suggestiveline = -log10(5e-08), genomewideline = -log10(5e-08),
          highlight = NULL, logp = TRUE, annotatePval = NULL,
          annotateTop = TRUE)

#####################
library(LOLA)

chromhmm<-loadRegionDB("/data/user_data/wanght/HNLong3/LOLA_base/chromhmm") 
locDER_ChromHMM<-runLOLA(rowRanges(CEN_low_DER),rowRanges(Entropy_SE),chromhmm,cores=1)
#locDER_ChromHMM_high<-runLOLA(DER_CEN_high,rowRanges(Entropy_SE),chromhmm,cores=1)
histone<-loadRegionDB("/data/user_data/wanght/HNLong3/LOLA_base/Histone_Broad_hg19")
locDER_histone<-runLOLA(rowRanges(CEN_low_DER),rowRanges(Entropy_SE),histone,cores=1)

save(locDER_ChromHMM,locDER_histone,file="LOLA_DER_cen_low.RData")

(apply(assay(DER_CEN_ent[,DER_CEN_ent$group == "Long-lived"]),1,median,na.rm=T) - 
    apply(assay(DER_CEN_ent[,DER_CEN_ent$group == "Elder Control"]),1,median,na.rm=T)) %>% summary
