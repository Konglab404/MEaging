setwd("~/data/Entropy/Result/")
library(stringr)
library(magrittr)
library(data.table)
library(SummarizedExperiment)
library(limma)

sampleinfo<-read.table("../Sampleinfo133")
load("../../DNAmAge/EpiDISH_176.RData")
Entropy_SE<-readRDS("Ent_SE_133_0.5_for3cats_rb.RDS")
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
###### Older vs Younger #####
# design_mat<-model.matrix(~0+(type=="Elder Control")+B+NK+CD4T+CD8T+Mono+Neutro+loc,
#                          data=blood[blood$type!="Long-lived",])
# # #design_mat<-model.matrix(~(type=="Elder Control"),data=cf[cf$type!="Long-lived",])
# colnames(design_mat)[1:2]<-c("Younger","Elder")
# contrast.matrix<-makeContrasts(contrasts = "Elder-Younger",levels = colnames(design_mat))
# fit_EvY<-lmFit(assay(Entropy_SE[,Entropy_SE$group!="Long-lived"]),design=design_mat)
# fit2_EvY<-contrasts.fit(fit_EvY,contrast.matrix)
# fit3_EvY<-eBayes(fit2_EvY)
# result_EvY<-topTable(fit3_EvY,coef=1,adjust="fdr",number = Inf)
# result_EvY$Delta<-apply(assay(Entropy_SE[,Entropy_SE$group=="Elder Control"]),1,mean,na.rm=T) - 
#                   apply(assay(Entropy_SE[,Entropy_SE$group=="Younger Control"]),1,mean,na.rm=T) 
# sum(result_EvY$Delta > 0.05 & result_EvY$adj.P.Val < 0.1, na.rm=T)
# 
# sum(result_OvY$logFC>0 & result_OvY$adj.P.Val<0.05,na.rm = T)
# OvY_hyper<-cbind(as.data.frame(rowRanges(Entropy_SE)),result_OvY)[result_OvY$adj.P.Val<0.05 & result_OvY$logFC>0,]
library(ppcor)
Ent<-assay(Entropy_SE[,Entropy_SE$Age<90])
# cf[cf$Age<90,]->EY
# EY$C<-as.numeric(EY$Age>70)
EY<-cbind(Entropy_SE[,Entropy_SE$Age<90]$Age,blood[colnames(Entropy_SE[,Entropy_SE$Age<90]),])
colnames(EY)[1]<-"Age"
EY$C<-as.numeric(EY$Age>70)

pcor_age<-function(x){ # x is num for lapply
  res<-data.frame(estimate=0,p.value=-1,statistic=0,n=0,gp=0,Method="N",Num=0)
  try({forTemp<-na.omit(data.frame(cbind((Ent[x,]),EY)))
       colnames(forTemp)[1]<-"Ent"
       res<-pcor.test(forTemp$Ent,forTemp$C,
                 forTemp[,c("B","NK","CD4T","CD8T","Mono","Neutro")], #forTemp[,c("Lymph","Mid","Gran")]
                 method = "kendall")
       res$Num<-x},silent = T)
  return(res)
}
library(pbapply)
library(parallel)
cl = makeForkCluster(4)
pblapply(X = 1:nrow(Ent),cl = cl,FUN = pcor_age)->Aging_DER
stopCluster(cl)
Aging_DER<-do.call("rbind",Aging_DER)
sum(Aging_DER$estimate>0 & Aging_DER$p.value<10e-3,na.rm = T) #
Aging_DER$p.adj<-p.adjust(Aging_DER$p.value, method = "fdr")
Aging_DER$qvalue<-qvalue(Aging_DER$p.value)$qvalues
Aging_DER$lfdr<-qvalue(Aging_DER$p.value)$lfdr
sum(Aging_DER$estimate>0 & Aging_DER$p.adj<0.05,na.rm = T) #121113
sum(Aging_DER$estimate<0 & Aging_DER$p.adj<0.05,na.rm = T) #51306
save(Aging_DER,file="tmp0630.RData")
# Aging_est<-apply(assay(Entropy_SE[Aging_DER$estimate>0 & Aging_DER$p.value<0.001,EC]),1,mean,na.rm=T) -
#            apply(assay(Entropy_SE[Aging_DER$estimate>0 & Aging_DER$p.value<0.001,YC]),1,mean,na.rm=T)
# sum(Aging_est>0.05) #3060


##### C vs Older #####
library(ppcor)
Ent<-assay(Entropy_SE[,Entropy_SE$Age>70])
# cf[cf$Age>70,]->EC
# EC$C<-as.numeric(EC$Age>=95)
ECs<-cbind(Entropy_SE[,Entropy_SE$Age>70]$Age,blood[colnames(Entropy_SE[,Entropy_SE$Age>70]),])
colnames(ECs)[1]<-"Age"
ECs$C<-as.numeric(ECs$Age>90)
pcor_Long<-function(x){ # x is num for lapply
  res<-data.frame(estimate=0,p.value=-1,statistic=0,n=0,gp=0,Method="N",Num=0)
  try({forTemp<-na.omit(data.frame(cbind((Ent[x,]),ECs)))
  colnames(forTemp)[1]<-"Ent"
  res<-pcor.test(forTemp$Ent,forTemp$C,
                 forTemp[,c("B","NK","CD4T","CD8T","Mono","Neutro")],
                 method = "kendall")  #Not Spearman
  res$Num<-x},silent = T)
  return(res)
}
library(pbapply)
library(parallel)
cl = makeForkCluster(6)
pblapply(X = 1:nrow(Ent),cl = cl,FUN = pcor_Long)->Long_DER
stopCluster(cl)
Long_DER<-do.call("rbind",Long_DER)
sum(Long_DER$estimate>0 & Long_DER$p.value<10e-3,na.rm = T) #
Long_DER$p.adj<-p.adjust(Long_DER$p.value, method = "fdr")
Long_DER$qvalue<-qvalue(Long_DER$p.value)$qvalues
Long_DER$lfdr<-qvalue(Long_DER$p.value)$lfdr
sum(Long_DER$estimate<0 & Long_DER$p.adj<0.05,na.rm = T) #62575
sum(Long_DER$estimate>0 & Long_DER$p.adj<0.05,na.rm = T) #41123

sum(Long_DER$estimate<0 & Long_DER$p.value<0.001 & Aging_DER$estimate>0 & Aging_DER$p.value<0.001,na.rm=T) #38217
sum(Long_DER$estimate<0 & Long_DER$p.adj<0.05 & Aging_DER$estimate>0 & Aging_DER$p.adj<0.05,na.rm=T) # 38923

##### Cen vs Elder #####
# design_mat<-model.matrix(~0+(type=="Long-lived")+B+NK+CD4T+CD8T+Mono+Neutro+loc,
#                          data=blood[blood$type!="Younger Control",])
# # #design_mat<-model.matrix(~(type=="Elder Control"),data=cf[cf$type!="Long-lived",])
# colnames(design_mat)[1:2]<-c("Elder","LLI")
# contrast.matrix<-makeContrasts(contrasts = "LLI-Elder",levels = colnames(design_mat))
# fit_LvE<-lmFit(assay(Entropy_SE[,Entropy_SE$group!="Younger Control"]),design=design_mat)
# fit2_LvE<-contrasts.fit(fit_LvE,contrast.matrix)
# fit3_LvE<-eBayes(fit2_LvE)
# result_LvE<-topTable(fit3_LvE,coef=1,adjust="fdr",number = Inf)
# result_LvE$Delta<-apply(assay(Entropy_SE[,Entropy_SE$group=="Long-lived"]),1,mean,na.rm=T) -
#                   apply(assay(Entropy_SE[,Entropy_SE$group=="Elder Control"]),1,mean,na.rm=T)
# sum(result_LvE$Delta < -0.1 & result_LvE$adj.P.Val < 0.05, na.rm=T)
# 
# sum(result_LvE$logFC<0 & result_LvE$adj.P.Val<0.05,na.rm = T)
#LvE_hyper<-cbind(as.data.frame(rowRanges(Entropy_SE)),result_LvE)[result_OvY$adj.P.Val<0.05 & result_OvY$logFC>0,]


#0####Longevity-specific 0Low-DER
Index<-Long_DER[Long_DER$estimate<0 & Long_DER$p.adj<0.05 & Aging_DER$estimate>0 & Aging_DER$p.adj<0.05,]$Num %>%
  na.omit
DER_CEN<-rowRanges(Entropy_SE)[Index,]
DER_CEN_ent<-(Entropy_SE)[Index,]

save(DER_CEN,file="DER_CEN_20230701.RData")

#0####Longevity-specific high-DER
Index<-Long_DER[Long_DER$estimate>0 & Long_DER$p.adj<0.05 & Aging_DER$estimate<0 & Aging_DER$p.adj<0.05,]$Num %>%
  na.omit
DER_CEN_high<-rowRanges(Entropy_SE)[Index,]
save(DER_CEN_high,file="DER_CEN_high_20230701.RData")

# #### Regression #####
# library(parallel)
# library(pbapply)
fit_simple<-lm(as.numeric(assay(Entropy_SE[1,]))~0+cf[,9]) #Model only Age
fit_multiple<-lm(as.numeric(assay(Entropy_SE[1,]))~0+cf[,1]+cf[,2]+cf[,3]+cf[,4]+cf[,5]+cf[,6]+cf[,9]) #Model with blood+age
# 
# rownames(cf)<-colnames(Entropy_SE)
# cf_Control<-cbind(cf[rownames(EY),],blood[rownames(EY),])
# 
# # Output: p-value for age; p-value for if blood can affect meth; OK flag
# cor_each<-function(i){ #perform linear regression for each CpG site and ages
#   res<-data.frame(p=0,blood_p=0,age_est=0,OK="N")
#   try({fit_simple<-lm(as.numeric(assay(Entropy_SE[i,rownames(EY)]))~cf_Control[,9])
#        fit_multiple<-lm(as.numeric(assay(Entropy_SE[i,rownames(EY)]))~cf_Control[,1]+cf_Control[,2]+
#                           cf_Control[,3]+cf_Control[,4]+cf_Control[,5]+cf_Control[,6]+cf_Control[,9])
#        compare<-anova(fit_multiple,fit_simple,test="LRT") #Anova to compare the two nested-models
#        compare$`Pr(>Chi)`[2]->res$blood_p                 #Lower p value suggest the blood is more associated with meth
#        summary(fit_simple)$coefficients[2,4]->res$p
#        summary(fit_simple)$coefficients[2,1]->res$age_est;res$OK<-'OK'
#        },silent = T)
#   return(res)
# }
# cor_each(1)
# 
# cat("Calculating");cat("\n")
# library(pbapply)
# library(parallel)
# cl = makeForkCluster(6)
# pblapply(X = 1:nrow(Entropy_SE),cl = cl,FUN = cor_each)->Aging_DER
# stopCluster(cl)
# DER_Aging<-do.call("rbind",Aging_DER)
# DER_Aging$loc<-Entropy_SE@rowRanges$loc
# 
# DER_Aging<-na.omit(DER_Aging)
# DER_Aging$q_age<-qvalue(DER_Aging$p)$qvalues
# DER_Aging$q_blood<-qvalue(DER_Aging$blood_p)$qvalues
# 
# sum(DER_Aging$q_age<0.05 & DER_Aging$q_blood>0.05 &DER_Aging$age_est>0)  
# #10296 aging-high-entropy 904 lowering
# 
# Aging_high_DER<-Entropy_SE[Entropy_SE@rowRanges$loc %in%
#                              DER_Aging[DER_Aging$q_age<0.05 & DER_Aging$q_blood>0.05 &DER_Aging$age_est>0,]$loc,]
# 
# cf_Control[,9]->Ages
# Aging_high_DER[,CEN]$Age->CENage
# # ttest_each<-function(i){ #perform linear regression for each CpG site and ages
# #   res<-data.frame(CEN=0,Expected=0,p=-1)
# #   try({fit_simple<-lm(as.numeric(assay(Aging_high_DER[i,c(EC,YC)]))~Ages)
# #        predict_entropy<-predict(fit_simple,data.frame(Ages=CENage))
# #        fit_ttest<-t.test(assay(Aging_high_DER[i,CEN]),predict_entropy)
# #        res<-data.frame(CEN=fit_ttest$estimate[1],Expected=fit_ttest$estimate[2],p=fit_ttest$p.value,row.names =
# #                          Aging_high_DER[i,]@rowRanges$loc)
# #   },silent = T)
# #   return(res)
# # }
# 
# ttest_each<-function(i){ #perform linear regression for each CpG site and ages
#   res<-data.frame(CvsExp=0,p=-1)
#   try({fit_simple<-lm(as.numeric(assay(Aging_high_DER[i,c(EC,YC)]))~Ages)
#   predict_entropy<-predict(fit_simple,data.frame(Ages=CENage))
#   tmp<-na.omit(data.frame(c=as.numeric(assay(Aging_high_DER[i,CEN])),es=predict_entropy))
#   fit_ttest<-t.test(tmp$c,tmp$es,paired = T)
#   res<-data.frame(CvsExp=fit_ttest$estimate[1],p=fit_ttest$p.value,row.names =
#                     Aging_high_DER[i,]@rowRanges$loc)
#   },silent = T)
#   return(res)
# }
# 
# ttest_each(8)
# 
# cl = makeForkCluster(8)
# pblapply(X = 1:nrow(Aging_high_DER),cl = cl,FUN = ttest_each)->CEN_der
# stopCluster(cl)
# DER_CEN<-do.call("rbind",CEN_der)
# 
# library(qvalue)
# DER_CEN$q<-p.adjust(DER_CEN$p,method = "BY")
# sum((DER_CEN$CEN - DER_CEN$Expected) < 0)
# 
# DER_CEN<-DER_CEN[DER_CEN$CvsExp < -0.2,]
# DER_CEN_gr<-rowRanges(Entropy_SE)[rowRanges(Entropy_SE)$loc %in% rownames(DER_CEN),]
# 
# save(DER_CEN,file="DER_CEN_20230301.RData")
################
library(LOLA)

chromhmm<-loadRegionDB("/data/user_data/wanght/HNLong3/LOLA_base/chromhmm") 
locDER_ChromHMM<-runLOLA(DER_CEN,rowRanges(Entropy_SE),chromhmm,cores=1)
locDER_ChromHMM_high<-runLOLA(DER_CEN_high,rowRanges(Entropy_SE),chromhmm,cores=1)
histone<-loadRegionDB("/data/user_data/wanght/HNLong3/LOLA_base/Histone_Broad_hg19")
locDER_histone<-runLOLA(DER_CEN,rowRanges(Entropy_SE),histone,cores=1)

save(locDER_ChromHMM,locDER_histone,file="LOLA_DER_cen_low.RData")

(apply(assay(DER_CEN_ent[,DER_CEN_ent$group == "Long-lived"]),1,median,na.rm=T) - 
    apply(assay(DER_CEN_ent[,DER_CEN_ent$group == "Elder Control"]),1,median,na.rm=T)) %>% summary

##############
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(annotatr)
load("DER_CEN_20230701.RData")

annos<-c('hg19_basicgenes')
annotations = build_annotations(genome = 'hg19', annotations = annos)
CEN_deranno<-annotate_regions(regions = DER_CEN,
                              annotations = annotations,ignore.strand = TRUE,quiet = FALSE)
CEN_deranno$annot %>% as.data.frame(stringsAsFactors=F)->CEN_deranno2
CEN_deranno2$loc<-CEN_deranno$loc
CEN_deranno2<-unique(CEN_deranno2[,c("loc","type","symbol")])
ggplot(rbind(data.frame(table(CEN_deranno2$type),type = "low-DERs"),
             data.frame(Var1="Distal",
                        Freq=38923-35929,
                        type="low-DERs"
                        )
             ),
       aes(fill=Var1,x=type,y=Freq))+geom_col(position = "fill")+
  coord_flip()+theme_classic()

backanno<-annotate_regions(regions = rowRanges(Entropy_SE),
                           annotations = annotations,ignore.strand = TRUE,quiet = FALSE)
backanno$annot %>% as.data.frame(stringsAsFactors=F)->backanno2
backanno2$loc<-backanno$loc
backanno2<-unique(backanno2[,c("loc","type","symbol")])
pros_back<-annotations[annotations$type == "hg19_genes_promoters"] %>% reduce

#######################################################3
annos<-c('hg19_basicgenes')
annotations = build_annotations(genome = 'hg19', annotations = annos)
CEN_deranno<-annotate_regions(regions = DER_CEN,
                              annotations = annotations,ignore.strand = TRUE,quiet = FALSE)
CEN_deranno<-cbind(CEN_deranno$loc,as.data.frame(CEN_deranno$annot))
CEN_deranno<-unique(CEN_deranno[,c("seqnames","start","type","symbol")])

CEN_deranno$symbol %>% unique() %>% length

CEN_deranno[grep(CEN_deranno$type,pattern = "promo"),]$symbol %>% 
  table %>% as.data.frame(stringsAsFactors=F)->promoter_table
promoter_table[promoter_table$Freq>=3,]$. %>% unique %>% cat
save(promoter_table,file="annotr_promoter_genes_low_DERs.RData")

Back_anno<-annotatePeak(peak = rowRanges(Entropy_SE),tssRegion = c(-3000,3000),
                          TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,annoDb = "org.Hs.eg.db")

Back_anno<-Back_anno@anno
Back_anno[grep(Back_anno$annotation,pattern = "Promo"),]$SYMBOL %>% unique %>% cat

#####Plot low DER
# library(ggplot2)
# ggplot(data.frame(age=Entropy_SE$Age,
#                   Entropy=apply(assay(Entropy_SE[Aging_DER[Aging_DER$estimate>0 & Aging_DER$p.adj<0.1,]$Num,]),
#                                 2,FUN=mean,na.rm=T)),
#        aes(x=age,y=Entropy))+
#   geom_point()+geom_smooth()+theme_bw()
