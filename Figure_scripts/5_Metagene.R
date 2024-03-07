setwd("~/data/Entropy/Result/")
Entropy_SE<-readRDS("Ent_SE_133_0.5_for3cats_rb.RDS")

cvs<-apply(assay(Entropy_SE),1,function(x){return(sd(x,na.rm = T)/mean(x,na.rm=T))})
quantile(cvs,seq(0,1,0.05))->qt


#PCA for outlier filter
pca.meth<-Entropy_SE[,] %>% assay %>% na.omit %>% t %>% prcomp
#pdf("111_samplePCA.pdf")
library(ggbiplot)
ggbiplot(pca.meth,obs.scale = 1,var.scale = 1,group = sampleinfo$group,
         ellipse = T,var.axes = F,ellipse.prob = 0.9) +
  scale_color_manual(values = c("Younger Control"="#61ccffff",
                                "Elder Control"="#b78214ff",
                                "Long-lived"="#f8130eff"))+
  theme(legend.direction = 'horizontal',legend.position = 'top') +
  geom_vline(xintercept = c(0), linetype = 6,color = "black") +
  geom_hline(yintercept = c(0),linetype = 6,color = "black") +
  theme_bw() + theme(panel.grid = element_line(colour = NA))
#save(pca.meth,file="PCA_133_ent.RData")
ggbiplot(pca.meth,obs.scale = 1,var.scale = 1,group = sampleinfo$batch,
         ellipse = T,var.axes = F,ellipse.prob = 0.99) + theme_classic()


dis_ent<-dist(Entropy_SE[,] %>% assay %>% na.omit %>% t,method = "manhattan") #recom：0.5-0.95
mds_x = cmdscale(dis_ent,k = 3)
mds_x = data.frame(mds_x)
xy<-cbind(mds_x,Entropy_SE$group)
colnames(xy)[4]<-"Group"

ggplot(xy, aes(x=X1,y=X3,colour = Group))+geom_point()+stat_ellipse(level = 0.5) +
  scale_color_manual(values = c("Younger Control"="#61ccffff",
                                "Elder Control"="#b78214ff",
                                "Long-lived"="#f8130eff"))+theme_bw()

ggplot(xy, aes(x=X1,y=X2,colour = Entropy_SE$batch))+geom_point()+stat_ellipse(level = 0.95)
save(xy,file = "Summary/MDS_back.RData")

library(Rtsne)
set.seed(19940731)
tsne_out <- Rtsne(
  (Entropy_SE[cvs<qt["95%"] & cvs>qt["50%"],] %>% assay %>% na.omit %>% t),
  dims = 2,
  pca = T,
  perplexity = 10,
  theta = 0.0,
  max_iter = 50
)
tsne_result = as.data.frame(tsne_out$Y)
colnames(tsne_result) = c("tSNE1","tSNE2")
ggplot(tsne_result,aes(tSNE1,tSNE2,color=Entropy_SE$group)) +
  geom_point()


###Distance to TSS from Fang Yuqi et al NAR 2023
dist_calc<-function(gr_in,gr_feature){
  
  dist_nearest=nearest(resize(gr_in,1,fix="center"),resize(gr_feature,1,fix="start"))
  
  gr_non_na=which(!is.na(dist_nearest))
  gr_feature=gr_feature[dist_nearest[gr_non_na]]
  gr_in=gr_in[gr_non_na]
  gr_in$dist=NA
  gr_in$gene=NA
  
  sgn <- as.integer(ifelse(strand(gr_feature)=="+",1,-1))
  gr_in$dist=sgn*(start(resize(gr_in,1,fix="center"))-start(resize(gr_feature,1,fix="start")))
  gr_in$gene=gr_feature$gene_name
  
  return(gr_in)
}

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
TSS_hg19<-promoters(TxDb.Hsapiens.UCSC.hg19.knownGene,upstream = 1,downstream = 1)
TSS_hg19$gene_name<-TSS_hg19$tx_name

Dist<-dist_calc(rowRanges(Entropy_SE),TSS_hg19)
save(Dist,file="Distance_2_TSS.RData")

CEN_mean<-apply(assay(Entropy_SE[,Entropy_SE$group == "Long-lived"]),1,mean,na.rm=T)
YC_mean<-apply(assay(Entropy_SE[,Entropy_SE$group == "Younger Control"]),1,mean,na.rm=T)
EC_mean<-apply(assay(Entropy_SE[,Entropy_SE$group == "Elder Control"]),1,mean,na.rm=T)

Dist$CEN<-CEN_mean
Dist$YC<-YC_mean
Dist$EC<-EC_mean

rm(CEN_mean,YC_mean,EC_mean,Entropy_SE);gc()

sm_len=100

Dist_5000<-Dist[Dist$dist %in% c(-3000:3000),]
Dist_5000$dist_smooth<-sm_len*(-Dist_5000$dist %/% sm_len)
plot_Dist<-aggregate(data.frame(CEN=Dist_5000$CEN,
                                EC=Dist_5000$EC,
                                YC=Dist_5000$YC),
                     by=list(Dist_5000$dist_smooth),FUN=mean,na.rm=T) %>%
  melt(measure.vars=2:4,value.name = "Ent",id.vars=1,variable.name="Type",)
ggplot(plot_Dist,aes(x=Group.1,y=Ent,color=Type))+geom_smooth(method = "gam",level=0.8)+theme_bw()


plot(plot_Dist$Group.1,plot_Dist$CEN,"l",col="black")
lines(plot_Dist$Group.1,plot_Dist$EC,"l",col="red")
lines(plot_Dist$Group.1,plot_Dist$YC,"l",col="blue")
################################################################## DER
load("DER_CEN_20230701.RData")
Dist_5000_low<-Dist[(Dist$dist %in% c(-3000:3000)) & (Dist$loc %in% DER_CEN$loc),]
Dist_5000_low$dist_smooth<-sm_len*(-Dist_5000_low$dist %/% sm_len)
plot_Dist_low<-aggregate(data.frame(CEN=Dist_5000_low$CEN,
                                EC=Dist_5000_low$EC,
                                YC=Dist_5000_low$YC),
                     by=list(Dist_5000_low$dist_smooth),FUN=mean,na.rm=T) %>%
  melt(measure.vars=2:4,value.name = "Ent",id.vars=1,variable.name="Type",)
ggplot(plot_Dist_low,aes(x=Group.1,y=Ent,color=Type))+geom_smooth(method = "gam",level=0.8)+theme_bw()

###
DER_CEN_ent<-Entropy_SE[Entropy_SE@rowRanges$loc %in% DER_CEN$loc,]
ES<-readRDS("Ent_SE_133_0.5_for3cats.RDS")

xp<-data.frame(YC = apply(assay(DER_CEN_ent[,DER_CEN_ent$group=="Younger Control"]),1,min,na.rm=T),
               EC = apply(assay(DER_CEN_ent[,DER_CEN_ent$group=="Elder Control"]),1,max,na.rm=T),
               CEN = apply(assay(DER_CEN_ent[,DER_CEN_ent$group=="Long-lived"]),1,median,na.rm=T))

cvs<-apply(assay(DER_CEN_ent),1,function(x){return(sd(x,na.rm = T)/mean(x,na.rm=T))})
quantile(cvs,seq(0,1,0.05))->qt
sds<-apply(assay(DER_CEN_ent),1,function(x){return(sd(x,na.rm=T))})

cvss<-sort(cvs)

dis_ent<-dist(DER_CEN_ent[names(cvss[1:20000]),] %>% assay %>% na.omit %>% t,method = "manhattan") #recom：cvs0.5-3.1448
mds_x = cmdscale(dis_ent,k = 3)
mds_x = data.frame(mds_x)
xy<-cbind(mds_x,Entropy_SE$group)
colnames(xy)[4]<-"Group"

ggplot(xy, aes(x=X1,y=X2,z=X3,colour = Group))+geom_point()+stat_ellipse(level = 0.95) +
  scale_color_manual(values = c("Younger Control"="#61ccffff",
                                "Elder Control"="#b78214ff",
                                "Long-lived"="#f8130eff"))+theme_bw()+axes_3

my_color = c("#61ccffff", "#b78214ff", "#f8130eff")
colors = my_color[as.numeric(Entropy_SE$group)]
p1 = scatterplot3d(xy[,1:3],color = colors,main="iris",pch = 16)
legend(p1$xyz.convert(8.5, 2.5, 5), legend = levels(iris$Species),
       col = my_color, pch = 16)
plot3d(xy[,1:3],color = colors,main="iris",pch = 16)

load("../../DNAmAge/EpiDISH_176.RData")
#assay(Entropy_SE)[assay(Entropy_SE)<0]<-0
cf<-EpiDISH_est$estF[str_replace(colnames(Entropy_SE),pattern = "Entropy",replacement = "sample_"),]
cf<-as.data.frame(cf)
DER_CEN_ent_raw@colData<-cbind(DER_CEN_ent_raw@colData,cf)
Agings<-DER_CEN_ent[,DER_CEN_ent_raw$group != "Long-lived"]
Ages<-Agings$Age
Loc<-Agings$loc
B<-Agings$B;NK<-Agings$NK;CD4T<-Agings$CD4T;CD8T<-Agings$CD8T;Mono<-Agings$Mono;Neutro<-Agings$Neutro
CENs<-DER_CEN_ent[,DER_CEN_ent_raw$Age>90]


#cf_Control[,9]->Ages
#Aging_high_DER[,CEN]$Age->CENage
ttest_each<-function(i){ #perform linear regression for each CpG site and ages
  res<-data.frame(CvsExp=0,p=-1)
  try({
  fit_simple<-lm(as.numeric(assay(Agings)[i,])~Ages+B+NK+CD4T+CD8T+Mono+Neutro)
  predict_entropy<-predict(fit_simple,data.frame(Ages=CENs$Age,
                                                 Loc=CENs$loc,
                                                 B=CENs$B,
                                                 NK=CENs$NK,
                                                 CD4T=CENs$CD4T,
                                                 CD8T=CENs$CD8T,
                                                 Mono=CENs$Mono,
                                                 Neutro=CENs$Neutro))
  tmp<-na.omit(data.frame(c=as.numeric(assay(CENs[i,])),es=predict_entropy))
  fit_ttest<-t.test(tmp$c,tmp$es,paired = T)
  res<-data.frame(CvsExp=median(tmp$c-tmp$es),p=fit_ttest$p.value,row.names =
                    Agings[i,]@rowRanges$loc)
  },silent = T)
  return(res)
}

cl = makeForkCluster(8)
pblapply(X = 1:nrow(Agings),cl = cl,FUN = ttest_each)->DER_CEN_diff
stopCluster(cl)
DER_CEN_diff<-do.call("rbind",DER_CEN_diff)
sum(DER_CEN_diff$CvsExp<0)
