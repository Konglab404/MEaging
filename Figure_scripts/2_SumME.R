setwd("~/data/Entropy/Result/")
library(stringr)
library(magrittr)
library(data.table)
library(SummarizedExperiment)
sampleinfo<-read.table("../Sampleinfo133")
sampleinfo$ID<-paste0("Entropy",sampleinfo$ID)
IDs<-sampleinfo$ID %>% as.character()
CEN<-sampleinfo[sampleinfo$Age>=95,]$ID %>% as.character()
YC<-sampleinfo[sampleinfo$Age <=70,]$ID %>% as.character()
EC<-setdiff(IDs,c(CEN,YC))
sampleinfo[sampleinfo$ID %in% YC,]$group <- "Younger Control"
sampleinfo[sampleinfo$ID %in% EC,]$group <- "Elder Control"
sampleinfo[sampleinfo$ID %in% CEN,]$group<-"Long-lived"
sampleinfo$group<-factor(sampleinfo$group,levels = c("Younger Control","Elder Control","Long-lived"))
#levels(sampleinfo$group)<-c("Younger Control","Elder Control","Long-lived")
rownames(sampleinfo)<-sampleinfo$ID

Entropy<-readRDS("Entropy_133samples.aggregate.RDS")


IDs<-c(sampleinfo$ID)
methRange<-Entropy$loc

Entropy<-Entropy[,..IDs]

CEN<-sampleinfo[sampleinfo$Age>=95,]$ID
YC<-sampleinfo[sampleinfo$Age <=70,]$ID
EC<-setdiff(IDs,c(CEN,YC))

methRange2<-GRanges(seqnames = str_split(methRange,pattern = ":|_",simplify = T)[,1],
                    IRanges(start = as.numeric(str_split(methRange,pattern = ":|_",simplify = T)[,2]),
                            end = as.numeric(str_split(methRange,pattern = ":|_",simplify = T)[,4])),
                    strand="*")
methRange2$loc<-methRange

rownames(sampleinfo)<-sampleinfo$ID
Entropy_SE<-SummarizedExperiment(assays = Entropy,colData = sampleinfo,rowRanges = methRange2)
save(Entropy_SE,file="Entropy_133_SummarizedExperiment.RData")

na_CEN<-apply(Entropy[,..CEN],1,FUN = function(x){return(sum(!is.na(x)))})
na_YC<-apply(Entropy[,..YC],1,FUN = function(x){return(sum(!is.na(x)))})
na_EC<-apply(Entropy[,..EC],1,FUN = function(x){return(sum(!is.na(x)))})

saveRDS(Entropy_SE[na_CEN>0.5*length(CEN) & na_EC>0.5*length(EC) & na_YC>0.5*length(YC),],
        file="Ent_SE_133_0.5_for3cats.RDS")
saveRDS(Entropy_SE[(na_CEN+na_EC+na_YC)>0.6*133,],
        file="Ent_SE_133_0.6_for133.RDS")
saveRDS(Entropy_SE[(na_CEN > 0.5*length(CEN)) & (na_EC+na_YC)>27,],
        file="Ent_SE_133_0.5_forCENandControls.RDS")
save(na_CEN,na_EC,na_YC,file="./NA_nums133.RData")


#####Summary #####
load("../../DNAmAge/EpiDISH_176.RData")
load("./NA_nums133.RData")
load("Entropy_133_SummarizedExperiment.RData")
library(limma)
Entropy<-removeBatchEffect(assay(Entropy_SE[na_CEN>0.5*length(CEN) & 
                                              na_EC>0.5*length(EC) & 
                                              na_YC>0.5*length(YC),]), #3642946
                           batch = sampleinfo$batch)

Entropy_each<-apply(Entropy,2,mean,na.rm=T)
names(Entropy_each)<-sampleinfo$ID
library(ggplot2)
sampleinfo$ME<-Entropy_each
library(ggsignif)
ggplot(data.frame(group=sampleinfo$group,ME=Entropy_each),aes(x=group,fill=group,y=ME))+
  geom_boxplot()+geom_point(position = position_jitter())+
  scale_fill_manual(values = c("Younger Control"="skyblue1", "Elder Control" = "gold4", "Long-lived"="red"))+
  geom_signif(comparisons = list(c(1,2),c(2,3)),color="black",test = "ks.test")+
  theme_classic()
ggplot(data.frame(age=sampleinfo$Age,ME=Entropy_each),aes(x=age,y=ME))+geom_point(aes(color=age<=80))+
  scale_color_manual(values = c("TRUE"="blue1", 
                                "FALSE" ="red"))+
  geom_smooth(aes(color=age<=80),method = "lm")+
  theme_classic()
cor.test(sampleinfo[sampleinfo$Age<95,]$Age,sampleinfo[sampleinfo$Age<95,]$ME,method = "pearson") #0.275 0.044
save(Entropy_each,file="ME_each_133.RData")

Ent_rm<-SummarizedExperiment(assays = Entropy,
                             rowRanges = rowRanges(Entropy_SE)[na_CEN>0.5*length(CEN) & 
                                                                 na_EC>0.5*length(EC) & 
                                                                 na_YC>0.5*length(YC),],
                             colData = sampleinfo)

saveRDS(Ent_rm,file="Ent_SE_133_0.5_for3cats_rb.RDS") 
#saveRDS(Ent_rm,file="Ent_SE_133_0.6_for133_rb.RDS") ##################

rm(Entropy,Entropy_SE,na_CEN,na_EC,na_YC,opd,methRange,testOrd,meth133,methRange2);gc()

CEN_mean<-apply(assay(Ent_rm[,CEN]),1,mean,na.rm=T)
YC_mean<-apply(assay(Ent_rm[,YC]),1,mean,na.rm=T)
EC_mean<-apply(assay(Ent_rm[,EC]),1,mean,na.rm=T)

hist(CEN_mean,breaks = 50,xlim = c(0,1))
hist(YC_mean,breaks = 50,xlim = c(0,1))
hist(EC_mean,breaks = 50,xlim = c(0,1))

ggplot(data=data.frame(value=c(CEN_mean,YC_mean,EC_mean),
                       type=rep(c("CEN","Young","EC"),each=6087544)))+
  geom_density(aes(x=value,color=type))+ylim(c(0,4))+theme_classic()+
  scale_color_manual(values = c("CEN"="red","EC"="darkgoldenrod","Young"="skyblue"))

save(CEN_mean,YC_mean,EC_mean,file="Mean_ME.RData")


#####Circos########
Entropy_SE<-readRDS("Ent_SE_133_0.5_for3cats_rb.RDS")
Circos_mean<-rowRanges(Entropy_SE)
Circos_mean$CEN<-CEN_mean
Circos_mean$YC<-YC_mean
Circos_mean$EC<-EC_mean

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
seqlengths(Circos_mean)<-seqlengths(TxDb.Hsapiens.UCSC.hg19.knownGene)[names(seqlengths(Circos_mean))]
chrs   =as.character(unique(seqnames(Circos_mean)))
widths =sapply(chrs,function(x,y) max(end(y[seqnames(y)==x,])),Circos_mean  )
all.wins=GRanges()
win.size=5000000; step.size=5000000
for(i in 1:length(chrs))
{
  # get max length of feature covered chromosome
  max.length=max(IRanges::end(Circos_mean[seqnames(Circos_mean)==chrs[i],])) 
  
  #get sliding windows with covered CpGs
  numTiles=floor(  (max.length-(win.size-step.size) )/step.size )+1
  numTiles=ifelse(numTiles<1, 1,numTiles)
  temp.wins=GRanges(seqnames=rep(chrs[i],numTiles),
                    ranges=IRanges(start=1+0:(numTiles-1)*step.size,
                                   width=rep(win.size,numTiles)) )
  all.wins=suppressWarnings(c(all.wins,temp.wins))
}

x<-findOverlaps(all.wins,Circos_mean)
Circos_mean$bins[x@to]<-x@from
Circ_binned<-aggregate(as.data.frame(Circos_mean)[,c("CEN","YC","EC")],list(Circos_mean$bins),FUN=mean)
all.wins2<-cbind(all.wins[Circ_binned$Group.1,],Circ_binned)
all.wins2$end<-all.wins2$end - 400000

dir.create("Circos")

#CEN 5mb methylation
my5m.cen<-data.frame(chrom=str_replace(all.wins2$seqnames,pattern = "chr",replacement = "hs"),
                     start=all.wins2$start,
                     end = all.wins2$end,
                     meth = all.wins2$CEN)
write.table(my5m.cen,file="Circos/5m_CEN_meth",sep="\t",row.names = F,quote = F,col.names = F)

#YC 5mb methylation
my5m.yc<-data.frame(chrom=str_replace(all.wins2$seqnames,pattern = "chr",replacement = "hs"),
                     start=all.wins2$start,
                     end = all.wins2$end,
                     meth = all.wins2$YC)
write.table(my5m.yc,file="Circos/5m_YC_meth",sep="\t",row.names = F,quote = F,col.names = F)

#EC 5mb methylation
my5m.ec<-data.frame(chrom=str_replace(all.wins2$seqnames,pattern = "chr",replacement = "hs"),
                     start=all.wins2$start,
                     end = all.wins2$end,
                     meth = all.wins2$EC)
write.table(my5m.ec,file="Circos/5m_EC_meth",sep="\t",row.names = F,quote = F,col.names = F)

##### Metagene #####
library(annotatr)
annos<-c('hg19_basicgenes')
annotations = build_annotations(genome = 'hg19', annotations = annos)
anno2<-annotations[!is.na(annotations$gene_id),]

Ent_anno<-annotate_regions(regions = rowRanges(Entropy_SE),
                      annotations = annotations,
                      ignore.strand = TRUE,
                      quiet = FALSE)
Ent_anno$annot<-Ent_anno$annot$type
Ent_anno<-unique(data.frame(seqnames=seqnames(Ent_anno),
                            start=start(Ent_anno),
                            end=end(Ent_anno),
                            type=Ent_anno$annot,loc=Ent_anno$loc)) %>% as.data.frame %>% as("GRanges")

MetaEnt<-data.frame()
for(i in 1:6){
  types<-unique(Ent_anno$type)[i]
  locs<-intersect(Ent_anno[Ent_anno$type == types,]$loc, Entropy_SE@rowRanges$loc)
  tmp<-data.frame(type=types,
                  ent = apply(assay(Entropy_SE[Entropy_SE@rowRanges$loc %in% locs,]),2,mean,na.rm=T),
                  group=sampleinfo$group)
  MetaEnt<-rbind(MetaEnt,tmp)
}
save(MetaEnt,file="Metagene.RData")
ggplot(MetaEnt,aes(fill=group,x=group,y=ent))+geom_boxplot()+facet_wrap(~type,nrow = 1)+
  geom_signif(comparisons = list(c(1,2),c(2,3)),test = ks.test,map_signif_level = c("***"=0.001,"**"=0.01,"*"=0.05))+
  theme_classic()+
  scale_fill_manual(values = c("Long-lived"="red","Elder Control"="darkgoldenrod","Younger Control"="skyblue"))

##### CpG OE promoter #####
Ent_anno<-annotate_regions(regions = rowRanges(Entropy_SE),
                           annotations = annotations,
                           ignore.strand = TRUE,
                           quiet = FALSE)
Ent_anno$annot<-Ent_anno$annot$id
load("promoter_CpGOE_annotatr_hg19.RData")
ProEnt<-data.frame()
for(i in 1:3){
  types<-unique(promoter$type)[i]
  locs<-intersect(Ent_anno[Ent_anno$annot %in% promoter[promoter$type == types,]$id,]$loc, Entropy_SE@rowRanges$loc)
  tmp<-data.frame(type=types,
                  ent = apply(assay(Entropy_SE[Entropy_SE@rowRanges$loc %in% locs,]),2,mean,na.rm=T),
                  group=sampleinfo$group)
  ProEnt<-rbind(ProEnt,tmp)
}
save(ProEnt,file="PromoterCpGOE.RData")
ggplot(ProEnt,aes(fill=group,x=group,y=ent))+geom_boxplot()+facet_wrap(~type,nrow = 1)+
  geom_signif(comparisons = list(c(1,2),c(2,3)),test = ks.test,map_signif_level = c("***"=0.001,"**"=0.01,"*"=0.05))+
  theme_classic()+
  scale_fill_manual(values = c("Long-lived"="red","Elder Control"="darkgoldenrod","Younger Control"="skyblue"))

##### Funcs Anno #####
annos<-c(annos<-c('hg19_cpgs', 'hg19_genes_intergenic',
                  'hg19_lncrna_gencode','hg19_enhancers_fantom'))
annotations = build_annotations(genome = 'hg19', annotations = annos)
anno2<-annotations[!is.na(annotations$gene_id),]

Ent_anno<-annotate_regions(regions = rowRanges(Entropy_SE),
                           annotations = annotations,
                           ignore.strand = TRUE,
                           quiet = FALSE)
Ent_anno$annot<-Ent_anno$annot$type
Ent_anno<-unique(data.frame(seqnames=seqnames(Ent_anno),
                            start=start(Ent_anno),
                            end=end(Ent_anno),
                            type=Ent_anno$annot,loc=Ent_anno$loc)) %>% as.data.frame %>% as("GRanges")

FuncEnt<-data.frame()
for(i in 1:7){
  types<-unique(Ent_anno$type)[i]
  locs<-intersect(Ent_anno[Ent_anno$type == types,]$loc, Entropy_SE@rowRanges$loc)
  tmp<-data.frame(type=types,
                  ent = apply(assay(Entropy_SE[Entropy_SE@rowRanges$loc %in% locs,]),2,mean,na.rm=T),
                  group=sampleinfo$group)
  FuncEnt<-rbind(FuncEnt,tmp)
}
save(FuncEnt,file="Funcgene.RData")
ggplot(FuncEnt,aes(fill=group,x=group,y=ent))+geom_boxplot()+facet_wrap(~type,nrow = 1)+
  geom_signif(comparisons = list(c(1,2),c(2,3)))+theme_classic()+
  scale_fill_manual(values = c("Long-lived"="red","Elder Control"="darkgoldenrod","Younger Control"="skyblue"))



##### Hi-C#####
library(rtracklayer)
compAB<-import("./HiC/GSE63525_GM12878_subcompartments.bed")
compAB[grep(compAB$name,pattern = "A")]$name<-"A"
compAB[grep(compAB$name,pattern = "B")]$name<-"B"

compA<-findOverlaps(rowRanges(Entropy_SE),compAB[compAB$score>0])
A_each<-apply(assay(Entropy_SE[unique(compA@from),]),2,mean,na.rm=T)
compB<-findOverlaps(rowRanges(Entropy_SE),compAB[compAB$score<0,])
B_each<-apply(assay(Entropy_SE[unique(compB@from),]),2,mean,na.rm=T)

sampleinfo$compA<-A_each
sampleinfo$compB<-B_each

library(ggsignif)
ggplot(sampleinfo,aes(x=group,y=A_each,fill=group))+geom_boxplot()+
  geom_signif(comparisons = list(c("Long-lived","Elder Control"),c("Elder Control","Younger Control")),test = ks.test)+
  scale_fill_manual(values = c("Younger Control"="skyblue1", "Elder Control" = "gold4", "Long-lived"="red"))+
  theme_classic()+theme(text=element_text(size=18))+ylab(label = "Compartment A")
ggplot(data.frame(age=sampleinfo$Age,GCA=A_each),aes(x=age,y=GCA))+geom_point(aes(color=age<=80))+
  scale_color_manual(values = c("TRUE"="blue1", 
                                "FALSE" ="red"))+geom_smooth(method = "lm",aes(color=age<=80))+
  theme_classic()
cor.test(sampleinfo[sampleinfo$Age<90,]$Age,
         sampleinfo[sampleinfo$Age<90,]$compA)

ggplot(sampleinfo,aes(x=group,y=B_each,fill=group))+geom_boxplot()+
  geom_signif(comparisons = list(c("Long-lived","Elder Control"),c("Elder Control","Younger Control")),test = ks.test)+
  scale_fill_manual(values = c("Younger Control"="skyblue1", "Elder Control" = "gold4", "Long-lived"="red"))+
  theme_classic()+theme(text=element_text(size=18))+ylab(label = "Compartment B")
ggplot(data.frame(age=sampleinfo$Age,GCB=B_each),aes(x=age,y=GCB))+geom_point(aes(color=age<=80))+
  scale_color_manual(values = c("TRUE"="blue1", 
                                "FALSE" ="red"))+
  theme_classic()
save(sampleinfo,file="HiC/HiC_result.RData")
cor.test(sampleinfo[sampleinfo$Age<90,]$Age,
         sampleinfo[sampleinfo$Age<90,]$compB)

######ATAC accissible#####
library(rtracklayer)
library(MethyAge2)
library(stringr)
dir.create("ATAC")
setwd("ATAC/") #Then upload GSE47753
GM12878_ATAC<-read.table("/data/user_data/wanght/HNLong/DMEAS/Another_81/Int3CG/ATAC/GSE47753_GM12878_ATACseq_50k_AllReps_ZINBA_pp08.bed",
                         header = F)
GM12878_ATAC<-GM12878_ATAC[,c(1,2,3)]
colnames(GM12878_ATAC)[1:3]<-c("seqnames","start","end")
GM12878_ATAC$strand<-"*";GM12878_ATAC$type<-"ATAC_acc"
GM12878_ATAC<-as(GM12878_ATAC,"GRanges")

ATAC<-findOverlaps(rowRanges(Ent_rm),GM12878_ATAC)
ATAC_each<-apply(assay(Ent_rm[unique(ATAC@from),]),2,mean,na.rm=T)
non_each<-apply(assay(Ent_rm[-unique(compA@from),]),2,mean,na.rm=T)

sampleinfo$ATAC<-ATAC_each
ggplot(sampleinfo,aes(x=group,y=ATAC,fill=group))+geom_boxplot()+
  geom_signif(comparisons = list(c("Long-lived","Elder Control"),c("Elder Control","Younger Control")),test = ks.test)+
  scale_fill_manual(values = c("Younger Control"="skyblue1", "Elder Control" = "gold4", "Long-lived"="red"))+
  theme_classic()+theme(text=element_text(size=18))
ggplot(data.frame(age=sampleinfo$Age,ATAC=sampleinfo$ATAC),aes(x=age,y=ATAC))+geom_point(aes(color=age<=80))+
  scale_color_manual(values = c("TRUE"="blue1", 
                                "FALSE" ="red"))+
  theme_classic()
cor.test(sampleinfo[sampleinfo$Age<90,]$Age,
         sampleinfo[sampleinfo$Age<90,]$ATAC)

sampleinfo$non<-non_each
ggplot(sampleinfo,aes(x=group,y=non,fill=group))+geom_boxplot()+
  geom_signif(comparisons = list(c("Long-lived","Elder Control"),c("Elder Control","Younger Control")))+
  scale_fill_manual(values = c("Younger Control"="skyblue1", "Elder Control" = "gold4", "Long-lived"="red"))+
  theme_classic()+theme(text=element_text(size=18))
cor.test(sampleinfo[sampleinfo$Age<90,]$Age,
         sampleinfo[sampleinfo$Age<90,]$non)

save(sampleinfo,file="ATAC_resualt.RData")
save(sampleinfo,file="../HiC/HiC_resualt.RData")
#####OGEE#####
setwd("../")
#load("OGEE/9606_unique_OGEE.RData")
OGEE<-read.csv("OGEE/OGEEv3_blood.csv")

annos<-c('hg19_basicgenes')
annotations = build_annotations(genome = 'hg19', annotations = annos)
anno2<-annotations[!is.na(annotations$gene_id),]

Ent_anno<-annotate_regions(regions = rowRanges(Entropy_SE),
                           annotations = annotations,
                           ignore.strand = TRUE,
                           quiet = FALSE)
Ent_anno$type<-Ent_anno$annot$type;Ent_anno$symbol<-Ent_anno$annot$symbol
Ent_anno$Enss<-"NE";Ent_anno[Ent_anno$symbol %in% OGEE$genes,]$Enss<-"E"
Ent_anno$annot<-"Ns"
Ent_anno<-unique(data.frame(seqnames=seqnames(Ent_anno),
                            start=start(Ent_anno),
                            end=end(Ent_anno),
                            type=Ent_anno$type,loc=Ent_anno$loc,Enss=Ent_anno$Enss)) %>% as.data.frame %>% as("GRanges")

EP<-Ent_anno[Ent_anno$type == "hg19_genes_promoters",]
EP_ent<-data.frame()
for(i in 1:2){
  types<-c("E","NE")[i]
  locs<-intersect(EP[EP$Enss == types,]$loc, Entropy_SE@rowRanges$loc)
  tmp<-data.frame(type=types,
                  ent = apply(assay(Entropy_SE[Entropy_SE@rowRanges$loc %in% locs,]),2,mean,na.rm=T),
                  group=sampleinfo$group)
  EP_ent<-rbind(EP_ent,tmp)
}
ggplot(EP_ent,aes(fill=group,x=group,y=ent))+geom_boxplot()+facet_wrap(~type,nrow = 1)+
  geom_signif(comparisons = list(c(1,2),c(2,3)),test = ks.test,map_signif_level = c("***"=0.001,"**"=0.01,"*"=0.05))+theme_classic()+
  scale_fill_manual(values = c("Long-lived"="red","Elder Control"="darkgoldenrod","Younger Control"="skyblue"))


#Ent_anno1[Ent_anno1$ENSEMBL %in% OGEE[OGEE$essential == "E",]$locus,]$Enss<-"E"



saveRDS(gene_each,file="./OGEE/gene_res.RDS")
