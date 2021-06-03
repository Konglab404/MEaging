library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(stringr)
library(magrittr)
setwd("/home/owht/KIZ/data/Entropy/Result/DMEAS/DifferenEntropy")
#data("tagMatrixList")

load("limma_DER_result_01_p001.RData")
DER_result_01_001p->DER_result;rm(DER_result_01_001p)
DER_result$seqnames<-str_split(rownames(DER_result),
                               pattern = ":|_",simplify=T)[,1]
DER_result$start<-str_split(rownames(DER_result),
                               pattern = ":|_",simplify=T)[,2]
DER_result$end<-str_split(rownames(DER_result),
                            pattern = ":|_",simplify=T)[,4]
DER_result<-DER_result[DER_result$start!="",]
DER_all<-as(DER_result,"GRanges")
DER_high<-as(DER_result[DER_result$logFC>0,],"GRanges")
DER_low<-as(DER_result[DER_result$logFC<0,],"GRanges")

DER_result$width<-as.numeric(DER_result$end) - as.numeric(DER_result$start)
DER_result$mid<-(as.numeric(DER_result$end) + as.numeric(DER_result$start)) %/%2

DER_all_anno<-annotatePeak(peak = DER_all,tssRegion = c(-2500,500),
                            TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                            annoDb = "org.Hs.eg.db")
upsetplot(DER_all_anno)
vennpie(DER_all_anno)
tagMatrix<-tagMatrixList[[4]]

######
DER_stat<-data.frame(x=rep(DER_all_anno@detailGenomicAnnotation %>% colnames,2),
                     type=c(rep("high",9),rep("low",9)))

plotAvgProf(tagMatrix = tagMatrix,xlim=c)
plotDistToTSS(DER_all_anno)
######Distance
Dist_TSS<-data.frame(High=cut(DER_high_anno_gr$distanceToTSS,breaks = c(-Inf,-500000,-50000,-5000,0,5000,50000,500000,Inf)) %>% table %>% as.numeric(),
                          Low=cut(DER_low_anno_gr$distanceToTSS,breaks = c(-Inf,-500000,-50000,-5000,0,5000,50000,500000,Inf)) %>% table %>% as.numeric())
Dist_TSS[,1]/colSums(Dist_TSS)[1]->Dist_TSS[,1]
Dist_TSS[,2]/colSums(Dist_TSS)[2]->Dist_TSS[,2]
library(reshape2)
Dist_TSS_melt<-melt(Dist_TSS)
Dist_TSS_melt$dist<-rep(paste0("Dist",1:8),2)
ggplot(Dist_TSS_melt,aes(x=dist,y=value,fill=variable))+
  geom_bar(stat="identity",position = "dodge")+
  theme_classic()+ylim(0,0.4)
###BED format
library(rtracklayer)
con=file("DER_high.bed")
export(object = DER_high,con = con,format = "bed")
con=file("DER_low.bed")
export(object = DER_low,con = con,format = "bed")

###Anno
DER_high_anno<-annotatePeak(peak = DER_high,tssRegion = c(-2500,500),
                           TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                           annoDb = "org.Hs.eg.db")
vennpie(DER_high_anno)
DER_low_anno<-annotatePeak(peak = DER_low,tssRegion = c(-2500,500),
                            TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                            annoDb = "org.Hs.eg.db")
vennpie(DER_low_anno)

DER_high_anno@anno->DER_high_anno_gr
DER_low_anno@anno->DER_low_anno_gr
DER_high_anno_gr[DER_high_anno_gr$annotation!="Distal Intergenic",]$SYMBOL %>% 
  unique
DER_all<-data.frame(genes=c(DER_high_anno_gr$SYMBOL,DER_low_anno_gr$SYMBOL),
                    diff=c(DER_high_anno_gr$logFC,DER_low_anno_gr$logFC),
                    anno=c(DER_high_anno_gr$annotation,DER_low_anno_gr$annotation),
                    methdiff=c(DER_high_anno_gr$Meth.diff,DER_low_anno_gr$Meth.diff))
DER_table<-DER_all[-grep(DER_all$anno,pattern = "D"),]$genes %>% 
  table %>% as.data.frame()
DER_tmp<-DER_all[-grep(DER_all$anno,pattern = "D"),]
DER_tmp<-aggregate(x=DER_tmp$diff,by=list(DER_tmp$genes),
                   FUN=function(x){return(sum(x>0))})
DER_table<-merge(x=DER_table,y=DER_tmp,by.x=1,by.y=1);rm(DER_tmp)
DER_table$highratio<-DER_table$x/DER_table$Freq

DER_all$anno<-str_replace(string = DER_all$anno,
                          pattern = "\\(.*\\)",replacement = "")
DER_stat<-aggregate(DER_all$diff,by=list(DER_all$anno),
                    FUN=function(x){return(sum(x>0))})
DER_tmp<-table(DER_all$anno) %>% as.data.frame()
DER_stat<-merge(x=DER_tmp,y=DER_stat,by.x="Var1",by.y="Group.1")
DER_stat$highratio<-DER_stat$x/DER_stat$Freq
DER_stat$Var1<-factor(DER_stat$Var1,levels=DER_stat$Var1[order(DER_stat$highratio)])

DER_table$low<-DER_table$Freq - DER_table$x
write.table(DER_all[DER_all$anno == "Promoter " &
                      DER_all$diff>0,]$genes %>% unique,
            col.names = F,row.names = F,quote=F,
            file="High_Promoter")
write.table(DER_all[DER_all$anno == "Promoter " &
                      DER_all$diff<0,]$genes %>% unique,
            col.names = F,row.names = F,quote=F,
            file="Low_Promoter")
write.table(DER_all[DER_all$methdiff>0 & 
                      DER_all$diff<0.1,]$genes %>% unique,
            col.names = F,row.names = F,quote = F,
            file="Hyper_low_DER_genes")

write.table(DER_table[,c(1,2)],
            file="DER_all_gene.rnk",quote=F,col.names = F,
            row.names = F,sep="\t")
write.table(DER_table[,c(1,3)],
            file="DER_high_gene.rnk",quote=F,col.names = F,
            row.names = F,sep="\t")
write.table(DER_table[,c(1,5)],
            file="DER_low_gene.rnk",quote=F,col.names = F,
            row.names = F,sep="\t")
#Density
library(ggplot2)
DER_entCEN<-data.frame(x=c(DER_result_01_001p[DER_result_01_001p$type == "High",]$entCEN,
                           DER_result_01_001p[DER_result_01_001p$type == "Low",]$entCEN),
                       type=c(rep("High",101853),rep("Low",18282)))
ggplot(DER_entCEN,aes(x,colour=type,fill=type))+
  geom_histogram(color="black",aes(y=0.05*..density..),
                 position = "identity",binwidth=0.05,alpha=0.4)+
  stat_density(geom = 'line',position = 'identity',aes(y=0.05* ..density..))+theme_classic()+
  theme(text=element_text(size=18))

DER_entSP<-data.frame(x=c(DER_result_01_001p[DER_result_01_001p$type == "High",]$entSP,
                           DER_result_01_001p[DER_result_01_001p$type == "Low",]$entSP),
                       type=c(rep("High",101853),rep("Low",18282)))
ggplot(DER_entSP,aes(x,colour=type,fill=type))+
  geom_histogram(color="black",aes(y=0.05*..density..),
                 position = "identity",binwidth=0.05,alpha=0.4)+
  stat_density(geom = 'line',position = 'identity',aes(y=0.05* ..density..))+theme_classic()+
  theme(text=element_text(size=18))


####For LOLA
High<-DER_high_anno_gr[,1]
High$name<-names(High);High$score<-0;High<-High[,2:3]
names(High)<-NULL
Low<-DER_low_anno_gr[,1]
Low$name<-names(Low);Low$score<-0;Low<-Low[,2:3]
names(Low)<-NULL
#####################ggplot
ggplot(DER_stat,aes(x=Var1,y=highratio))+geom_bar(stat = "identity",
                                                position = position_dodge())+
  theme_classic()+ylim(c(0,1))

load("limma_methdiff_DER_result_p001.RData")
Meth.diff2<-Meth.diff2[rownames(Meth.diff2) %in% 
                         rownames(DER_result),]

####Prepare for liftover
tmp<-data.frame(chrom=DER_result[DER_result$logFC<0,]$seqnames,
                start=DER_result[DER_result$logFC<0,]$mid - 55,
                end=DER_result[DER_result$logFC<0,]$mid + 55)
write.table(tmp,file="./Liftover/DER_low_55bp.bed",col.names = F,row.names = F,
            quote=F,sep="\t")
tmp<-data.frame(chrom=DER_result[DER_result$logFC>0,]$seqnames,
                start=DER_result[DER_result$logFC>0,]$mid - 55,
                end=DER_result[DER_result$logFC>0,]$mid + 55)
write.table(tmp,file="./Liftover/DER_high_55bp.bed",col.names = F,row.names = F,
            quote=F,sep="\t")
LifO<-read.table("Liftover/Liftover_result")
LifO$Species<-rownames(LifO)
LifO<-melt(LifO)
LifO$Species<-rep(c("B5","B4","B3","B2","B1","B0"),2)
LifO$Species<-factor(LifO$Species)
ggplot(LifO,aes(x=Species,y=value,color=variable,group=variable))+
  geom_line() + theme_classic()

#####For Circos
high<-as.data.frame(DER_high)
high$seqnames<-str_replace(high$seqnames,pattern = "chr",
                           replacement = "hs")
write.table(high[,c(1,2,3,6)],quote=F,sep = "\t",
            row.names = F,col.names = F,
            file = "../../Circos/data/High_DER")
low<-as.data.frame(DER_low)
low$seqnames<-str_replace(low$seqnames,pattern = "chr",
                           replacement = "hs")
write.table(low[,c(1,2,3,6)],quote=F,sep = "\t",
            row.names = F,col.names = F,
            file = "../../Circos/data/low_DER")
10185

load("DER_backg_10417699.RData")
backg$loc<-paste0(backg@seqnames,"_",backg@ranges@start)
DER_high_anno_gr$loc<-paste0(DER_high_anno_gr@seqnames,"_",DER_high_anno_gr@ranges@start)
DER_low_anno_gr$loc<-paste0(DER_low_anno_gr@seqnames,"_",DER_low_anno_gr@ranges@start)
backg_real<-backg[!backg$loc %in% c(DER_high_anno_gr$loc, DER_low_anno_gr),]

boxplot(backg_real@ranges@width,
        DER_high_anno_gr@ranges@width,
        DER_low_anno_gr@ranges@width,
        outline=F)

################

#GSEA preRank
DER_all<-data.frame(genes=c(DER_high_anno_gr$SYMBOL,DER_low_anno_gr$SYMBOL),
                    ens=c(DER_high_anno_gr$ENSEMBL,DER_low_anno_gr$ENSEMBL),
                    diff=c(DER_high_anno_gr$logFC,DER_low_anno_gr$logFC),
                    anno=c(DER_high_anno_gr$annotation,DER_low_anno_gr$annotation),
                    methdiff=c(DER_high_anno_gr$Meth.diff,DER_low_anno_gr$Meth.diff))
DER_table<-DER_all[-grep(DER_all$anno,pattern = "D"),]$genes %>% 
  table %>% as.data.frame()
DER_tmp<-DER_all[-grep(DER_all$anno,pattern = "D"),]
DER_tmp<-aggregate(x=DER_tmp$diff,by=list(DER_tmp$genes),
                   FUN=function(x){return(sum(x>0))})
DER_table<-merge(x=DER_table,y=DER_tmp,by.x=1,by.y=1);rm(DER_tmp)
DER_table$low<-DER_table$Freq - DER_table$x

write.table(DER_table[,c(1,3)],row.names = F,col.names = F,sep = " ",
            quote=F,file="GSEA/perRank_DER_high.rnk")
write.table(DER_table[,c(1,5)],row.names = F,col.names = F,sep = " ",
            quote=F,file="GSEA/perRank_DER_low.rnk")

colnames(DER_table)[c(1,3)]<-c("ens","high")
DER_table$ens<-as.character(DER_table$ens)

##########OGEE 4 DERs
load("~/KIZ/data/from_web/OGEE/9606_unique_OGEE.RData")

Egene<-OGEE[OGEE$essential == "E",]$locus
NEgene<-OGEE[OGEE$essential == "NE",]$locus
DER_table$essential<-"NE"
DER_table[DER_table$ens %in% Egene,]$essential<-"E"
tmp<-matrix(c(sum(DER_table[DER_table$essential == "E",]$high),
              sum(DER_table[DER_table$essential == "E",]$low),
              sum(DER_table[DER_table$essential == "NE",]$high),
              sum(DER_table[DER_table$essential == "NE",]$low)),ncol=2)
fisher.test(tmp) 
# #
# p-value = 1.193e-07
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.2943381 0.5652579
# sample estimates:
#   odds ratio 
# 0.4061337 
# #
save(tmp,DER_table,file="OGEE/OGEE.RData")
save(DER_all_anno,file="DER_anno_ChIPseeker.RData")


par(mfrow=c(1,2))
boxplot(DER_table[DER_table$ens %in% Egene,]$high,DER_table[DER_table$ens %in% NEgene,]$high,
        outline=F)
boxplot(DER_table[DER_table$ens %in% Egene,]$low,DER_table[DER_table$ens %in% NEgene,]$low,
        outline=F)
dev.off()



############ DER 4 nc aging genes

load("DER_anno_ChIPseeker.RData")
str_replace_all(DER_all$anno,pattern ="\\(.*\\)",replacement = "")->DER_all$anno

DER_all_table<-DER_all[-grep(DER_all$anno,pattern = "D"),]$genes %>% 
  table %>% as.data.frame()
DER_tmp<-DER_all[-grep(DER_all$anno,pattern = "D"),]
DER_tmp<-aggregate(x=DER_tmp$diff,by=list(DER_tmp$genes),
                   FUN=function(x){return(sum(x>0))})
DER_all_table<-merge(x=DER_all_table,y=DER_tmp,by.x=1,by.y=1);rm(DER_tmp)
DER_all_table$low<-DER_table$Freq - DER_all_table$x

DER_pro_table<-DER_all[DER_all$anno == "Promoter ",]$genes %>% 
  table %>% as.data.frame()
DER_tmp<-DER_all[DER_all$anno == "Promoter ",]
DER_tmp<-aggregate(x=DER_tmp$diff,by=list(DER_tmp$genes),
                   FUN=function(x){return(sum(x>0))})
DER_pro_table<-merge(x=DER_pro_table,y=DER_tmp,by.x=1,by.y=1);rm(DER_tmp)
DER_pro_table$low<-DER_pro_table$Freq - DER_pro_table$x

DER_gb_table<-DER_all[-grep(DER_all$anno,pattern = "D|P"),]$genes %>% 
  table %>% as.data.frame()
DER_tmp<-DER_all[-grep(DER_all$anno,pattern = "D|P"),]
DER_tmp<-aggregate(x=DER_tmp$diff,by=list(DER_tmp$genes),
                   FUN=function(x){return(sum(x>0))})
DER_gb_table<-merge(x=DER_gb_table,y=DER_tmp,by.x=1,by.y=1);rm(DER_tmp)
DER_gb_table$low<-DER_gb_table$Freq - DER_gb_table$x

age_gene<-read.table("/home/owht/KIZ/data/from_web/NC_2015_geneexp_age/nc9570_age_exp",
                     header = T,stringsAsFactors = F)


age_gene_der_all<-merge(DER_all_table,age_gene,by.x=1,by.y="NEW.Gene.ID",all.y=T)
age_gene_der_all[is.na(age_gene_der_all)]<-0
all_melt<-melt(age_gene_der_all,measure.vars = c(3,4))
ggplot(all_melt,aes(x=variable,fill=Direction,y=value))+geom_boxplot()
t.test(age_gene_der_all[age_gene_der_all$Direction == "+",]$Freq,
       age_gene_der_all[age_gene_der_all$Direction=="-",]$Freq)

age_gene_der_pro<-merge(DER_pro_table,age_gene,by.x=1,by.y="NEW.Gene.ID")
pro_melt<-melt(age_gene_der_pro,measure.vars = c(3,4))
ggplot(pro_melt,aes(x=variable,fill=Direction,y=value))+geom_boxplot()
t.test(age_gene_der_pro[age_gene_der_pro$Direction == "+",]$x,
       age_gene_der_pro[age_gene_der_pro$Direction=="-",]$x)
t.test(age_gene_der_pro[age_gene_der_pro$Direction == "+",]$low,
       age_gene_der_pro[age_gene_der_pro$Direction=="-",]$low)

###########################
DER_all_anno_gr<-DER_all_anno@anno
DER_high_gr<-DER_all_anno_gr[DER_all_anno_gr$en.diff>0,]
DER_low_gr<-DER_all_anno_gr[DER_all_anno_gr$en.diff<0,]

library(rtracklayer)
library(MethyAge2)
chromHMM<-import.bed("~/KIZ/data/Basic/RefNome/Human/bed/GM12878/GM12878_chromHMM.bed")
insul<-chromHMM[grep(chromHMM$name,pattern = "Insu"),]
rm(chromHMM)

low_insulator<-DMSAnnotate(DER_low_gr,insul)
write.csv(low_insulator,file="Metascape/DER_low_insulator.csv")












