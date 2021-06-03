setwd("/home/owht/KIZ/data/Entropy/Result/DMEAS/DifferenEntropy/EWAS")
setwd("/media/owht/YNWA-wht//KIZ/data/Entropy/Result/DMEAS/DifferenEntropy/EWAS/")
library(ggplot2)
library(rtracklayer)
library(magrittr)

#preparation
times<-2000 #5000 times sampling
size<-18282 # 18282 low-DECSs
load("../DER_backg_10417699.RData")
samplings<-matrix(nrow=times,
                  ncol=size,data=0)
for(i in 1:times){
  samplings[i,]<-sample(x = 1:length(backg),replace = F,size = size)
}

#load data
load("/home/owht/Data/KIZ/data/AP/Dataset/probe_Liftover_result.RData")
load("/home/owht/Data/KIZ/data/AP/Dataset/Clean_33disease.RData")
cginfo_gr<-cginfo_gr[unique(associa$Probe.id)]
associa<-associa[,c(1,2,7,8,12)]

#shuffle
results<-matrix(nrow = nrow(disease),
                ncol=times,data=0)
rownames(results)<-disease$disease
for(d in 1:nrow(disease)){
  ass<-associa[associa$disease == disease[d,]$disease,]$Probe.id
  tmp_ass<-cginfo_gr[ass]
  cat(paste0("disease is ",disease[d,]$disease,"\n"))
  for(t in 1:times){
    tmp_back<-backg[samplings[t,]]
    results[d,t]<-intersect(tmp_ass,tmp_back) %>% length()
    if(t%%200 == 0){cat(paste0("200!","\n"))}
  }
}

#load the low DECSs
low<-import.bed("../bed_fasta/low_DECSs.bed")

real<-rep(0,33)
for(d in 1:nrow(disease)){
  ass<-associa[associa$disease == disease[d,]$disease,]$Probe.id
  tmp_ass<-cginfo_gr[ass]
  real[d]<-intersect(tmp_ass,low) %>% length()
}

save(samplings,results,real,file="sampling_shuffle_33diseases_low.RData")

#
load("sampling_shuffle_33diseases_low.RData")
par(mfrow = c(4,4))
for(d in 1:33){
  hist(results[d,],main=disease[d,]$disease,xlim = c(0,max(real[d],max(results[d,]))+1),
       xlab = "Simulated hits")
  abline(v = real[d])
}


result2<-data.frame(disease=disease$disease,
                    pvalue = apply(1-(real > results),
                                   MARGIN = 1,
                                   FUN = function(x){return(sum(x))})/2000,
                    nums = real)
result2<-merge(result2,disea_probes,by.x = 1,by.y = 1)
result2$ratio<-(result2$nums)/(result2$Freq)
result2$disease<-factor(result2$disease,
                        levels = result2$disease[order(result2$ratio,decreasing = T)])
ggplot(result2,aes(x = disease,y=ratio,fill = nums))+
  geom_bar(stat = "identity")+theme_classic()+
  scale_fill_gradient(low = "grey",high="brown")

save(result2,file="low_shuffle_result.RData")

#
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
load("low_shuffle_result.RData")
dir.create("disease")
ath<-intersect(low,
               cginfo_gr[associa[associa$disease == "atherosclerosis",]$Probe.id]) #atherosclerosis
ath<-annotatePeak(peak = ath,tssRegion = c(-2500,500),
                  TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                  annoDb = "org.Hs.eg.db")

#Circos
result2<-merge(result2,disease,by.x=1,by.y=2)
rownames(result2)<-result2$disease
result2<-result2[result2$pvalue<0.05,]
data("UCSC.HG19.Human.CytoBandIdeogram")

RCircos.Set.Core.Components(cyto.info=UCSC.HG19.Human.CytoBandIdeogram,
                            tracks.inside=2,tracks.outside=0,chr.exclude = c("chrX","chrY"))
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
xx<-data.frame()
for (d in 1:12){
  ath<-intersect(low,
                 cginfo_gr[associa[associa$disease == result2[d,]$disease,]$Probe.id]) 
  tmp<-table(ath@seqnames) %>% as.data.frame(row.names = 1) %>% t
  xx<-rbind(xx,tmp)
}
rownames(xx)<-result2$disease
pheatmap(xx,annotation_row = result2[,5:6])

#
associa$hit<-"no"
hi<-findOverlaps(cginfo_gr,low)
associa[associa$Probe.id %in% 
          names(cginfo_gr[hi@from]),]$hit<-"yes"
table(associa$Correlation, associa$hit) %>% fisher.test()
# Fisher's Exact Test for Count Data
# 
# data:  .
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.3687849 0.5447039
# sample estimates:
# odds ratio 
#  0.4490267 
ned<-table(associa$Correlation, associa$hit) %>% melt
ggplot(ned,aes(x=Var1,y=value,fill=Var2))+geom_bar(stat = "identity")+scale_y_log10()+
  coord_flip()+theme_classic()
associa$CHROM<-cginfo_gr[associa$Probe.id]@seqnames %>% as.character()
associa$POS<-cginfo_gr[associa$Probe.id]@ranges@start
associa$P.value<-as.numeric(associa$P.value)
associa$CHROM<-str_split(associa$CHROM,pattern = "chr",simplify = T)[,2]
associa$CHROM<-as.numeric(associa$CHROM)
manhattan(na.omit(associa[associa$hit == "yes",]),chr = "CHROM",bp = "POS",p="P.value",
          )
associa[associa$hit == "yes",] -> EWAS_overlap_lowDECS
write.csv(EWAS_overlap_lowDECS,file="../../../../FigTable/tables/EWAS.csv")
