setwd("~/data/Entropy/Result/")
library(stringr)
library(magrittr)
library(data.table)
library(SummarizedExperiment)
library(limma)
library(EpiDISH)
library(reshape2)

sampleinfo<-read.table("../Sampleinfo133")

load("../../DNAmAge/EpiDISH_176.RData")
cf<-EpiDISH_est$estF[str_replace(colnames(Entropy_SE),pattern = "Entropy",replacement = "sample_"),]
cf<-as.data.frame(cf)
cf$type<-Entropy_SE$group
cf$Age<-Entropy_SE$Age
cf_EC<-cf[cf$type %in% c("Long-lived","Elder Control"),]
cf_EC$type<-as.character(cf_EC$type)
cf_EC<-melt(cf_EC[,c(1:6,8)],measure.vars=1:6,id.vars=7,factorsAsStrings = TRUE)
cf_EC$variable<-as.character(cf_EC$variable)
#levels(cf_EC$variable)<-c
ggplot(cf_EC,aes(x=variable,fill=type,y=value))+geom_boxplot()+
  stat_compare_means(aes(group=type),label = "p.signif",method = "kruskal.test")+
  theme_bw()+theme(text=element_text(size=18))


Entropy_SE<-readRDS("Ent_SE_133_0.5_for3cats_rb.RDS")
DER_CEN_SE<-Entropy_SE[rowRanges(Entropy_SE)$loc %in% DER_CEN$loc, Entropy_SE$group %in% c("Long-lived","Elder Control")]
cf_EC<-cf[str_replace(colnames(DER_CEN_SE),pattern = "Entropy",replacement = "sample_"),]

tmpassay<-assay(DER_CEN_SE)

library(impute)
tmpassay<-impute.knn(tmpassay ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500)
tmpassay$data[tmpassay$data<0]<-0
tmpassay$data[tmpassay$data>0]<-1

fake<-CellDMC(tmpassay$data,as.numeric(DER_CEN_SE$group)-2,as.matrix(cf_EC[,1:6]),adjPMethod = "none")
x = data.frame(type=names(apply(Cell_DER$dmct,2,function(x){return(sum(x==-1))})),
               num = apply(Cell_DER$dmct,2,function(x){return(sum(x==-1)+1)})) 
ggplot(x,aes(x=type,y=num))+geom_bar(stat = "identity")+scale_y_log10(breaks = c(1:10,10,20,30,40,50,60,70,80,90,100,
                                                                                 200,300,400,500,600,700,800,900,1000,
                                                                                 2000,3000,4000,5000,6000,7000,8000,9000,10000))+
  theme_bw()
####################
DER_CEN[Cell_DER$dmct[,2] == -1,]->B_DER
DER_CEN[Cell_DER$dmct[,3] == -1,]->NK_DER
DER_CEN[Cell_DER$dmct[,4] == -1,]->CD4T_DER
DER_CEN[Cell_DER$dmct[,5] == -1,]->CD8T_DER
DER_CEN[Cell_DER$dmct[,6] == -1,]->Mono_DER
DER_CEN[Cell_DER$dmct[,7] == -1,]->Neutro_DER
save(B_DER,NK_DER,CD4T_DER,CD8T_DER,Mono_DER,Neutro_DER,file="Cell_DERs.RData")
library(annotatr)
annos<-c('hg19_cpgs')
annotations = build_annotations(genome = 'hg19', annotations = annos)
CpG_res<-data.frame()
der_names<-c("B_DER","NK_DER","CD4T_DER","CD8T_DER","Mono_DER","Neutro_DER")
for(i in 1:6){
  tmp<-annotate_regions(regions = list(B_DER,NK_DER,CD4T_DER,CD8T_DER,Mono_DER,Neutro_DER)[[i]],
                        annotations = annotations,ignore.strand = TRUE,quiet = FALSE)
  tmp<-data.frame(table(tmp$annot$type),type = der_names[i])
  CpG_res<-rbind(CpG_res,tmp)
}
ggplot(CpG_res,aes(fill=Var1,x=type,y=Freq))+geom_col(position = "fill")+coord_flip()+theme_classic()

annos<-c('hg19_basicgenes')
annotations = build_annotations(genome = 'hg19', annotations = annos)
Genic_res<-data.frame()
der_names<-c("B_DER","NK_DER","CD4T_DER","CD8T_DER","Mono_DER","Neutro_DER")
for(i in 1:6){
  tmp<-annotate_regions(regions = list(B_DER,NK_DER,CD4T_DER,CD8T_DER,Mono_DER,Neutro_DER)[[i]],
                        annotations = annotations,ignore.strand = TRUE,quiet = FALSE)$annot %>% 
    as.data.frame(stringsAsFactors=F)
  tmp<-unique(tmp[,c("seqnames","start","type","symbol")])
  tmp<-rbind(data.frame(table(tmp$type),type = der_names[i]),
             data.frame(Var1="Distal",
                        Freq=(length(list(B_DER,NK_DER,CD4T_DER,CD8T_DER,Mono_DER,Neutro_DER)[[i]])-
                             (annotate_regions(regions = list(B_DER,NK_DER,CD4T_DER,CD8T_DER,Mono_DER,Neutro_DER)[[i]],
                                           annotations = annotations,ignore.strand = TRUE,quiet = FALSE)$loc %>% unique() %>% length)),
                        type=der_names[i]
             )
  )
  Genic_res<-rbind(Genic_res,tmp)
}
ggplot(Genic_res,aes(fill=Var1,x=type,y=Freq))+geom_col(position = "fill")+coord_flip()+theme_classic()

save(CpG_res,Genic_res,file="./CellDMR/Annotations_CellDER.RData")
############################

for(i in 1:6){
  tmp<-annotate_regions(regions = list(B_DER,NK_DER,CD4T_DER,CD8T_DER,Mono_DER,Neutro_DER)[[i]],
                        annotations = annotations,ignore.strand = TRUE,quiet = FALSE)
  tmp$annot[tmp$annot$type == "hg19_genes_promoters",]$symbol %>% table %>% as.data.frame(stringsAsFactors=F)->tmp_pro
}

neutro<-annotate_regions(regions = Neutro_DER,
                 annotations = annotations,ignore.strand = TRUE,quiet = FALSE)
neutro2<-as.data.frame(neutro$annot);neutro2$loc<-neutro$loc
neutro2_pro<-unique(neutro2[neutro2$type == "hg19_genes_promoters",c("seqnames","start","symbol")])
neutro2_pro<-table(neutro2_pro$symbol) %>% as.data.frame(stringsAsFactors=F)
neutro2_pro[neutro2_pro$Freq>=3,]$Var1 %>% unique %>% cat

DER_CEN_grlist<-list(B_DER,NK_DER,CD4T_DER,CD8T_DER,Mono_DER,Neutro_DER)
library(LOLA)
histone<-loadRegionDB("/data/user_data/wanght/HNLong3/LOLA_base/chromhmm") 


######## partitioned LDSC #####
annos<-c('hg19_basicgenes')
annotations = build_annotations(genome = 'hg19', annotations = annos)
anno2<-annotations[!is.na(annotations$gene_id),]
for(i in 1:6){
  tmp<-annotate_regions(regions = list(B_DER,NK_DER,CD4T_DER,CD8T_DER,Mono_DER,Neutro_DER)[[i]],
                        annotations = annotations,ignore.strand = TRUE,quiet = FALSE)$annot %>% 
    as.data.frame(stringsAsFactors=F)
  tmp<-anno2[anno2$gene_id %in% tmp$gene_id,] %>% na.omit %>% reduce
  export(tmp,format = "BED",con=paste0("../ldsc/Resources/",der_names[i],".bed"))
}

files<-list.files("../ldsc",pattern = "results",full.names = T)

res<-data.frame(Cell_DMR=list.files("../ldsc",pattern = "results"),
                Enrichment=0,p=0)
res$Cat<-str_split(res$Cell_DMR,pattern = "_|\\.",simplify = T)[,2]
res$Cell_DMR<-str_split(res$Cell_DMR,pattern = "_|\\.",simplify = T)[,3]
for(i in 1:12){
  fi<-read.table(files[i],header = T)
  fi$Enrichment[2]->res[i,]$Enrichment
  fi$Enrichment_p[2]->res[i,]$p
}
ggplot(res,aes(x=Enrichment,y=-log(p)))+
  geom_point(aes(fill=Cell_DMR,colour=Cat,size=-log(p)),shape=21)+theme_bw()
