###MAV
setwd("/data1/wanght/Single_cell")
library(Matrix)
library(Seurat)
library(dplyr)
library(Seurat)
#library(useful)
library(cowplot)
library(harmony)
library(annotatr)

PBMC_SC <- readRDS(file = "./PNAS19/PNAS2019_processed_seurat.rds") #Japan
PBMC_45 <- readRDS(file = "./NG2018/NG2018_Seurat_processed.RDS") #NG 

#Japan PNAS data
dir.create("Japan_PBMC_SC")
IDs<-unique(PBMC_SC$ID)

resMAV<-data.frame(gene="t")

for(id in IDs){
  expr<-PBMC_SC[,PBMC_SC$ID == id]@assays$RNA@data
  #expr <- log2(expr + 1)
  expr <- expr[rowMeans(expr > 0.01) > 0.01,]
  m <- rowMeans(expr)
  lm <- log(m)
  var <- apply(expr,1,var)
  logvar <- log(var)
  hypervar_var <- resid(loess(var~lm))
  hypervar_logvar <- resid(loess(logvar~lm))
  res <- data.frame(mean=m,var=var,hypervar_var=hypervar_var,hypervar_logvar=hypervar_logvar)
  resT<-data.frame(gene=rownames(res),MAV=res$hypervar_var)
  colnames(resT)[2]<-id
  resMAV<-merge(resMAV,resT,by='gene',all=T)
  
  #saveRDS(res,file=paste0('./Japan_PBMC_SC/MAV_',id,"_",unique(PBMC_SC@meta.data$Age[PBMC_SC@meta.data$ID == id]),"_yrs.RDS"))
}

resMAV_quant<-resMAV
for(i in 2:ncol(resMAV_quant)){
  quantile(na.omit(resMAV[,i]),prob=seq(0,1,0.01))->x
  findInterval(resMAV[,i],x)->resMAV_quant[,i]
}
save(resMAV_quant,resMAV,file="Japan_PBMC_SC/MAV_all_res.RData")

####
load("Japan_PBMC_SC/MAV_all_res.RData")
resMAV_quant<-na.omit(resMAV_quant)
setwd("~/data/Entropy/Result/")
load("DER_CEN_20230701.RData")
annos<-c('hg19_basicgenes')
annotations = build_annotations(genome = 'hg19', annotations = annos)
CEN_deranno<-annotate_regions(regions = DER_CEN,
                              annotations = annotations,ignore.strand = TRUE,quiet = FALSE)
CEN_deranno<-cbind(CEN_deranno$loc,as.data.frame(CEN_deranno$annot))
CEN_deranno<-unique(CEN_deranno[,c("seqnames","start","type","gene_id")])

CEN_deranno[grep(CEN_deranno$type,pattern = "promo"),]$gene_id %>% 
  table %>% as.data.frame(stringsAsFactors=F)->promoter_table
promoter_table[promoter_table$Freq>=1,]$. %>% unique -> Gene_symbol
ensembls <- mapIds(org.Hs.eg.db, keys = Gene_symbol, keytype = "ENTREZID", column="ENSEMBL")

low_der_MAV<-resMAV[resMAV$gene %in% ensembls,]
cens<-low_der_MAV[,2:8] %>% apply(1,mean,na.rm=T)
controls<-low_der_MAV[,9:13] %>% apply(1,mean,na.rm=T)
wilcox.test(cens,controls) #p=4.888e-15

mean(cens,na.rm=T)
mean(controls,na.rm=T)
std.error(cens)
std.error(controls)
ggplot(data=data.frame(type=c("CEN","Controls"),
                  MAV=c(mean(cens,na.rm=T),mean(controls,na.rm=T)),
                  SE=c(std.error(cens),std.error(controls))),aes(x=type,y=MAV))+
  geom_point(stat = "identity")+
  geom_errorbar(aes(ymax=MAV+SE,ymin=MAV-SE))+coord_flip()+theme_bw()

mean(cens,na.rm=T)/mean(controls,na.rm=T) #0.5672661

other_der_MAV<-resMAV[!resMAV$gene %in% ensembls,]

other_cens<-other_der_MAV[,2:8] %>% apply(1,mean,na.rm=T)
other_controls<-other_der_MAV[,9:13] %>% apply(1,mean,na.rm=T)
ggplot(data=data.frame(type=c("CEN","Controls"),
                       MAV=c(mean(other_cens,na.rm=T),mean(other_controls,na.rm=T)),
                       SE=c(std.error(other_cens),std.error(other_controls))),aes(x=type,y=MAV))+
  geom_point(stat = "identity")+
  geom_errorbar(aes(ymax=MAV+SE,ymin=MAV-SE))+coord_flip()+theme_bw()
mean(other_cens,na.rm=T)/mean(other_controls,na.rm=T) #0.9448796
wilcox.test(other_cens,other_controls)

low_der_MAV<-resMAV[!resMAV$gene %in% ensembls,]
control<-low_der_MAV[,2:8] %>% apply(1,mean,na.rm=T)
cens<-low_der_MAV[,9:13] %>% apply(1,mean,na.rm=T)
t.test(control,cens)

set.seed(19940731)
res_pert<-rep(0,1000)
for(i in 1:1000){
  tmp<-resMAV[sample(1:nrow(resMAV),size = length(Gene_symbol)),]
  tmp_cens<-tmp[,2:8] %>% apply(1,mean,na.rm=T)
  tmp_controls<-tmp[,9:13] %>% apply(1,mean,na.rm=T)
  res_pert[i]<-mean(tmp_cens,na.rm=T)/mean(tmp_controls,na.rm=T)
}
hist(res_pert,breaks = 50,xlim = c(0.5,1))
sum(res_pert< 0.5672661)/1000
abline(v = 0.5672)

sampleinfo<-read.table("PNAS19/Sampleinfo.txt",header = T)








##############################################
#NG data
setwd("~/data/Entropy/Result/")
dir.create("NG_PBMC_SC")
IDs<-unique(PBMC_45$Donor_ID)
IDs<-IDs[-grep(IDs,pattern = "LL")]

resMAV<-data.frame(gene="t")

for(id in IDs){
  expr<-PBMC_45[,PBMC_45$Donor_ID == id]@assays$RNA@data
  #expr <- log2(expr + 1)
  expr <- expr[rowMeans(expr > 0.1) > 0.01,]
  m <- rowMeans(expr)
  lm <- log(m)
  var <- apply(expr,1,var)
  logvar <- log(var)
  hypervar_var <- resid(loess(var~lm,))
  hypervar_logvar <- resid(loess(logvar~lm))
  res <- data.frame(mean=m,var=var,hypervar_var=hypervar_var,hypervar_logvar=hypervar_logvar)
  resT<-data.frame(gene=rownames(res),MAV=res$hypervar_var)
  colnames(resT)[2]<-id
  resMAV<-merge(resMAV,resT,by='gene',all=T)
  
  #saveRDS(res,file=paste0('./Japan_PBMC_SC/MAV_',id,"_",unique(PBMC_SC@meta.data$Age[PBMC_SC@meta.data$ID == id]),"_yrs.RDS"))
}

resMAV_quant<-resMAV
for(i in 2:ncol(resMAV_quant)){
  quantile(na.omit(resMAV[,i]),prob=seq(0,1,0.1))->x
  findInterval(resMAV[,i],x)->resMAV_quant[,i]
}
save(resMAV_quant,resMAV,file="./NG_PBMC_SC/MAV_all_res.RData")

####
unique(PBMC_45@meta.data[,c("Donor_ID","Age")])->sampleinfo
rownames(sampleinfo)<-sampleinfo$Donor_ID
sampleinfo<-sampleinfo[IDs,]

cell_num<-as.data.frame(table(PBMC_45$Donor_ID))

CEN_deranno[grep(CEN_deranno$type,pattern = "promo"),]$gene_id %>% 
  table %>% as.data.frame(stringsAsFactors=F)->promoter_table
promoter_table[promoter_table$Freq>=3,]$. %>% unique -> Gene_symbol
ensembls <- mapIds(org.Hs.eg.db, keys = Gene_symbol, keytype = "ENTREZID", column="ENSEMBL")

Aging_der_MAV<-resMAV[resMAV$gene %in% ensembls,]
meanAging_MAV<-apply(Aging_der_MAV[,2:46], 2,mean,na.rm=T) %>% as.data.frame()
MAV_Aging<-merge(sampleinfo,meanAging_MAV,by.x="Donor_ID",by.y="row.names")

MAV_Aging<-merge(cell_num,MAV_Aging,by.x="Var1",by.y="Donor_ID")

colnames(MAV_Aging)[4]<-"value"
quantile(MAV_Aging$Freq,seq(0,1,0.1)) #10%:337.6 %90:760.8
MAV_A<-MAV_Aging[(MAV_Aging$Age %in% (29:69)) & (MAV_Aging$Freq>340 & MAV_Aging$Freq<760), ]
plot(MAV_A$Age,MAV_A$value)
cor.test(MAV_A$Age,MAV_A$value,method = "kendall",continuity = T,exact = F,alternative = "greater")
aov(value~Age,data=MAV_A) %>% summary
ggplot(MAV_A,aes(x=as.factor(Age),group=as.factor(Age),y=value))+geom_violin(width=1.2)+geom_boxplot(width=0.1)+theme_bw()
ggplot(MAV_A,aes(x=Age,y=value))+geom_point()+geom_smooth(method = "lm")+theme_bw()

########################## CellDMR ###########
setwd("~/data/Entropy/Result/")



library(SingleR)
library(celldex)
ref.data<-HumanPrimaryCellAtlasData(ensembl = T)
pred<-SingleR(test = as.matrix(PBMC_SC@assays$RNA@data),
              ref=ref.data,
              labels = ref.data$label.main,
              method = "cell")
