setwd("~/data/Entropy/")
dir.create("Result")
#library(bsseq)
library(stringr)
library(magrittr)
library(data.table)

#library(BSgenome.Hsapiens.UCSC.hg38)

sample133<-read.table("Sampleinfo133",header = T)
sample133_files<-paste0("Ent89/Entropy",sample133$ID)
sample133_files %in% list.files("Ent89",full.names = T) %>% sum #all 133
Ent_files <- sample133_files

Entropy1<-fread(Ent_files[1])
Entropy1<-Entropy1[!startsWith(Entropy1$ID,"_"),]
Entropy1$loc = paste0(Entropy1$seqnames,":",Entropy1$ID)
Entropy1<-Entropy1[!startsWith(Entropy1$ID,"_"),]
ent_aggre11<-data.table(loc=Entropy1$loc,Entropy340=(Entropy1$ent))

for(i in 2:133){
  tmp<-fread(Ent_files[i])
  if(ncol(tmp) == 8){
    tmp<-tmp[!startsWith(tmp$segment,"_"),2:8]
    tmp$loc = paste0(tmp$Chr,":",tmp$segment)
    tmp<-data.table(loc=tmp$loc,tmp=(tmp$Entropy))
  }else{
    tmp<-tmp[!startsWith(tmp$ID,"_"),]
    tmp$loc = paste0(tmp$seqnames,":",tmp$ID)
    tmp<-data.table(loc=tmp$loc,tmp=(tmp$ent))
  }

  colnames(tmp)[2]<-str_split(Ent_files[i],pattern = "/",simplify = T)[,2]
  
  ent_aggre11<-merge.data.table(ent_aggre11,tmp,by="loc",all.x=T,all.y = T)
  cat(i);cat(" ")
  gc()
}
saveRDS(ent_aggre11,file="Result/Entropy_133samples.aggregate.RDS")

#na_nums<-apply(ent_aggre11[,2:133],1,FUN = function(x){return(sum(!is.na(x)))})
#sum(apply(test1[,2:177],1,max,na.rm=T)>0)
# [1] 115

#########Only samples in WGBS manu#####
ent133<-readRDS("./Result/Entropy_133samples.aggregate.RDS")
ent133<-ent_aggre11
rm(ent_aggre11);gc()

####### 176 samples #########
sample133<-read.table("Sampleinfo133",header = T)
sample133$ID<-paste0("Entropy",sample133$ID)
IDs<-sample133$ID
CEN<-sample133[sample133$Age>=95,]$ID
YC<-sample133[sample133$Age <=70,]$ID
EC<-setdiff(IDs,c(CEN,YC))

entCEN<-ent133[,..CEN]
na_CEN<-apply(entCEN,1,FUN = function(x){return(sum(!is.na(x)))})
rm(entCEN);gc()

entYC<-ent133[,..YC]
na_YC<-apply(entYC,1,FUN = function(x){return(sum(!is.na(x)))})
rm(entYC);gc()

entEC<-ent133[,..EC]
na_EC<-apply(entEC,1,FUN = function(x){return(sum(!is.na(x)))})
rm(entEC);gc()

sum(na_CEN>0.5*length(CEN) & na_EC>0.5*length(EC) & na_YC>0.5*length(YC)) #6087544 EC 14 YC 30
na_YC[na_CEN>0.5*length(CEN) & na_EC>0.5*length(EC) & na_YC>0.5*length(YC)] %>% summary 
na_EC[na_CEN>0.5*length(CEN) & na_EC>0.5*length(EC) & na_YC>0.5*length(YC)] %>% summary 


sum((na_CEN+na_EC+na_YC)>0.6*133) #EC median 10 YC median 30 14084997
na_YC[(na_CEN+na_EC+na_YC)>0.6*133] %>% summary 
na_EC[(na_CEN+na_EC+na_YC)>0.6*133] %>% summary 

sum((na_CEN > 0.5*length(CEN)) & (na_EC+na_YC)>27) #EC median 9 YC median 29 15770197
na_YC[(na_CEN > 0.5*length(CEN)) & (na_EC+na_YC)>27] %>% summary 
na_EC[(na_CEN > 0.5*length(CEN)) & (na_EC+na_YC)>27] %>% summary 


#############133 samples meth##########
sample133<-read.table("Sampleinfo133",header = T)
sample133_files<-paste0("Ent89/Entropy",sample133$ID)
sample133_files %in% list.files("Ent89",full.names = T) %>% sum #all 133
Ent_files <- sample133_files

Meth1<-fread(Ent_files[1])
Meth1<-Meth1[!startsWith(Meth1$segment,"_"),2:8]
Meth1$loc = paste0(Meth1$Chr,":",Meth1$segment)
ent_aggre11<-data.table(loc=Meth1$loc,Meth340=(Meth1$Level))

for(i in 2:133){
  tmp<-fread(Ent_files[i])
  if(ncol(tmp) == 8){
    tmp<-tmp[!startsWith(tmp$segment,"_"),2:8]
    tmp$loc = paste0(tmp$Chr,":",tmp$segment)
    tmp<-data.table(loc=tmp$loc,tmp=(tmp$Level))
  }else{
    tmp<-tmp[!startsWith(tmp$ID,"_"),]
    tmp$loc = paste0(tmp$seqnames,":",tmp$ID)
    tmp<-data.table(loc=tmp$loc,tmp=(tmp$meth))
  }
  
  colnames(tmp)[2]<-str_split(Ent_files[i],pattern = "/",simplify = T)[,2]
  
  ent_aggre11<-merge.data.table(ent_aggre11,tmp,by="loc",all.x=T,all.y = T)
  cat(i);cat(" ")
  gc()
}
saveRDS(ent_aggre11,file="Result/Meth_133samples.aggregate.RDS")

