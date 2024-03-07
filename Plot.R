setwd("/home/owht/Data/KIZ/data/Entropy2_rev/Region_entropy")
sampleinfo<-read.table("../PE/fis_meth_rb_with_info.txt")
sampleinfo<-sampleinfo[,c(6:9)]  
load("TE_enh_insu.RData")
rownames(sampleinfo) == Results$ID
sampleinfo<-cbind(sampleinfo,Results)
sampleinfo$Group<-"LLI"
sampleinfo[sampleinfo$Age<=70,]$Group<-"Younger control"
sampleinfo[sampleinfo$Age>70 & 
             sampleinfo$Age<90,]$Group<-"Elder control"
sampleinfo$Group<-factor(sampleinfo$Group,levels = c("Younger control",
                                                     "Elder control",
                                                     "LLI"))

library(reshape2)
plots<-melt(sampleinfo[,c(6:12)],measure.vars = 1:6)
library(ggpubr)
ggboxplot(plots, x="Group", y="value", 
          color="Group", palette="npg",facet.by = "variable")+
  stat_compare_means(comparisons = 
                       list(c("Younger control","Elder control"),
                            c("Elder control","LLI"),
                            c("LLI","Younger control")),
                     method = "wilcox.test")
