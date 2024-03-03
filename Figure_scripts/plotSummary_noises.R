#####The template scripts for plotting boxplot & dotplots #####
## We used ggplot2 and a ggplot-derived package ggpubr.

##@@ MPR complexity @@##

### Boxplot
PE<-read.table("MPR_meth_rb_with_info.txt")
PE$Group<-"LLI"
PE[PE$Age<=70,]$Group<-"Younger control"
PE[PE$Age>70 & PE$Age<90,]$Group<-"Elder control"
ggboxplot(PE,x="Group",y="MPR",color="Group",add = "jitter")+
  stat_compare_means(comparisons=list(c("LLI","Elder control"),
                                      c("Elder control","Younger control")))+
  theme_classic()

### dotplot
cor.test(PE[PE$group == "F1SP",]$ndemb3,PE[PE$group == "F1SP",]$Age)
ggplot(PE,aes(x=Age,y=MPR,color=group))+
  geom_point()+geom_smooth(method = "lm")+
  scale_color_manual(values = c("C"="red","F1SP"="grey4"))+
  theme_classic()

