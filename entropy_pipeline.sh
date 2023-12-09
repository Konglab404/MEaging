ID=$1

mkdir ~/data/Entropy/DMEAS/sample${ID}
cd ~/data/Entropy/DMEAS/sample${ID}
python2 ~/data/Entropy/DMEAS/test.spli.py ~/data/Entropy/CpG_rep/CpG_context_${ID}_1.trimmed_bismark_bt2_pe.txt.gz
mkdir ~/data/Entropy/DMEAS/sample${ID}/bisdir
mv ~/data/Entropy/DMEAS/sample${ID}/*.bis ~/data/Entropy/DMEAS/sample${ID}/bisdir
cd ~/data/Entropy/DMEAS/

perl ~/data/Entropy/DMEAS/DMEAS_3_cpg.pl -d gw -s ms ~/data/Entropy/DMEAS/sample${ID}/bisdir ~/data/Entropy/DMEAS/hg19gp -o ~/data/Entropy/DMEAS/sample${ID}/
rm -r ~/data/Entropy/DMEAS/sample${ID}/bisdir




````
library(stringr)
files<-list.files(pattern = "3CG.txt",recursive=T)
res<-data.frame()
for(i in 1:22){
  tmp<-read.table(files[i],skip=1)
  tmp$Chr<-paste0("chr",str_split(files[i],pattern="\\/",simplify=T)[1])
  colnames(tmp)<-c("segment","reads","Level","Entropy","MethyCG","totalCG","Chr")
  res<-rbind(res,tmp)
}
write.table(res,file=paste0("Entropy",str_split(getwd(),pattern="sample",simplify=T)[2]))
`````
