file=$1

seqtk sample -s100 ../Conversion/fqs/${file}_1.fq.gz 0.1 | gzip -1 > ../Conversion/fqs/${file}_1_sub.fq.gz
seqtk sample -s100 ../Conversion/fqs/${file}_2.fq.gz 0.1 | gzip -1 > ../Conversion/fqs/${file}_2_sub.fq.gz
bitmapperBS --search ./T2T_CHM13v2.0/hs1.fa --seq1 ../Conversion/fqs/${file}_1_sub.fq.gz --seq2 ../Conversion/fqs/${file}_2_sub.fq.gz -t 6 --bam -o ${file}.bam
rm ../Conversion/fqs/${file}_1*
rm ../Conversion/fqs/${file}_2*

sambamba_1.0.1 sort --tmpdir=./tmp -o ${file}.sorted.bam -t 4 ./${file}.bam
rm ${file}.bam*

metheor tag -i ${file}.sorted.bam -o ${file}_XM.bam -g ./T2T_CHM13v2.0/hs1.fa
metheor mhl -i ${file}_XM.bam -o ./MHL/${file}.mhl -d 5 -p 3
rm ${file}.sorted.bam*
rm ${file}_XM.bam*
