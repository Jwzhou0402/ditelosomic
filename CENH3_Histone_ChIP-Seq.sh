
#######ChIP-seq
##CENH3
for i in `cat ./CENH3_CSDt.txt`
do

bsub -J bwa -n 10 -R span[hosts=1] -o %J.out -e %J.err -q normal "fastp -p -w 15 -l 30 -i ~/CHIP_seq/CENH3/RawData/"$i"_R1.fq.gz -I ~/CHIP_seq/CENH3/RawData/"$i"_R2.fq.gz -o ~/CHIP_seq/CENH3/00_qc/"$i"_clean.R1.fastq.gz -O ~/CHIP_seq/CENH3/00_qc/"$i"_clean.R2.fastq.gz -h ~/CHIP_seq/CENH3/00_qc/"$i".html;\

bwa mem -t 20 ~/CSGL/wheat.CSGL.fa ~/CHIP_seq/CENH3/00_qc/"$i"_clean.R1.fastq.gz ~/CHIP_seq/CENH3/00_qc/"$i"_clean.R2.fastq.gz | awk '{if(\$0~/^@/||\$5>=20) {print \$0}}' | samtools sort -@ 10 -o ~/CHIP_seq/CENH3/02_mapping/"$i"_mapq20_sort.bam;\

samtools rmdup ~/CHIP_seq/CENH3/02_mapping/"$i"_mapq20_sort.bam ~/CHIP_seq/CENH3/03_Rmdup/"$i"_mapq20_sort_rmdup.bam;\

samtools flagstat ~/CHIP_seq/CENH3/03_Rmdup/"$i"_mapq20_sort_rmdup.bam > ~/CHIP_seq/CENH3/05_flagstat/"$i"_flagstat.txt;\

java -jar /public/home/jingwzhou/miniconda3/share/picard-2.26.6-0/picard.jar CollectInsertSizeMetrics I= ~/CHIP_seq/CENH3/02_mapping/"$i"_mapq20_sort.bam O= ~/CHIP_seq/CENH3/04_CollectInsertSizeMetrics/len_distri_"$i".txt H= ~/CHIP_seq/CENH3/04_CollectInsertSizeMetrics/04_CollectInsertSizeMetrics/len_distri_"$i".pdf

bamToBed -i ~/CHIP_seq/CENH3/03_Rmdup/"$i"_mapq20_sort_rmdup.bam | sort -k 1,1 -k 2,2n > ~/CHIP_seq/CENH3/06_bamToBed/"$i"_mapq20_sort_rmdup.bed

bamCoverage -b ~/CHIP_seq/CENH3/03_Rmdup/"$i"_mapq20_sort_rmdup.bam --binSize 10 --normalizeUsing RPKM -p 10 -o ~/CHIP_seq/CENH3/03_Rmdup/"$i"_bin10_RPKM.bw;\

bamCompare -b1 ~/CHIP_seq/CENH3/03_Rmdup/"$i"_mapq20_sort_rmdup.bam -b2 ~/DNA_seq/03_rmdu/final_bam/euploid/"$i"_sort_dedupilcate.bam --scaleFactorsMethod None --normalizeUsing RPKM --operation log2 --outFileName ./"$i".CENH3_VS_"$i".reDNA.RPKM.log2.bin10.bw --binSize 10 -p 15"

sleep 10
done

####ReadsCounts
for sample in ./*_mapq20_sort_rmdup.bed; do
index=$(basename $sample |sed 's/_mapq20_sort_rmdup.bed//')
prefix=$(dirname $sample)

bsub  -J samtools -n 1 -o 1kb_${index}-%J.out -e 1kb_${index}-%J.err -R span[hosts=1] "perl ./bed2-1Kb_wheat_CSGL.pl ${prefix}/${index}_mapq20_sort_rmdup.bed ./${index}_CENH3_1kb_ReadsCounts.txt;\
perl ./bed2-10Kb_wheat_CSGL.pl ${prefix}/${index}_mapq20_sort_rmdup.bed ./${index}_CENH3_10kb_ReadsCounts.txt;\
perl ./bed2-100Kb_wheat_CSGL.pl ${prefix}/${index}_mapq20_sort_rmdup.bed ./${index}_CENH3_100kb_ReadsCounts.txt"

sleep 10
done


#-------------------------------------Dt7DS-CENH3-Chip---Dt7DS_chr7D_CENH3_1kb_ReadsCounts.txt
# Dt7DS_chr7D_CENH3_1kb_ReadsCounts.txt
###merge###
CS_Dt7DS_Dt7DL_1kb_chr7D_CENH3 <- read.table("./merge_CSX5_Dt7DS_Dt7DL_1kb_chr7D.txt", head=T,fill = TRUE, sep = "\t")
x_min1 <- min(CS_Dt7DS_Dt7DL_1kb_chr7D_CENH3[1536:3734,"Start"])
x_max1 <- max(CS_Dt7DS_Dt7DL_1kb_chr7D_CENH3[1536:3734,"Start"])
y_min1 <- min(CS_Dt7DS_Dt7DL_1kb_chr7D_CENH3[1536:3734,"CS"])
y_max1 <- max(CS_Dt7DS_Dt7DL_1kb_chr7D_CENH3[1536:3734,"CS"])

x_min2 <- min(CS_Dt7DS_Dt7DL_1kb_chr7D_CENH3[4464:6950,"Start"])
x_max2 <- max(CS_Dt7DS_Dt7DL_1kb_chr7D_CENH3[4464:6950,"Start"])
y_min2 <- min(CS_Dt7DS_Dt7DL_1kb_chr7D_CENH3[4464:6950,"CS"])
y_max2 <- max(CS_Dt7DS_Dt7DL_1kb_chr7D_CENH3[4464:6950,"CS"])

pdf('CS_Dt7DS_Dt7DL_1kb_chr7D_CENH3.pdf', width=4, height=3)
ggplot(CS_Dt7DS_Dt7DL_1kb_chr7D_CENH3, aes(x = Start)) +
  geom_line(aes(y = CS, color = "CS")) +
  geom_line(aes(y = Dt7DS, color = "Dt7DS")) +
  geom_line(aes(y = Dt7DL, color = "Dt7DL")) +
  labs(title = "CS and CSDt7DS read Counts of chr7D within 1kb window",
       x = "Position", y = "Read Counts") +
  annotate("rect", xmin = x_min1, xmax = x_max1,
           ymin = y_min1, ymax = y_max1,
           alpha = 0.2)+
  annotate("rect", xmin = x_min2, xmax = x_max2,
           ymin = y_min2, ymax = y_max2,
           alpha = 0.2)+
  theme_classic()+
  scale_color_manual(values = c("CS" = "#f0a699", "Dt7DS" = "#e9a96b", "Dt7DL" = "#219ebc"))
dev.off()
print(paste0())

#r_squared
r <- read.table("./R2_merge_CS_Dt7DS_1kb_chr7D_CENH3.txt", head=T,fill = TRUE, sep = "\t")
r <- cor(r$CS, r$Dt7DS)
r_squared <- r^2
print(paste("R-squared: ", round(r_squared, 4)))

CS_Dt7DL_1kb_r <- read.table("./R2_merge_CS_Dt7DL_1kb_chr7D_CENH3.txt", head=T,fill = TRUE, sep = "\t")
r <- cor(CS_Dt7DL_1kb_r$CS, CS_Dt7DL_1kb_r$Dt7DL)
r_squared <- r^2
print(paste("R-squared: ", round(r_squared, 4)))


##Correlation between chromosome length and centromere size and linear regression equation
library(ggplot2)
library(viridis)
library(viridisLite)
library(data.table)
library(dplyr)
library(ggpubr)

data <- read.table("./CentromereSize_ChromosomeLength.txt",header = TRUE)
pdf('CentromereSize_ChromosomeLength.pdf', width=4, height=3)
ggplot(data, aes( x = Centromere, y = Chromosome, color = group)) + theme_classic()+ geom_point( aes( color = group), size=1)+ geom_smooth(method = 'lm', formula = y ~ x, se = T) + stat_cor(data=data, method = "pearson")
dev.off()






library(ggplot2)
library(ggpubr)
library(ggsignif)
library(tidyverse)
library(ggprism)

#TelocentricOther
df1 <- read.delim("./TelocentricOther.txt")
pdf('Fig4E_TelocentricOther.pdf', width=3.2, height=2.3)
p1 <-ggplot(df1,aes(x=group,y=values))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+
  geom_boxplot(aes(fill=group),
               outlier.colour="white",size=0.8)+
  theme(panel.background =element_blank(),
        axis.line=element_line(),
        legend.position="none",plot.title = element_text(size=14))+
  geom_jitter(width = 0.2)+
  ggtitle("boxplot")+
  geom_signif(comparisons = list(c("Telocentric","Other")),
                     map_signif_level = T,
                     test = t.test,
                     y_position = c(1),
                     tip_length = c(c(0.05,0.05)),
                     size=0.8,color="black")
p1
dev.off()

#Fig4E_SVregion_WithoutSVregion
df2 <- read.delim("./SVregion_WithoutSVregion.txt")
pdf('SVregion_WithoutSVregion.pdf', width=4, height=3)
p2 <-ggplot(df2,aes(x=group,y=Pearson))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+
  geom_boxplot(aes(fill=group), 
               outlier.colour="white",size=0.8)+
  theme(panel.background =element_blank(),
        axis.line=element_line(),
        legend.position="none",plot.title = element_text(size=14))+
  geom_jitter(width = 0.2)+
  ggtitle("boxplot")+ 
  geom_signif(comparisons = list(c("SV region","Without SV region")),
                     map_signif_level = T,
                     test = t.test,
                     y_position = c(1),
                     tip_length = c(c(0.05,0.05)),
                     size=0.8,color="black")
p2
dev.off()





#########Histone ChIP seq
###rep1_rep2
for i in `cat ./Histone_ChIP.txt`
do
bsub -J Histone -n 15 -o ../report/%J."$i".out -e ../report/%J."$i".err -R span[hosts=1] -q smp "fastp -p -w 15 -l 30 -i ./"$i"_1.fastq.gz -I ./"$i"_2.fastq.gz -o ./"$i"_1_clean.fq.gz -O ./"$i"_2_clean.fq.gz -h ./"$i".html;\

bowtie2 --no-unal --threads 20 --sensitive -k 3 -q --phred33 -x ~/CSGL/bowtie2-build/wheat.CSGL -1 ./"$i"_1_clean.fq.gz  -2 ./"$i"_2_clean.fq.gz | awk '{if(\$0~/^@/||\$5>=20) {print \$0}}' | samtools sort -@ 15 -o ./"$i"_q20_sort.bam;\

samtools index -c -@ 10 ./"$i"_q20_sort.bam;\

java -jar picard.jar MarkDuplicates \
      I= ./"$i"_q20_sort.bam \
      O= ./"$i".final.bam \
   REMOVE_DUPLICATES=true \
      M= ./"$i"_sort_metrics.txt;\

samtools index -c -@ 10 ./"$i".final.bam;\

samtools flagstat ~/CHIP_seq/H3K4/rmdup/"$i".final.bam > ~/CHIP_seq/H3K4/flagstat/"$i".flagstat;\

bamToBed -i ./"$i".final.bam | sort -k 1,1 -k 2,2n > ./"$i".final.bed;\

bamCoverage -b ./"$i".final.bam --binSize 150 --normalizeUsing RPKM -p 8 -o "$i"_RPKM_150bin.bw"

sleep 10
done

###rep0_100kb_ReadCounts
for sample in /*_rep1.bam; do
index=$(basename $sample |sed 's/_rep1.bam//')
prefix=$(dirname $sample)

bsub  -J perl -n 1 -o ReadsCounts-${index}-%J.out -e ReadsCounts-${index}-%J.err -R span[hosts=1] "samtools merge ${index}_rep0.bam ${index}_rep1.bam ${index}_rep2.bam;\
samtools index -c ${index}_rep0.bam;\

bamCoverage -b ./${index}_rep0.bam --binSize 150 --normalizeUsing RPKM -p 8 -o ${index}_rep0_RPKM_150bin.bw

bamToBed -i ./${index}_rep0.bam | sort -k 1,1 -k 2,2n > ./${index}_rep0.bed;\

perl ./bed2-100Kb_wheat_CSGL.pl ${prefix}/${index}_rep0.bed ./${index}_ChIP_100kb_ReadsCounts.txt;\

perl ./bed2-10Kb_wheat_CSGL.pl ${prefix}/${index}_rep0.bed ./${index}_ChIP_10kb_ReadsCounts.txt;\

perl ./bed2-1Kb_wheat_CSGL.pl ${prefix}/${index}_rep0.bed ./${index}_ChIP_1kb_ReadsCounts.txt"
done
