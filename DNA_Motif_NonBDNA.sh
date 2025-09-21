#####01_bwa index
bsub -J fastp -n 8 -o 01_index-%J.out -e 01_index-%J.err -R span[hosts=1] "bwa index -a bwtsw wheat.CSGL.fa"

#####02_qc
fastqc -o ~/CHIP_seq/H3K4/QC_after/trimmomatic/ ~/CHIP_seq/H3K4/QC_after/trimmomatic/*paired.fq.gz;\

for sample in ./*_1.fq.gz; do
index=$(basename $sample |sed 's/_1.fq.gz//')
prefix=$(dirname $sample)
bsub -J fastp -n 8 -o 01_qc-%J.out -e 01_qc-%J.err -R span[hosts=1] "fastp -p -w 15 -l 30 -i ${prefix}/${index}_1.fq.gz -I ${prefix}/${index}_2.fq.gz -o ${prefix}/${index}_clean.R1.fastq.gz -O ${prefix}/${index}_clean.R2.fastq.gz -h ${prefix}/${index}.html"
done

#####03_mapping_dedupilcate_normalized
for sample in ./*_1_clean.fq.gz; do

index=$(basename $sample |sed 's/_1_clean.fq.gz//')
prefix=$(dirname $sample)

bsub -J ${index} -n 8 -o %J.${index}.out -e %J.${index}.err -R span[hosts=1] "bwa mem -t 10 ./bwa_index/wheat.CSGL.fa ${prefix}/${index}_1_clean.fq.gz ${prefix}/${index}_2_clean.fq.gz | awk '{if(\$0~/^@/||\$5>=20) {print \$0}}' | samtools sort -@ 20  -o ${index}_sort.bam;\

java -jar /public/home/jingwzhou/miniconda3/share/picard-2.26.6-0/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${prefix}/${index}_sort.bam O=${prefix}/${index}_sort_dedupilcate.bam METRICS_FILE=${prefix}/${index}_sort_metrics.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=tmp${prefix}/${index} ;\

samtools index -c -@ 10 ${prefix}/${index}_sort_dedupilcate.bam;\

samtools flagstat ${prefix}/${index}_sort_dedupilcate.bam > ${index}_sort_dedupilcate.flagstat;\

bamToBed -i ${prefix}/${index}_sort_dedupilcate.bam | sort -k 1,1 -k 2,2n > ${prefix}/${index}_sort_dedupilcate.bed;\

bamCoverage -of bigwig -b ${prefix}/${index}_dedupilcate.bam --binSize 50 --normalizeUsing RPGC --effectiveGenomeSize 14607800139 -p 10 -o ${index}_1Xnormalized_50.bw;\

bamCoverage -b ${prefix}/${index}_dedupilcate.bam --binSize 10 --normalizeUsing RPKM -p 10 -o ${index}_RPKMnormalized_10.bw"
sleep 10
done


#####04_CNV
###04.1_cnvnator
source /public/home/software/.bashrc
module load CNVnator

for sample in ./*_sort_dedupilcate.bam; do

index=$(basename $sample |sed 's/_sort_dedupilcate.bam//')
prefix=$(dirname $sample)

bsub -J cnvnator -n 10 -o cnvnator-${index}_2000-%J.out -e cnvnator-${index}_2000-%J.err -R span[hosts=1] -q normal "cnvnator -genome ~/DNA_seq/12_CNV/split/wheat.CSGL.fa -root ${index}_2000.root -tree ~/DNA_seq/11_CALLSNP/06_rmdup/${index}_final.bam;\
cnvnator -genome ~/DNA_seq/12_CNV/split/wheat.CSGL.fa -root ${index}_2000.root -his 2000 -d ~/DNA_seq/12_CNV/split;\
cnvnator cnvnator -genome ~/DNA_seq/12_CNV/split/wheat.CSGL.fa -root ${index}_2000.root -stat 2000;\
cnvnator -root ${index}_2000.root -eval 2000 > ${index}_2000.eval.ratio;\
cnvnator -root ${index}_2000.root -partition 2000;\
cnvnator -root ${index}_2000.root  -call 2000 > ${index}_2000.cnv;\
cnvnator2VCF.pl ${index}_2000.cnv > ${index}_2000.cnv.vcf"
done

###04.2_mosdepth
for sample in ~/DNA_seq/03_rmdu/*_dedupilcate.bam; do

index=$(basename $sample |sed 's/_dedupilcate.bam//')
prefix=$(dirname $sample)

bsub  -J ${index} -n 8 -o %J.${index}.out -e %J.${index}.err -R span[hosts=1] "mosdepth -n -t 4 -F 1796 -T 1,2,3,4,5 --by 100000 ${prefix}/${index}_100Kb ${prefix}/${index}_dedupilcate.bam;\

less ${index}_100Kb.regions.bed.gz | cut -f 1,2,4 | awk 'BEGIN{print "Chr\tStart\tCounts"}1' > ${index}_100Kb.regions_cov.txt"
sleep 10
done

##04.3_plot
source /public/home/software/.bashrc
module load RStudio/4.2

setwd("D:") 
library(ggplot2)
library(viridis)
library(viridisLite)
library(data.table)
##############Dt1AL_100Kb_215-219.regions_cov.txt
data1 <- read.delim("./Dt1AL_100Kb.regions_cov.txt") 
data1_less2 <- subset(data1, Counts <=80) 
pdf(paste0("Dt1AL_100Kb.pdf"), width=1.2, height=1)
ggplot(data1_less2, aes(x = Start, y = Counts)) +
  geom_area(fill = "#3C92C6", alpha = 1) +
  scale_color_viridis("Depth", option = "viridis") +
  theme_classic() +
  theme(axis.line = element_line(size = 0.5))+
  ylab("Read Depth") +
  xlab("Position") +
  ggtitle("Dt1AL") +
  theme(plot.title = element_text(hjust = 0.5, size = 4, face = "bold"),
        axis.title.x = element_text(size = 3, face = "bold"),
        axis.title.y = element_text(size = 3, face = "bold"),
        axis.text = element_text(size = 4, face = "bold"),
        legend.key.height = unit(1,'cm')) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1))) 
dev.off()
print

########detection of chromosomal structural variations in wheat ditelosomic stocks
cols<-c('#6EAC8B','#455E9B')
pal<-colorRampPalette(cols)
image(x=1:20,y=1,z=as.matrix(1:20),col=pal(20))
Dt1BS_subset_less2 <- subset (Dt1BS_sort_100Kb.regions_cov, Counts <= 7) 
pdf(paste0("Dt1BS_100Kb_240108.pdf"), width=3, height=1.5)
ggplot(Dt1BS_subset_less2[CSDT_1BS_subset_less2$Chr=="chr7A",], aes(x = Start/1000000, y = Counts, color = Counts))+
  geom_point(size = 0.05)+
  scale_color_gradientn(colors = pal(20))+ 
  theme_classic()+
  ylab("Read Depth")+
  xlab("Position")+
  ggtitle("Dt1BS read depth of Chr7A within 100Kb window")+
  theme(plot.title = element_text(hjust = 0.5,size = 4, face = "bold"),
        axis.title.x = element_text(size = 5, face = "bold"),
        axis.title.y = element_text(size = 5, face = "bold"),
        axis.text = element_text(size = 5, face = "bold"),
        legend.key.height = unit(0.5,'cm'),
        legend.key.width = unit(0.2, 'cm') 
  ) +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE))
dev.off()
print(paste0())
