###Install dependencies
conda install bwa
conda install -c bioconda fastp
conda install -c bioconda fastqc
conda install -c bioconda bedtools
conda install -c bioconda StringTie
conda install -c bioconda python2.7
conda install -c bioconda featureCounts
conda install -c bioconda picard-tools-2.2.1
conda install -c bioconda deeptools
conda install -c bioconda mosdepth
conda install -c bioconda seqtk
conda install -c bioconda macs2
conda install -c bioconda bowtie2
conda install -c bioconda seqkit
conda install -c bioconda subread
conda install -c bioconda ViennaRNA
conda install -c bioconda SeqPrep
conda install -c bioconda Bismark
conda install -c bioconda trim-galore
conda install -c bioconda trimmomatic
conda install -c bioconda phantompeakqualtools
conda install -c bioconda idr
conda install -c bioconda RNAfold
conda install -c bioconda STAR
conda install -c bioconda gfa

#bedToBigBed
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
chmod +x bedToBigBed
mv bedToBigBed /miniconda/bin
mv bedToBigBed ./miniconda3/envs/py2/bin

###Install R package
install.packages("ggplot2")
install.packages("viridis")
install.packages("viridisLite")
install.packages("data.table")
install.packages("dplyr")
install.packages("ggpubr")
install.packages("ggsignif")
install.packages("tidyverse")
install.packages("ggprism")
install.packages("DESeq2")
install.packages("patchwork")
install.packages("readxl")
install.packages("pheatmap")
install.packages("ggupset")

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("nucleR")
done

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("GenomicFeatures", "ChIPseeker"))
done
