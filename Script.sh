#!/bin/bash
#Download SRA file, Cut adapter sequence CTGTAGGCACCATCAAT  and remove sequence shorter than 3
cd /Users/lanqingying/Downloads/sratoolkit.2.8.2-1-mac64/bin
./prefetch SRR
./fastq-dump SRR.sra
~/Library/Python/2.7/bin/cutadapt -a -m 3 XXXXX SRR.fastq>SRR_trimmed.fastq

#remove trRNA
export PATH=$PATH:~/Downloads/bowtie2-2.2.9 
cd /Users/lanqingying/Documents/SRA/Bowtie2
bowtie2-build -f file.fa trRNA
bowtie2 -x trRNA -U SRR_trimmed.fastq -S output.sam —-un SRR_removed.fastq —-quiet
rm output.sam

#Length distribution
python LengthHist.py

#map to reference genome
export PATH=$PATH:~/Downloads/bowtie2-2.2.9 
cd /Users/lanqingying/Documents/SRA/Bowtie2
bowtie2 -x Ecoli -U SRR_removed.fastq -S SRR.sam  —-quiet

export PATH=/usr/local/samtools-1.6/bin:$PATH
samtools view -Sb SRR.sam > SRR.bam
samtools sort SRR.bam >SRR_sorted.bam
bedtools bamtobed -i SRR_sorted.bam>SRR_sorted.bed
#align to 5' end
bedtools genomecov -ibam SRR_sorted.bam -5  -bg  > SRR_5.bedgraph
#align to 3' end
bedtools genomecov -ibam SRR_sorted.bam -3  -bg  > SRR_3.bedgraph
#align to center 
python CenterAlign.py
bedtools sort -i SRR_c.bedgraph > SRR_center.bedgraph
rm SRR_c.bedgraph

#Plot the alignment
Python plot.py



