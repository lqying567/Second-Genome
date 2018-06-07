#!/bin/bash
echo hello world
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
#Creat a seq length distribution histgram for a fastq file
import numpy as np
import os
import matplotlib.pyplot as plt
from Bio import SeqIO
os.chdir('/Users/lanqingying/Documents/SRA/Bowtie2')
s=[];t=[];
handle=open('Test_removed.fastq')
while True:
    try:
        seq_id=handle.next().strip("\n")
        seq=handle.next().strip('\n')
        index=handle.next().strip('\n')
        quality=handle.next().strip('\n')
        s.append(len(seq))
    except StopIteration:
        break
handle.close()

#for record in SeqIO.parse('input.fasq','fastq'):
 #   t.append(len(record.seq))

p=plt.hist(s,bins=range(min(s),max(s)+1))
plt.xlabel('Length')
plt.show(p)

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
#Creat a read alignment plot
import numpy as np
import matplotlib.pyplot as plt
import os
os.chdir('/Users/lanqingying/Documents/SRA/sra')
n=5000000
bigl=[0]*n
#bigl=np.zeros(n,dtype=np.int)
lf=100; rt=50
with open('SRR1734430_5.bedGraph') as data:
    for line in data:
        s=line.split('\t')
        bigl[int(s[1])]=float(s[3])
arr=np.empty((0,lf+rt+1),float)       
with open('Annotation_sorted.gtf') as f:
    for line in f:
        s=line.split('\t')
        if s[2]=='CDS':
            if s[6]=='+':   #plus strand
                start=int(s[3]); stop=int(s[4])
                l=bigl[(stop-lf):(stop+rt+1)]

            else:           #minus strand
                start=int(s[4]); stop=int(s[3])
                l=bigl[(stop-rt):(stop+lf+1)][::-1]

            if sum(l)>0:
                    l=np.array(l)/float(sum(l))
                    arr=np.append(arr,[l],axis=0)
        else:
            continue
       
x=range(-lf,rt+1)
y=np.mean(arr,axis=0)
plt.plot(x,y,marker='.',linestyle='-')
plt.xlabel('Distance to stop codon')
plt.show()




