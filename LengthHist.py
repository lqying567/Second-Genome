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
