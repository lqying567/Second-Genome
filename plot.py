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

