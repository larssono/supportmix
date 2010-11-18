import pylab, fileReader, sys, numpy as np
from scipy.linalg import svd

files=fileReader.concurrentFileReader(*sys.argv[1:])

subjects=files.next()
snpLabels=[]        #stores snp labels from in files
snpLocations=[]     #stores physical location from files
snpVals=[]
for i, (snpName, snpLocation, snps) in enumerate(files):
    snpLabels.append(snpName)
    snpLocations.append(float(snpLocation))
    snpVals.append(fileReader.nucleotides2Haplotypes(sum(snps, [])))

snps=np.asarray(snpVals)

[u,s,vt]=svd(snps,0)

nPops=map(len, subjects)
colors=pylab.cm.copper(np.linspace(0,1,len(subjects)))
colors[-1,:]=[1,0,0,1]

idx0=0
for i, sub in enumerate(subjects):
    idx1=len(sub)+idx0
    pylab.plot(vt[1,idx0:idx1], vt[2,idx0:idx1], '.', markersize=10, color=colors[i,:]) 
    idx0=idx1
pylab.xlabel('PC1')
pylab.ylabel('PC2')
pylab.legend([l.split('.')[0] for l in  sys.argv[1:]])
pylab.show()
