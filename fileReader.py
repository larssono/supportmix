import gzip
import numpy as np
from collections import deque
        
def concurrentFileReader(*args):
    """Given a list of files returns common lines one at a time as determined by first word.
    
    First call returns the first line, i.e the column headers.
    Subsequent calls returns one line at a time where row headers are similar."""
    fps= map(gzip.open, args)  #open input files
    #fps= map(open, args)  #open input files
    lineDeques=[]  #Create storage for read lines
    lineLabels=[]  #Create storage for labels in readLines

    for i, fp in enumerate(fps): 
        lineDeques.append(deque())
        lineLabels.append(dict())
    
    subjectLists=[np.asarray(fp.readline().strip().split()[2:], dtype=np.str_) for fp in fps]
    yield subjectLists   #First time called give subjects

    try:
        while True:
            multiReadLine(fps, lineDeques, lineLabels)
            foundRow = findCommonRow(lineLabels)
            while foundRow=='':   #not found common row
                multiReadLine(fps, lineDeques, lineLabels)
                foundRow = findCommonRow(lineLabels)
            out=[]
            for fileDict, fileDeque in zip(lineLabels, lineDeques):  #Output the common value
                line = fileDeque.popleft()
                del fileDict[line.split(None, 1)[0]]
                while not line.startswith(foundRow):
                    line = fileDeque.popleft()
                    del fileDict[line.split(None, 1)[0]]
                out.append(line)
            out = [l.split() for l in out]     #Split line into parts
            snpNames=out[0][0]
            snpLocations=out[0][1]
            snps=[l[2:] for l in out]          #Extract values
            yield snpNames, snpLocations, snps 
    except EOFError:
        return
        
    
def multiReadLine(fps, lineDeques, lineLabels):
    """Reads one line from each file in fps and stores at end of each
    list stored in lineDeques.  Raises EOFError when one file reaches its end.
    """
    nFilesAtEnd=0
    for i, fp in enumerate(fps): #Read next line
        str=fp.readline()
        if str != '':
            lineDeques[i].append(str.strip())
            lineLabels[i][str.split(None, 1)[0]] = 1
        else:
            nFilesAtEnd+=1
    if nFilesAtEnd==len(fps) or (np.array(map(len, lineDeques))==0).any():
        raise EOFError
    
def findCommonRow(lineLabels):
    for label in lineLabels[0]:
        if sum([label in otherLabels for otherLabels in lineLabels[1:]])==len(lineLabels)-1:
            return label
    return ''

ALLELES={0:'A', 1:'T', 2:'G', 3:'C'}
def __findAlleles(snps):
    """Given list of nucleotides returns most common and least common"""
    alleleFreq=[snps.count(nucl) for nucl in ['A','T','G','C']]
    idx=np.argsort(alleleFreq)
    return ALLELES[idx[-2]], ALLELES[idx[-1]]

def nucleotides2Haplotypes(snps):
    """Given a list of nucleotide labels returns -1, 1
    encoding for major and minor allele.

    params:
       snps: list of characters, e.g. ['G', 'A', 'G', 'A', 'A', 'A']

    return: numpy array of [-1,1]  e.g. [-1, 1,-1, 1, 1, 1]
    """
    majorAllele,minorAllele=__findAlleles(snps)
    #Filter and return
    outSNPs=np.ones(len(snps), np.short)
    for i in range(len(snps)):
        outSNPs[i]=(snps[i]==majorAllele)*-2 + 1
    return outSNPs

def nucleotides2SNPs(snps):
    """Given list nucleotide labels returns -1,0,1 using every two
    nucleotides to encode -1 for minor allele homozygote, 0
    heterozygote and 1 for major allele homozygote.

    params:
       snps: list of characters, e.g. ['G', 'A', 'G', 'A', 'A', 'A']

    return: numpy array of [-1,0,1]  e.g. [0, 0 , 1]
    """
    majorAllele,minorAllele=__findAlleles(snps)
    outSNPs=np.ones(len(snps)/2, np.short)
    #Filter and return
    for i in range(0,len(snps),2):
        outSNPs[i/2]=(snps[i]==snps[i+1])*(-2*(snps[i]==majorAllele)+1)
    return outSNPs


if __name__ == '__main__':
    f=concurrentFileReader('file1.gz', 'file2.gz', 'file3.gz')
    for l in f:
        print l


