import gzip
import numpy as np
from collections import deque

class snpData(object):
    """Keeps track of SNP data for subjects or Haplotypes."""
    
    def __init__(self, ):
        """
        """
        
def readHapMapFile(file):
    """Reads first and second column of file storing SNP names and position in bp
    Return: subjectList, snpNames, snpLocations, snps 
    """
    fp=gzip.open(file); 
    subjectList=np.asarray(fp.readline().strip().split()[2:], dtype=np.str_)
    snps = np.asarray([l.split() for l in fp], dtype=np.str_)
    snpNames=snps[:,0]
    snpLocations=snps[:,1].astype(int)
    snps=snps[:,2:]
    return subjectList, snpNames, snpLocations, snps 

def concurrentFileReader(*args):
    """Given a list of files returns common lines one at a time.
    
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
            for fileDeque in lineDeques:  #Output the common value
                line = fileDeque.popleft()
                while not line.startswith(foundRow):
                    line = fileDeque.popleft()
                out.append(line)
            print map(len, lineDeques),
            out = [l.split() for l in out]     #Split line into parts
            snpNames=[l[0] for l in out]       #Extract row headers
            snpLocations=[l[1] for l in out]
            snps=[l[2:] for l in out]          #Extract values
            yield snpNames, snpLocations, snps 
    except EOFError:
        print "eof reached"
        return
        
    
def multiReadLine(fps, lineDeques):
    """Reads one line from each file in fps and stores at end of each
    list stored in lineDeques.  Raises EOFError when one file reaches its end.
    """
    nFilesAtEnd=0
    for i, fp in enumerate(fps): #Read next line
        str=fp.readline()
        if str != '':
            lineDeques[i].append(str.strip())
        else:
            nFilesAtEnd+=1
    if nFilesAtEnd==len(fps) or (np.array(map(len, lineDeques))==0).any():
        raise EOFError
    


def findCommonRow(lineDeques):
    for line in lineDeques[0]:
        print len(lineDeques[0]),
        try:
            label = line.split(None, 1)[0]
            found=1
            for otherDeque in lineDeques[1:]:
                for otherLine in otherDeque:
                    if  label == otherLine.split(None, 1)[0]:
                        found +=1
                        break
            if found==len(lineDeques):
                return label
        except IndexError:
            return ''
    return ''

    
def nucleotides2SNPs(snps):
    """Recieves an array of nucleotiedes and converts to numpy array of 0,1"""
    nSNPs, nSubj = snps.shape
    out=np.zeros((nSNPs, nSubj), dtype=int)
#     translation={'A':4,'C':3,'G':2, 'T':1}
#     for i in range(nSNPs):
#         for j in range(nSubj):
#             out[i,j]=translation[snps[i,j]]
#     return out
    
    for i, snp in enumerate(snps):
        genotypes=list(set(snp))
        if (snp==genotypes[0]).sum() >nSubj/2.:
            majorAllele=genotypes[0]
        else:
            majorAllele=genotypes[1]            
        out[i,:]=np.asarray(snp)==majorAllele
    out[out==0] = -1    
    return out

def findCommonIdx(snpNames1, snpNames2 ):
    """
    """
    commonSNPs=set(snpNames1).intersection(snpNames2)
    idx1=np.zeros(len(commonSNPs), np.int); idx2=np.zeros(len(commonSNPs), np.int)
    for i, snp in enumerate(commonSNPs):
        idx1[i] = np.nonzero(snpNames1==snp)[0]
        idx2[i] = np.nonzero(snpNames2==snp)[0]
    return idx1, idx2


if __name__ == '__main__':
    f=concurrentFileReader('file1.gz', 'file2.gz', 'file3.gz')
    for l in f:
        print l
