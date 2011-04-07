import gzip
import numpy as np
from collections import deque
from popgen import nucleotides2Haplotypes,nucleotides2SNPs

def openfile(file):
    """Returns a file object by opening either a gz or uncompressed file. """
    if file.endswith('.gz'):
        return gzip.open(file)
    else:
        return open(file)


def csv2Dict(file, sep='\t'):
    """Given a tab delimited file returns the first column as list and
    remaining columns as dictionary of lists.
    Parameters:
    `file` input file
    `sep` - optional parameter incase file is not delimited by tab
    """
    fp=openfile(file)
    head=fp.readline().rstrip().split(sep)
    head=[l.strip('"') for l in head]
    head= dict(zip(head[1:], range(1, len(head)))) #translation between names and columns
    vals=[line.split(sep) for line in fp]
    phenotypes={}
    subjects=[val[0].strip('"') for val in vals]
    for key, i in head.items():
        phenotypes[key]=np.asarray([val[i].strip('"') for val in vals]) 
    return subjects, phenotypes

def csv2DictDict(file, sep='\t'):
    """Given a tab delimited file returns a dict of dicts where the
      first key references a column and the second key the row.

      Example:"  id       age    gender
                 NML-001   24     M
                 NML-002   48     F"

          fileDict=csv2DictDict(file)
          fileDict['NML-001']['age']  -> 24
          fileDict['NML-002']['gender'] -> F"""    
    fp=openfile(file)
    head=fp.readline().rstrip().split(sep)
    head=[l.strip('"') for l in head]
    head= dict(zip(head[1:], range(1, len(head)))) #translation between names and columns
    phenotypes={}
    for line in fp:
        line = line.split(sep) 
        subject=line[0].strip('"')
        phenotypes[subject]={}
        for key, i in head.items(): 
            phenotypes[subject][key]=line[i].strip('"')
    return phenotypes


        
def concurrentFileReader(*args, **kwargs):
    """Given a list of files returns common lines one at a time as determined by
    first word (or key word)
    
    First call returns the first line, i.e the column headers.
    Subsequent calls returns one line at a time where row headers are similar."""
    fps= map(openfile, args)  #open input files
    key= kwargs.get('key', 0)
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
            multiReadLine(fps, lineDeques, lineLabels, key)
            foundRow = findCommonRow(lineLabels)
            while foundRow=='':   #not found common row
                multiReadLine(fps, lineDeques, lineLabels, key)
                foundRow = findCommonRow(lineLabels)
            out=[]
            for fileDict, fileDeque in zip(lineLabels, lineDeques):  #Output the common value
                line = fileDeque.popleft()
                del fileDict[line.split(None, key+1)[key]]
                while line.split(None, key+1)[key] != foundRow:
                    line = fileDeque.popleft()
                    del fileDict[line.split(None, key+1)[key]]
                out.append(line)
            out = [l.split() for l in out]     #Split line into parts
            snpNames=out[0][0]
            snpLocations=out[0][1]
            snps=[l[2:] for l in out]          #Extract values
            yield snpNames, snpLocations, snps 
    except EOFError:
        return
        
    
def multiReadLine(fps, lineDeques, lineLabels, key):
    """Reads one line from each file in fps and stores at end of each
    list stored in lineDeques.  Raises EOFError when one file reaches its end.
    """
    nFilesAtEnd=0
    for i, fp in enumerate(fps): #Read next line
        str=fp.readline()
        if str != '':
            lineDeques[i].append(str.strip())
            lineLabels[i][str.split(None, key+1)[key]] = 1
        else:
            nFilesAtEnd+=1
    if nFilesAtEnd==len(fps) or (np.array(map(len, lineDeques))==0).any():
        raise EOFError


def findCommonRow(lineLabels):
    for label in lineLabels[0]:
        if sum([label in otherLabels for otherLabels in lineLabels[1:]])==len(lineLabels)-1:
            return label
    return ''


if __name__ == '__main__':
    f=concurrentFileReader('file1.gz', 'file2.gz', 'file3.gz')
    for l in f:
        print l


