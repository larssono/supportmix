import gzip
import numpy as np
from collections import deque, defaultdict
from popgen import nucleotides2Haplotypes, nucleotides2SNPs

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
    first word (or key word(s))

    First call returns the first line, i.e the column headers.
    Subsequent calls returns one line at a time where row headers are similar.

    Parameters:
    `args` list of fileNames 
    `key` - optional parameter of column number to use to match files on (default 0)
            can also be a list of columns to use
    `nHeaders` - optional parameters for the number of header rows output 
                 before the actual data (default 1)
    `nLabels` - optional parameter of number of columns with labels returned as 
                sepereate paramters when generator is called (default 2)

    Returns:
       1) First invocation returns a list of nHeaders (default 1) items -- one item
          per header line in file.
       2) consequitive invocations return
           labels - tuple of first nLabels (default 2) column values 
           data - list of data values in columns nLabels+1 to end of line
    """
    fps= map(openfile, args)  #open input files
    key= kwargs.get('key', [0])
    nHeaders= kwargs.get('nHeaders', 1)
    nLabels= kwargs.get('nLabels', 2)
    if type(key)==type(0):
        key=[key]
    maxKEY=max(key)+1
    lineDeques=[]  #Create storage for read lines
    lineLabels=defaultdict(int)  #Create storage for labels in readLines
    nFILES=len(fps)
    for i, fp in enumerate(fps): 
        lineDeques.append(deque())
    
    headerList=[]
    for i in range(nHeaders):
        headerList.append([np.asarray(fp.readline().strip().split()[nLabels:], dtype=np.str_) for fp in fps])
    if len(headerList)>0: yield headerList   #First time called give subjects

    try:
        while True:
            foundRow=multiReadLine(fps, lineDeques, lineLabels, key, nFILES, maxKEY)
            while foundRow=='':   #no common row read yet
                foundRow=multiReadLine(fps, lineDeques, lineLabels, key, nFILES, maxKEY)
            out=[]
            for i, fileDeque in enumerate(lineDeques):  #Output the common value
                line = fileDeque.popleft()
                label='_'.join([item for k,item in enumerate(line.split(None, maxKEY)) if k in key])
                if i==0: del lineLabels[label]
                while label != foundRow:
                    line = fileDeque.popleft()
                    label='_'.join([item for k,item in enumerate(line.split(None, maxKEY)) if k in key])
                    if i==0: del lineLabels[label]
                out.append(line)
            out = [l.split() for l in out]     #Split line into parts
            labelRows=out[0][0:nLabels]
            data=[l[nLabels:] for l in out]          #Extract values
            yield labelRows, data
    except EOFError:
        return
        
    
def multiReadLine(fps, lineDeques, lineLabels, key, nFILES, maxKEY):
    """Reads one line from each file in fps and stores at end of each
    list stored in lineDeques.  Raises EOFError when one file reaches its end.

    Returns the row value at the key position that is present in all files
    """
    found=''
    nFilesAtEnd=0
    for i, fp in enumerate(fps): #Read next line
        line=fp.readline()
        if line != '':
            lineDeques[i].append(line.strip())
            label='_'.join([item for k,item in enumerate(line.split(None, maxKEY)) if k in key])
            lineLabels[label] +=1
            if lineLabels[label]==nFILES:
                found=label
        else:
            nFilesAtEnd+=1
    if nFilesAtEnd==len(fps) or (np.array(map(len, lineDeques))==0).any():
        raise EOFError
    return found


if __name__ == '__main__':
    f=concurrentFileReader('/home/lom/current_projects/human_genome_data/HumanPhasedData/hgdp/hgdp.chr22.bgl.phased.gz',
                           '/home/lom/current_projects/human_genome_data/HumanPhasedData/qatar/qatarFlt.chr22.bgl.phased.gz',
                           '/home/lom/current_projects/human_genome_data/HumanPhasedData/hapmap3/ceu/hapmap3_r2_b36_fwd.consensus.qc.poly.chr22_ceu.phased.gz')
    print f.next()
    for s, d in f:
        print (s[0], s[1], d)

