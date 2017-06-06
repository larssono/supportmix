import  sys,  numpy as np
sys.path.append('..')

def nchoose2(n):
    return n*(n-1)/2

def fst(vals0, vals1, isNorm=True):
    """Calculates fst from two matrices where minorAllele is stored as 1 and major allele as -1"""
    nij=np.asarray([vals0.shape[1], vals1.shape[1]], np.float)
    nAll=nij.sum()
    xi=((vals0==1).sum(1)+(vals1==1).sum(1))/nAll
    xij=[(vals0==1).sum(1)/nij[0], (vals1==1).sum(1)/nij[1]]
    if isNorm:
        top=0; bottom=0
        for j in range(2):
            top+=nchoose2(nij[j])*np.sum(2*nij[j]/(nij[j]-1)*xij[j]*(1-xij[j]))
            bottom+=nchoose2(nij[j])
        return 1-top/bottom/np.sum(2*nAll/(nAll-1)*xi*(1-xi))
    return 1-(np.sum((1-xij[0])*xij[0])+np.sum((1-xij[1])*xij[1]))/2/ np.sum((1-xi)*xi)

def fstFile(file1, file2, isNorm=True):
    """Calculates fst based on data stored in two tab delimited files.
    Each file must contain the following columns: snp names; snp position; snp values. """
    import fileReader
    files=fileReader.concurrentFileReader(file1, file2, key=1)
    subjects=files.next()[0]
    vals=[]
    for ([snpName, snpLocation], snps) in files:
        vals.append(nucleotides2Haplotypes(sum(snps, [])))
    vals=np.asarray(vals, np.float)
    nSamples=len(subjects[0])
    return fst(vals[:,:nSamples], vals[:,nSamples:], isNorm)

def nucleotides2Haplotypes(pop, warnings=True):
    """Given a matrix of nucleotide labels like {'A','T','G','C'}
    returns matrix of -1, 1 encoding for minor and major allele respectively

    params:
       snps: array of characters, e.g.
       [['G', 'A', 'G', 'A', 'A', 'A'],
        ['G', 'A', 'G', 'A', 'T', 'T']] with haplotypes (subjects in across columns)

    return: numpy array of [-1,1]  e.g. [[-1, 1,-1, 1, 1, 1], [-1, 1, -1, 1, 1, 1]]
    """
    if type(pop)==type([]):
        majorAllele=findAlleles(pop, warnings)[0]
        outSNPs=np.ones(len(pop), np.short)
        for i in range(len(pop)):
            outSNPs[i]=(pop[i]==majorAllele)*-2 + 1
        return outSNPs
    majorAllele=[findAlleles(snps)[0] for snps in pop]
    return np.asarray([(subject==majorAllele)*-2+1 for subject in pop.T]).T


 
def nucleotides2SNPs(snps, warnings=True):
    """Given list nucleotide labels returns -1,0,1 using every two
    nucleotides to encode -1 for minor allele homozygote, 0
    heterozygote and 1 for major allele homozygote.
    TODO: Fix this so it also works with matrices

    params:
       snps: list of characters, e.g. ['G', 'A', 'G', 'A', 'A', 'A']

    return: numpy array of [-1,0,1]  e.g. [0, 0 , 1]
    """
    majorAllele,minorAllele=findAlleles(snps, warnings)
    outSNPs=np.ones(len(snps)/2, np.short)
    #Filter and return
    for i in range(0,len(snps),2):
        outSNPs[i/2]=(snps[i]==snps[i+1])*(-2*(snps[i]==majorAllele)+1)
    return outSNPs

class geneticMap(object):
    """keeps track of genetic Map locations and returns closest genetic map
    location given a snp location. """

    def __init__(self,file ):
        """Builds conversion object for mapping base-pair position to genetic map position.

        Parameters:
        - `file` - tab delimited file with three columns: position [bp], combined rate [cM/Mb], Genetic Map [cM].
                   Example files for human genome release 36 can be found at:
                      ftp://ftp.hapmap.org/hapmap/recombination/2008-03_rel22_B36/rates/
                   and for relase 37 at:
                      ftp://ftp.hapmap.org/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz
        """
        fp=open(file)
        self.m=np.asarray([np.asarray(l.split())[[0,2]] for l in fp.readlines()[1:]], np.float)

    def pos2gm(self, pos):
        """Converts position in bp to position in centiMorgans"""
        m=self.m
        pos=np.asarray(pos)
        results=np.empty(pos.shape, dtype=np.float)

        i=m[:,0].searchsorted(pos)  #Find probable locations in map
        i[i==len(m)]=len(m)-1       #Correct those that are beyond the end of the genetic map
        idx=(m[i,0]==pos)|(i==len(m)-1)|(i==0)   #Fill results with those that match exactly
        results[idx]=m[i[idx], 1];
        idx=np.logical_not(idx)
        #Fill those that require interpolation
        results[idx]= (m[i[idx],1]-m[i[idx]-1,1]) / (m[i[idx],0]-m[i[idx]-1,0])*(pos[idx]-m[i[idx]-1,0]) + m[i[idx]-1,1]
        return results



def findAlleles(snps, warnings=True):
    """Given list of nucleotides returns most common and least common"""
    alleles = list(set(snps))
    alleleFreq=[np.sum(snps==nucl) for nucl in alleles]
    idx=np.argsort(alleleFreq)
    if len(idx)>2 and warnings:
        print "There are more than two alleles in the input files: %s" %','.join(alleles)
    elif len(idx)<2:
        return alleles[idx], 'Q'
    return alleles[idx[-1]], alleles[idx[-2]]


if __name__ == '__main__':
    print fstFile(*sys.argv[1:])
