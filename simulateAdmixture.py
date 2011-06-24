import sys
import numpy as np
import scipy.stats as stats
from optparse import OptionParser, OptionValueError

import fileReader
import supportConfig
import popgen

USAGE="""%prog [options] ancestralFile1 ancestralFile2 [...]

where the ancestralFiles contain phased samples in tab delimited
format snps in rows and samples in columns.  One row of column headers
and two columns of header information (i.e. rsId and position).  The
files have to be from the same chromosome indicated by including
chr[0-9]* in the name.
"""
F_GM=supportConfig.createDataPath('data/hapmap2/genetic_map_chr%s_b36.txt')


def poissonMating(pop1, pop2, snpPos, mapFile,  nGens=1, percentPop1=0.2):
    """Depricated use poissonMultiMate instead."""
    return poissonMultiMating([pop1, pop2], snpPos, mapFile, nGens, percentPop1)[0:2]

def poissonMultiMating(pops, snpPos, mapFile, nGens, percentPop1=0.5):
    """Performs mating between individuals in pops using a Poisson model for the recombination events.
    
    Arguments:
    - `pops`: List of two numpy array each with the same nSNPs x nOffspring dimensions
    - `snpPos`: List of positions of all nSNPs basepairs in bp position
    - `mapFile`: File with translations from bp positions to genetic map positions
    - `nGens`: Number of generations to simulate (defualt=1)
    - `percentPop1`: Percentage of alleles from first population in pops (default=.5)
    """
    #Set up output variables
    nPops=len(pops)
    nSNPs, nOffspring=pops[0].shape
    outPop=np.empty_like(pops[0])
    outPopOrig=np.empty(pops[0].shape, dtype=np.byte)
    #Read map and calculate recombination distances
    gm=popgen.geneticMap(mapFile)
    gmPos=gm.pos2gm(snpPos)
    dM=np.diff(gmPos)/100.*nGens #morgans x generations
    #Determine probabilites of differen ancestral populations
    alphas=np.ones(nPops)*(1-percentPop1)/(nPops-1)
    alphas[0]=percentPop1
    print alphas,
    #Step through each offspring and assign haploid genomes
    for i in range(nOffspring):
        recombPos=gmPos[dM>np.random.uniform(size=len(gmPos)-1)] #Determine bp positions of recomb.
        j=0
        for pos in recombPos:  #Step through locations where switch happens
            origin=np.nonzero(np.random.uniform()<=np.cumsum(alphas))[0].min() #0=pop1, 1=pop2 ...
            while gmPos[j]<pos:
                outPop[j, i] = pops[origin][j,i]
                outPopOrig[j,i]=origin
                j+=1
        origin=int(np.nonzero(np.random.uniform()<=np.cumsum(alphas))[0].min()) #0=pop1, 1=pop2 ...
        while j<nSNPs:  #Fill in end of haplotype
            outPop[j, i] = pops[origin][j,i]
            outPopOrig[j,i]=origin
            j+=1
    np.set_printoptions(precision=2)
    print  ' -->', np.asarray([np.sum(outPopOrig==i)/float(outPopOrig.size) for i in range(nPops)]), '%  in', np.sum(np.diff(outPopOrig, axis=0)!=0, axis=0)+1, 'blocks across individuals'
    return outPop, outPopOrig, gmPos


def saveHaplotypes(filename, subjectNames, snpNames, snpPos, snpVals):
    """Saves file in same format at HapMap3 phased data"""
    fp=open(filename, 'w')
    fp.write('rsID	position_b36	%s\n' %'\t'.join(subjectNames))
    for i, name in enumerate(snpNames):
        fp.write('%s\t%s\t' %(name, snpPos[i]))
        for val in snpVals[i,:]:
            fp.write(str(val)+'\t')
        fp.write('\n')
    fp.close()

def saveTped(filename, subjectNames, snpNames, snpPos, mapPos, snpVals, chrom):
    """Saves file in tped and tfam format"""
    nSubs=len(subjectNames)
    with open(filename+'.tfam', 'w') as tFamFp:
        pedPadding="0 %s 0 0 0 0 0"
        for subject in subjectNames[::2]:
            tFamFp.write("%s\n"%pedPadding%subject.replace('_a', ''))
    with open(filename+'.tped', 'w') as fp:
        for i, (rsId, rsPos, rsMap) in enumerate(zip(snpNames, snpPos, gmPos)):
            fp.write('%s %s %0.12g %i' %(chrom, rsId, rsMap, rsPos))
            fp.write(' %s'*nSubs %tuple(snpVals[i, :]))
            fp.write('\n')

def readFiles(files, fileType='beagle', chrom=None):
    nFiles=len(files)
    if fileType=='beagle':
        files=fileReader.concurrentFileReader(*files)
        subjects=files.next()[0]
    elif fileType=='tped':
        tfams=[f.replace('.tped', '.tfam') for f in files]
        tfams=[fileReader.openfile(f) for f in tfams]
        subjects=[]
        for f in tfams:
            subs=[[l.split(None, 1)[0]+'_a',l.split(None, 1)[0]+'_b']  for l in f]
            subjects.append(np.asarray(sum(subs, [])))
        files=fileReader.concurrentFileReader(*files, nHeaders=0, key=[0,1], nLabels=4)
    else:
        sys.stderr.write('ERROR: Filetype has to be either beagle or tped')
        sys.exit()
    snpNames=[]; snpPos=[];  pops=[[] for i in range(nFiles)]
    for s, l in files:
        if fileType=='tped':
            if chrom!=None and chrom!=s[0]:
                continue
            s=[s[1], s[3]]
        snpNames.append(s[0])
        snpPos.append(int(s[1]))
        for i in range(nFiles):
            pops[i].append(l[i])
    nSNPs=len(snpNames)
    pops=map(np.asarray, pops)
    nPops=[l.shape[1] for l in pops]
    return pops,  nPops, subjects, nSNPs, snpPos, snpNames


def permuteIdx(n, nOffspring):
    """Permutes the order of indexes
    
    Arguments:
    - `n`: number of samples to be permuted
    """
    idx=2*np.random.permutation(n/2)  #Permutes every other sample
    idx=np.vstack((idx, idx+1)).flatten(1) #adds matching haploid genome
    return np.hstack((np.random.permutation(idx[:nOffspring]), idx[nOffspring:])) #re-permutes the mate individuals
    


def optSaveFileProcessor(option, opt_str, value, parser):
    try: 
        parser.values.saveFiles=value.split(',')
    except:
        raise OptionValueError('label specification to %s incorrect' %opt_str)


if __name__ == '__main__':
    parser = OptionParser(usage=USAGE)
    parser.add_option('-c', '--chromosome', type='string', dest='chrom', default='22',
                      help='Chromsome being analysed (defualt 22)', metavar='N')
    parser.add_option('-s', '--save', type='string', dest='saveFiles', action='callback', 
                      callback=optSaveFileProcessor, metavar='files',
                      help='prefixes of outfiles (example ceu_chr22,yri_chr22,admixed_ceu_yri_chr22')
    parser.add_option('-n', type='int', dest='nOffspring', default=6,
                      help='number of offspring to create (defualt 6)', metavar='N')
    parser.add_option('-g', '--nGens', type='int', dest='nGens', default=8,
                      help='number of generations since start of admixture', metavar='N')
    parser.add_option('-a', '--alpha', type='float', dest='alpha', default=.5,
                      help='Ratio of populations', metavar='N')
    parser.add_option('-f', '--fileType', type='str', dest='fileType', default='beagle',
                      help='Input Filetype, can be either "beagle" or "tped"')
    parser.add_option('--geneticMap', type='str', dest='geneticMapFile', default=None,
                      help='genetic map for specified chromosome')
    (options, args) = parser.parse_args()

    fileNames=args
    if options.saveFiles==None:
        options.saveFiles=['anc%i_chr%s'%(i+1, options.chrom) for i in range(len(fileNames))]
        options.saveFiles.append('admixed_chr%s'%options.chrom)
    if options.fileType not in ['beagle', 'tped']:
        sys.stderr.write('Input files has to be in plink "tped" or "beagle" meaning simplified haplotype format')
        sys.exit(1)
    if options.geneticMapFile==None: #and options.fileType=='beagle':
        options.geneticMapFile=F_GM%options.chrom

    pops, nPops, subjects, nSNPs, snpPos, snpNames = readFiles(fileNames,
                                                               fileType=options.fileType,
                                                               chrom=options.chrom)
    #determine which will mate with each other.
    indexes=[permuteIdx(n, options.nOffspring) for n in nPops]

    matePairs=[pops[i][:, index[:options.nOffspring]] for i, index in enumerate(indexes)]
    mateSubjects=[subjects[i][index[:options.nOffspring]] for i, index in enumerate(indexes)]
    ancestralPops=[pops[i][:, index[options.nOffspring:]] for i, index in enumerate(indexes)]
    ancestralSubjects=[subjects[i][index[options.nOffspring:]] for i, index in enumerate(indexes)]

    #Running the mating
    admixedPop, admixedOrigin, gmPos=poissonMultiMating(matePairs, snpPos, options.geneticMapFile,
                                                        options.nGens, percentPop1=options.alpha)

    #Save populations
    for i in range(len(nPops)):
        saveTped(options.saveFiles[i], ancestralSubjects[i], snpNames, snpPos, gmPos, ancestralPops[i], options.chrom)
    saveTped(options.saveFiles[-1], ['ind']*options.nOffspring, snpNames, snpPos, gmPos, admixedPop, options.chrom)
    saveTped(options.saveFiles[-1]+'_origin', ['ind']*options.nOffspring, snpNames, snpPos, gmPos, admixedOrigin, options.chrom)

