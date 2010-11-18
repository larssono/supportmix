import sys, numpy as np, scipy.stats as stats
sys.path.append('../')
import fileReader

def poissonMating(pop1, pop2, snpPos, mapFile,  nGens=1, percentPop1=0.2):
    """Caries out mating in similar manner to Hapmix """
    nSNPs, nOffspring=pop1.shape
    outPop=np.empty_like(pop1)
    outPopOrig=np.empty(pop1.shape, dtype=np.byte)
    #Read map and calculate recombination distances
    map=np.array([l.strip().split() for l in open(mapFile).readlines()[1:]], np.float)
    dM=np.diff(map[:,2])/100*nGens  #morgans x  generations
    
    for i in range(nOffspring):
        alpha=percentPop1
        #alpha=stats.beta.rvs(BETA_BETA,BETA_ALPHA)  #Determine percentage of POP1 
        recombPos=map[dM>np.random.uniform(size=len(map)-1), 0] #Determine bp positions of recomb.
        j=0
        for pos in recombPos:  #Step through locations where switch happens
            origin=int(np.random.uniform()>alpha)  #0 = pop1, 1=pop2
            while snpPos[j]<pos:  
                if origin==0:
                    outPop[j, i] = pop1[j,i]
                    outPopOrig[j,i]=0
                else:
                    outPop[j, i] = pop2[j,i]
                    outPopOrig[j,i]=1
                j+=1
        origin=int(np.random.uniform()>alpha)  #0 = pop1, 1=pop2
        while j<nSNPs:  #Fill in end of haplotype
            if origin==0:
                outPop[j, i] = pop1[j,i]
                outPopOrig[j,i]=0
            else:
                outPop[j, i] = pop2[j,i]
                outPopOrig[j,i]=1
            j+=1
    return outPop, outPopOrig

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

def readFiles(files):
    nFiles=len(files)
    files=fileReader.concurrentFileReader(*files)
    subjects=files.next()
    snpNames=[]; snpPos=[];  pops=[[] for i in range(nFiles)]
    for l in files:
        snpNames.append(l[0])
        snpPos.append(int(l[1]))
        for i in range(nFiles):
            pops[i].append(l[2][i])
    nSNPs=len(snpNames)
    pops=map(np.asarray, pops)
    nPops=[l.shape[1] for l in pops]
    return pops,  nPops, subjects, nSNPs, snpPos, snpNames

if __name__ == '__main__':
    CHR='chr22'
    DATA_DIR='../../../human_genome_data/'
    POP_FILES=[DATA_DIR+'hapmap3/CEU/TRIOS/hapmap3_r2_b36_fwd.consensus.qc.poly.%s_ceu.phased.gz'%CHR,
               DATA_DIR+'hapmap3/YRI/TRIOS/hapmap3_r2_b36_fwd.consensus.qc.poly.%s_yri.phased.gz'%CHR]
    F_GM=DATA_DIR+'/hapmap2/genetic_map_%s_b36.txt'%CHR
    BETA_ALPHA=12
    BETA_BETA=3
    NOFFSPRING=6
    NGENS=8

    pops, nPops, subjects, nSNPs, snpPos, snpNames = readFiles(POP_FILES)
    pop1, pop2=pops
    nPop1, nPop2=nPops
    idxPop1=np.random.permutation(nPop1) #Pick mating pairs
    idxPop2=np.random.permutation(nPop2)
    admixedPop, admixedOrigin=poissonMating(pop1[:,idxPop1[:NOFFSPRING]],
                                            pop2[:,idxPop2[:NOFFSPRING]],
                                            snpPos, F_GM, NGENS)
    #Save populations
    saveHaplotypes('ancestral_ceu.%s.csv'%CHR, subjects[0][idxPop1[NOFFSPRING:]], snpNames, 
                   snpPos, pop1[:,idxPop1[NOFFSPRING:]])
    saveHaplotypes('ancestral_yri.%s.csv'%CHR, subjects[1][idxPop2[NOFFSPRING:]], snpNames, 
                   snpPos,pop2[:,idxPop2[NOFFSPRING:]])
    saveHaplotypes('admixed_ceu_yri.%s.csv'%CHR, ['ind']*NOFFSPRING, snpNames, snpPos, admixedPop)
    saveHaplotypes('admixed_ceu_yri_origin.%s.csv'%CHR, ['ind']*NOFFSPRING, snpNames, snpPos, admixedOrigin)

