import sys, numpy as np, scipy.stats as stats
sys.path.append('../')
import fileReader

POP_FILES=['../data_hapmap3/hapmap3_r2_b36_fwd.consensus.qc.poly.chr22_ceu.phased.gz',
           '../data_hapmap3/hapmap3_r2_b36_fwd.consensus.qc.poly.chr22_yri.phased.gz']#,
           #'../data_hapmap3/hapmap3_r2_b36_fwd.consensus.qc.poly.chr22_chd.unr.phased.gz']

F_GM='genetic_map_chr22_b36.txt'
BETA_ALPHA=12
BETA_BETA=3
NOFFSPRING=6
NGENS=8


def poissonMating(pop1, pop2, snpPos, nGens=1, percentPop1=0.2):
    """Caries out mating in similar manner to Hapmix """
    nSNPs, nOffspring=pop1.shape
    outPop=np.empty_like(pop1)
    outPopOrig=np.empty(pop1.shape, dtype=np.byte)
    #Read map and calculate recombination distances
    map=np.array([l.strip().split() for l in open(F_GM).readlines()[1:]], np.float)
    dM=np.diff(map[:,2])/100*nGens  #morgans x  generations
    
    for i in range(nOffspring):
        #alpha=percentPop1
        alpha=stats.beta.rvs(BETA_BETA,BETA_ALPHA)  #Determine percentage of POP1 
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
            
if __name__ == '__main__':
    files=fileReader.concurrentFileReader(POP_FILES[0], POP_FILES[1])
    subjects=files.next()
    snpNames=[]; snpPos=[];  pop1=[];  pop2=[]
    for l in files:
        snpNames.append(l[0])
        snpPos.append(int(l[1]))
        pop1.append(l[2][0])
        pop2.append(l[2][1])
    nSNPs=len(snpNames)
    pop1=np.asarray(pop1)
    pop2=np.asarray(pop2)
    nSNPs, nPop1 = pop1.shape
    nSNPs, nPop2 = pop2.shape
    idxPop1=np.random.permutation(nPop1) #Pick mating pairs
    idxPop2=np.random.permutation(nPop2)
    admixedPop, admixedOrigin=poissonMating(pop1[:,idxPop1[:NOFFSPRING]],
                                            pop2[:,idxPop2[:NOFFSPRING]],
                                            snpPos, NGENS)
    #Save populations
    saveHaplotypes('ancestral1.csv', subjects[0][idxPop1[NOFFSPRING:]], snpNames, 
                   snpPos, pop1[:,idxPop1[NOFFSPRING:]])
    saveHaplotypes('ancestral2.csv', subjects[1][idxPop2[NOFFSPRING:]], snpNames, 
                   snpPos,pop2[:,idxPop2[NOFFSPRING:]])
    saveHaplotypes('admixed.csv', ['ind']*NOFFSPRING, snpNames, snpPos, admixedPop)
    saveHaplotypes('admixed_origin.csv', ['ind']*NOFFSPRING, snpNames, snpPos, admixedOrigin)

