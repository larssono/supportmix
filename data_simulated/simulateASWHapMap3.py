import sys
import numpy as np
import scipy.stats as stats
sys.path.append('../')
import fileReader

F_CEU='../data_hapmap3/hapmap3_r2_b36_fwd.consensus.qc.poly.chr22_ceu.phased.gz'
F_YRI='../data_hapmap3/hapmap3_r2_b36_fwd.consensus.qc.poly.chr22_yri.phased.gz'
F_GM='genetic_map_chr22_b36.txt'

def poissonMating(pop1, pop2, snpPos, nGens=1, percentPop1=0.2):
    BETA_ALPHA=12
    BETA_BETA=3
    nSNPs, nOffspring=pop1.shape
    outPop=np.empty_like(pop1)
    outPopOrig=np.empty(pop1.shape, dtype=np.byte)
    #Read map and calculate recombination distances
    map=np.array([l.strip().split() for l in open(F_GM).readlines()[1:]], np.float)
    dM=np.diff(map[:,2])/100*nGens  #morgans x  generations
    
    for i in range(nOffspring):
        #alpha=percentPop1
        alpha=stats.beta.rvs(BETA_BETA,BETA_ALPHA)[0]  #Determine percentage of POP1 
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
    fp=open(filename, 'w')
    fp.write('rsID	position_b36	%s\n' %'\t'.join(subjectNames))
    for i, name in enumerate(snpNames):
        fp.write('%s\t%s\t' %(name, snpPos[i]))
        for val in snpVals[i,:]:
            fp.write(str(val)+'\t')
        fp.write('\n')
    fp.close()
            
if __name__ == '__main__':
    nOffspring=12
    nGenerations=8
    files=fileReader.concurrentFileReader(F_CEU, F_YRI)
    subjects=files.next()
    snpNames=[]; snpPos=[];  ceu=[];  yri=[]
    for l in files:
        snpNames.append(l[0])
        snpPos.append(int(l[1]))
        ceu.append(l[2][0])
        yri.append(l[2][1])
    nSNPs=len(snpNames)
    ceu=np.asarray(ceu)
    yri=np.asarray(yri)
    nSNPs, nCEU = ceu.shape
    nSNPs, nYRI = yri.shape
    #Pick mating pairs
    idxCeu=np.random.permutation(nCEU)
    idxYri=np.random.permutation(nYRI)
    ancCeu=ceu[:,idxCeu[:nOffspring]]
    ancYri=yri[:,idxYri[:nOffspring]]
    trainCeu=ceu[:,idxCeu[nOffspring:]]
    trainYri=yri[:,idxYri[nOffspring:]]

    admixedPop, admixedOrigin=poissonMating(ancCeu, ancYri, snpPos, nGenerations)

    #Save populations
    saveHaplotypes('anc_ceu.csv', subjects[0][idxCeu[nOffspring:]], snpNames, snpPos,trainCeu)
    saveHaplotypes('anc_yri.csv', subjects[1][idxYri[nOffspring:]], snpNames, snpPos,trainYri)
    saveHaplotypes('admixed.csv', ['ind']*nOffspring, snpNames, snpPos, admixedPop)
    saveHaplotypes('admixed_origin.csv', ['ind']*nOffspring, snpNames, snpPos, admixedOrigin)




