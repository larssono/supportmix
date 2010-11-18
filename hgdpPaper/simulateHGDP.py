import sys;sys.path.extend(['..'])
import simulateAdmixture as mate, numpy as np
import popgen 

CHR='chr1'
HGDP_FILE='../../../human_genome_data/HumanPlinkData/HGDP/phasedBeagle/hgdp.%(CHR)s.bgl.phased.gz'%locals()
GM_FILE='../../../human_genome_data/hapmap2/genetic_map_%(CHR)s_b36.txt'%locals()
POPLABEL_FILE='../../../human_genome_data/HumanPlinkData/HGDP/samples_for_nstd27_5.csv'
SMALLEST_POP=20
NGENS=5
NOFFSPRING=2

#Create map from name to population
popDict=dict([(l.split(',')[2], l.split(',')[7]) for l in open(POPLABEL_FILE)])
popDict.pop('Sample ID')

#Read monlithic file containing all populations 
pops, nPops, subjects, nSNPs, snpPos, snpNames = mate.readFiles([HGDP_FILE])

#Divide said file into individual populations
popNames=[]
for pop in set(popDict.values()):
    indexes=[]
    for i, sub in enumerate(subjects[0]):
        if popDict[sub[:-2]] == pop:
            indexes.append(i)
    popNames.append(pop)
    pops.append(pops[0][:, indexes])
    subjects.append(subjects[0][indexes])
subjects.pop(0) #Remove unsorted subjects
pops.pop(0)     #Remove unsorted populations
nPops=np.asarray([l.shape[1] for l in pops])


#Go through and mate pairwise populations
popIdx=np.nonzero(nPops>SMALLEST_POP)[0]  #Populations with enough samples
for  i, idxPop1 in enumerate(popIdx):
    pop1Name=popNames[idxPop1].lower().replace(' ', '_')
    pop1=pops[idxPop1]
    subjects1=subjects[idxPop1]
    mate.saveHaplotypes('hgdp/ancestral_hgdp_%s.%s.csv'%(pop1Name, CHR), subjects1[NOFFSPRING*2:], 
                        snpNames, snpPos, pop1[:,NOFFSPRING*2:])
    for idxPop2 in popIdx[i+1:]:
        pop2Name=popNames[idxPop2].lower().replace(' ', '_')
        pop2=pops[idxPop2]
        subjects2=subjects[idxPop2]
        admixedPop, admixedOrigin=mate.poissonMating(pop1[:,:NOFFSPRING*2], pop2[:,:NOFFSPRING*2],
                                                     snpPos, mapFile=GM_FILE,nGens=NGENS, percentPop1=.5)
        mate.saveHaplotypes('hgdp/admixed_hgdp_%s_%s.%s.csv'%(pop1Name, pop2Name, CHR), ['ind']*2*NOFFSPRING, 
                            snpNames, snpPos, admixedPop)
        mate.saveHaplotypes('hgdp/admixed_origin_hgdp_%s_%s.%s.csv'%(pop1Name, pop2Name, CHR), 
                            ['ind']*2*NOFFSPRING, snpNames, snpPos, admixedOrigin)

        vals=popgen.nucleotides2Haplotypes(np.hstack((pop1,pop2)))
        fst=popgen.fst(vals[:, :len(subjects1)], vals[:, len(subjects1):])
        print 'ancestral_hgdp_%s.%s.csv\tancestral_hgdp_%s.%s.csv\t%g' %(pop1Name,CHR,  pop2Name, CHR,fst)
 
