import sys;sys.path.extend(['..'])
import simulateAdmixture as mate, numpy as np
import popgen 
from variables import *

HGDP_FILE='data/HGDP_raw_data/phasedBeagle/hgdp.%(CHR)s.bgl.phased.gz'%locals()
POPLABEL_FILE='data/HGDP_raw_data/samples_for_nstd27_5.csv'
SMALLEST_POP=20

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
popIdx=np.nonzero(nPops>SMALLEST_POP)[0]  #Populations with enough samples

##############################################################
#Go through and mate pairwise populations and calculate Fst
##############################################################
fpFst=open(FILEFST, 'w')
for  i, idxPop1 in enumerate(popIdx):
    pop1Name=popNames[idxPop1].lower().replace(' ', '_')
    pop1=pops[idxPop1]
    subjects1=subjects[idxPop1]
    mate.saveHaplotypes(FILEANCESTRAL%pop1Name, subjects1[NOFFSPRING*2:], 
                        snpNames, snpPos, pop1[:,NOFFSPRING*2:])
    print '2 Admixed:', pop1Name
    for idxPop2 in popIdx[i+1:]:
        pop2Name=popNames[idxPop2].lower().replace(' ', '_')
        pop2=pops[idxPop2]
        subjects2=subjects[idxPop2]
        admixedPop, admixedOrigin=mate.poissonMating(pop1[:,:NOFFSPRING*2], pop2[:,:NOFFSPRING*2],
                                                     snpPos, mapFile=GM_FILE, nGens=NGENS, percentPop1=.5)
        mate.saveHaplotypes(FILE2ADMPOPS %(pop1Name, pop2Name), ['ind']*2*NOFFSPRING, 
                            snpNames, snpPos, admixedPop)
        mate.saveHaplotypes(FILE2ADMPOPSORIGIN%(pop1Name, pop2Name), 
                            ['ind']*2*NOFFSPRING, snpNames, snpPos, admixedOrigin)
        #Calculate Fst in a slow manner. 
        vals=popgen.nucleotides2Haplotypes(np.hstack((pop1,pop2)))
        fst=popgen.fst(vals[:, :len(subjects1)], vals[:, len(subjects1):])
        fpFst.write('%s\t%s\t%g\n' %(pop1Name, pop2Name,fst))
fpFst.close() 

##############################################################
#Go through and mate threeway between CEU-YRI-X and JPT-CHD-X
##############################################################
popYri=pops[popNames.index('Yoruba')]
popCeu=pops[popNames.index('French')]
popChd=pops[popNames.index('Yoruba')]
popJpn=pops[popNames.index('Bedouin')]

#Mix admixed with third population
for i, (admixedPop1, origin1) in enumerate([mate.poissonMating(popYri[:,:NOFFSPRING*2], popCeu[:,:NOFFSPRING*2], snpPos, mapFile=GM_FILE, nGens=NGENS, percentPop1=.5),
                                            mate.poissonMating(popChd[:,:NOFFSPRING*2], popJpn[:,:NOFFSPRING*2], snpPos, mapFile=GM_FILE, nGens=NGENS, percentPop1=.5)]):
    for idx in popIdx:
        if ((i==0 and idx in [popNames.index('Yoruba'), popNames.index('French')]) or
            (i==1 and idx in [popNames.index('Yoruba'), popNames.index('Bedouin')])):
            continue
        pop3Name=popNames[idx].lower().replace(' ', '_')
        print '3 Admixed:', pop3Name
        pop3=pops[idx]
        admixedPop, admixedOrigin=mate.poissonMating(admixedPop1, pop3[:,:NOFFSPRING*2], snpPos, mapFile=GM_FILE, nGens=NGENS, percentPop1=.7)
        #Fix origin 
        admixedOrigin[admixedOrigin==1]=2
        admixedOrigin[admixedOrigin==0]=origin1[admixedOrigin==0]
        if i==0:
            mate.saveHaplotypes(FILE3ADMPOPSCEUYRI %pop3Name, ['ind']*2*NOFFSPRING, snpNames, snpPos, admixedPop)
            mate.saveHaplotypes(FILE3ADMPOPSCEUYRIORIGIN  %pop3Name, ['ind']*NOFFSPRING, snpNames, snpPos, admixedOrigin)
        elif i==1:
            mate.saveHaplotypes(FILE3ADMPOPSYORBED  %pop3Name, ['ind']*2*NOFFSPRING, snpNames, snpPos, admixedPop)
            mate.saveHaplotypes(FILE3ADMPOPSYORBEDORIGIN %pop3Name, ['ind']*NOFFSPRING, snpNames, snpPos, admixedOrigin)

##############################################################
#Simulate HGDP samples with varying numbers of generations
##############################################################
alpha=0.5
for pop1Name, pop2Name in POPS:
    idxPop1=popNames.index(pop1Name.capitalize())
    idxPop2=popNames.index(pop2Name.capitalize())
    pop1=pops[idxPop1]
    pop2=pops[idxPop2]
    for nGens in SAMPLEGENERATIONS:
        print 'Generations', pop1Name, pop2Name, nGens
        admixedPop, admixedOrigin=mate.poissonMating(pop1[:,:NOFFSPRING*2], pop2[:,:NOFFSPRING*2],
                                                     snpPos, mapFile=GM_FILE,nGens=nGens, percentPop1=alpha)
        mate.saveHaplotypes(FILEGENSADMIX%(pop1Name, pop2Name, nGens), ['ind']*2*NOFFSPRING, 
                            snpNames, snpPos, admixedPop)
        mate.saveHaplotypes(FILEGENSADMIXORIGIN%(pop1Name, pop2Name, nGens), 
                            ['ind']*2*NOFFSPRING, snpNames, snpPos, admixedOrigin)
        if nGens==5:
            for a in ALPHAS:
                admixedPop, admixedOrigin=mate.poissonMating(pop1[:,:NOFFSPRING*2], pop2[:,:NOFFSPRING*2],
                                                             snpPos, mapFile=GM_FILE,nGens=nGens, percentPop1=a)
                mate.saveHaplotypes(FILEALPHAADMIX%(pop1Name, pop2Name, nGens,a), 
                                    ['ind']*2*NOFFSPRING, snpNames, snpPos, admixedPop)
                mate.saveHaplotypes(FILEALPHAADMIXORIGIN%(pop1Name, pop2Name, nGens,a), 
                                    ['ind']*2*NOFFSPRING, snpNames, snpPos, admixedOrigin)
