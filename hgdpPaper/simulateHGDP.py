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
OUTFILE_FST='hgdp_fst.csv'

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

# ##############################################################
# #Go through and mate pairwise populations and calculated Fst
# ##############################################################
# fpFst=open(OUTFILE_FST, 'w')
# for  i, idxPop1 in enumerate(popIdx):
#     pop1Name=popNames[idxPop1].lower().replace(' ', '_')
#     pop1=pops[idxPop1]
#     subjects1=subjects[idxPop1]
#     mate.saveHaplotypes('hgdp_ancestral/ancestral_hgdp_%s.%s.csv'%(pop1Name, CHR), subjects1[NOFFSPRING*2:], 
#                         snpNames, snpPos, pop1[:,NOFFSPRING*2:])
#     for idxPop2 in popIdx[i+1:]:
#         pop2Name=popNames[idxPop2].lower().replace(' ', '_')
#         pop2=pops[idxPop2]
#         subjects2=subjects[idxPop2]
#         admixedPop, admixedOrigin=mate.poissonMating(pop1[:,:NOFFSPRING*2], pop2[:,:NOFFSPRING*2],
#                                                      snpPos, mapFile=GM_FILE,nGens=NGENS, percentPop1=.5)
#         mate.saveHaplotypes('hgdp2/admixed_hgdp_%s_%s.%s.csv'%(pop1Name, pop2Name, CHR), ['ind']*2*NOFFSPRING, 
#                             snpNames, snpPos, admixedPop)
#         mate.saveHaplotypes('hgdp2/admixed_origin_hgdp_%s_%s.%s.csv'%(pop1Name, pop2Name, CHR), 
#                             ['ind']*2*NOFFSPRING, snpNames, snpPos, admixedOrigin)
#         #Calculate Fst in a slow manner. 
#         vals=popgen.nucleotides2Haplotypes(np.hstack((pop1,pop2)))
#         fst=popgen.fst(vals[:, :len(subjects1)], vals[:, len(subjects1):])
#         fpFst.write('ancestral_hgdp_%s.%s.csv\tancestral_hgdp_%s.%s.csv\t%g\n' %(pop1Name,CHR,  pop2Name, CHR,fst))
# fpFst.close() 

##############################################################
#Go through and mate threeway between CEU-YRI-X and JPT-CHD-X
##############################################################
popYri=pops[popNames.index('Yoruba')]
popCeu=pops[popNames.index('French')]
popChd=pops[popNames.index('Han')]
popJpn=pops[popNames.index('Japanese')]
#Mix admixed with third population
for i, (admixedPop1, origin1) in enumerate([mate.poissonMating(popYri[:,:NOFFSPRING*2], popCeu[:,:NOFFSPRING*2], snpPos, mapFile=GM_FILE, nGens=NGENS, percentPop1=.5),
                                          mate.poissonMating(popChd[:,:NOFFSPRING*2], popJpn[:,:NOFFSPRING*2], snpPos, mapFile=GM_FILE, nGens=NGENS, percentPop1=.5)]):
    for idx in popIdx:
        if ((i==0 and idx in [popNames.index('Yoruba'), popNames.index('French')]) or
            (i==1 and idx in [popNames.index('Han'), popNames.index('Japanese')])):
            continue
        pop3Name=popNames[idx].lower().replace(' ', '_')
        print pop3Name
        pop3=pops[idx]
        admixedPop, admixedOrigin=mate.poissonMating(admixedPop1, pop3[:,:NOFFSPRING*2], snpPos, mapFile=GM_FILE, nGens=NGENS, percentPop1=.7)
        #Fix origin 
        admixedOrigin[admixedOrigin==1]=2
        admixedOrigin[admixedOrigin==0]=origin1[admixedOrigin==0]
        if i==0:
            mate.saveHaplotypes('hgdp3/admixed_hgdp_yoruba_french_%s.%s.csv' %(pop3Name, CHR), ['ind']*2*NOFFSPRING, snpNames, snpPos, admixedPop)
            mate.saveHaplotypes('hgdp3/admixed_hgdp_origin_yoruba_french_%s.%s.csv' %(pop3Name, CHR), ['ind']*NOFFSPRING, snpNames, snpPos, admixedOrigin)
        elif i==1:
            mate.saveHaplotypes('hgdp3/admixed_hgdp_han_japanese_%s.%s.csv' %(pop3Name, CHR), ['ind']*2*NOFFSPRING, snpNames, snpPos, admixedPop)
            mate.saveHaplotypes('hgdp3/admixed_hgdp_origin_han_japanese_%s.%s.csv' %(pop3Name, CHR), ['ind']*NOFFSPRING, snpNames, snpPos, admixedOrigin)

asdf
##############################################################
#Simulate HGDP samples with varying numbers of generations
##############################################################
popIdx=[popNames.index(name) for name in ['Yoruba', 'French', 'Han', 'Bedouin']]
alpha=0.5
for  i, idxPop1 in enumerate(popIdx):
    pop1Name=popNames[idxPop1].lower().replace(' ', '_')
    pop1=pops[idxPop1]
    for idxPop2 in popIdx[i+1:]:
        pop2Name=popNames[idxPop2].lower().replace(' ', '_')
        pop2=pops[idxPop2]
        for nGens in [1,5,10,50,100,200]:
            admixedPop, admixedOrigin=mate.poissonMating(pop1[:,:NOFFSPRING*2], pop2[:,:NOFFSPRING*2],
                                                         snpPos, mapFile=GM_FILE,nGens=nGens, percentPop1=alpha)
            mate.saveHaplotypes('hgdp_generations/admixed_hgdp_%s_%s_%ggens.%s.csv'%(pop1Name, pop2Name, nGens, CHR), ['ind']*2*NOFFSPRING, 
                                snpNames, snpPos, admixedPop)
            mate.saveHaplotypes('hgdp_generations/admixed_origin_hgdp_%s_%s_%ggens.%s.csv'%(pop1Name, pop2Name, nGens, CHR), 
                                ['ind']*2*NOFFSPRING, snpNames, snpPos, admixedOrigin)


#simulate two populations with different alpha
