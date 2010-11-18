import simulateAdmixture as mate, numpy as np

DATA_DIR='../../../human_genome_data/'
CHR='chr22'
POP1='yri'
POP2='mkk'

POP_FILES=(DATA_DIR+'data_hapmap3/hapmap3_r2_b36_fwd.consensus.qc.poly.%s_%s.phased.gz'%(CHR,POP1), 
           DATA_DIR+'data_hapmap3/hapmap3_r2_b36_fwd.consensus.qc.poly.%s_%s.phased.gz'%(CHR, POP2))
F_GM=DATA_DIR+'genetic_map_%s_b36.txt'%CHR
NGENS=8
NOFFSPRING=6

#Read populations
pops, nPops, subjects, nSNPs, snpPos, snpNames = mate.readFiles(POP_FILES)
idxPop1=np.random.permutation(nPops[0]) #Pick mating pairs
idxPop2=np.random.permutation(nPops[1])
#Mix the first two populations
admixedPop, admixedOrigin=mate.poissonMating(pops[0][:,idxPop1[:NOFFSPRING]],
                                             pops[1][:,idxPop2[:NOFFSPRING]],
                                             snpPos, NGENS, .5)

# #Save populations
mate.saveHaplotypes('ancestral_%s.%s.csv'%(POP1, CHR), subjects[0][idxPop1[NOFFSPRING:-10]], snpNames, 
                    snpPos, pops[0][:,idxPop1[NOFFSPRING:-10]])
mate.saveHaplotypes('ancestral_%s.%s.csv'%(POP2, CHR), subjects[1][idxPop2[NOFFSPRING:-12]], snpNames, 
                    snpPos,pops[1][:,idxPop2[NOFFSPRING:-12]])
mate.saveHaplotypes('admixed_%s_%s.%s.csv'%(POP1, POP2, CHR), ['ind']*NOFFSPRING, snpNames, snpPos, admixedPop)
mate.saveHaplotypes('admixed_origin_%s_%s.%s.csv'%(POP1, POP2, CHR), ['ind']*NOFFSPRING, snpNames, snpPos, admixedOrigin)

