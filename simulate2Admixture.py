import simulateAdmixture as mate, numpy as np

DATA_DIR='../../human_genome_data/'
CHR='chr1'
POP1='yoruba'
POP2='karitiana'

POP_FILES=(DATA_DIR+'hapmap3/YRI/TRIOS/hapmap3_r2_b36_fwd.consensus.qc.poly.%s_%s.phased.gz'%(CHR,POP1), 
           DATA_DIR+'hapmap3/MKK/hapmap3_r2_b36_fwd.consensus.qc.poly.%s_%s.phased.gz'%(CHR, POP2))
POP_FILES=(DATA_DIR+'HumanPhasedData/hgdp/yoruba/hgdp_yoruba.%s.bgl.phased.gz'%CHR,
           DATA_DIR+'HumanPhasedData/hgdp/karitiana/hgdp_karitiana.%s.bgl.phased.gz'%CHR)
F_GM=DATA_DIR+'hapmap2/genetic_map_%s_b36.txt'%CHR
NGENS=20
NOFFSPRING=16

#Read populations
pops, nPops, subjects, nSNPs, snpPos, snpNames = mate.readFiles(POP_FILES)

idxPop1=np.random.permutation(nPops[0]) #Pick mating pairs
idxPop2=np.random.permutation(nPops[1])
#Mix the first two populations
admixedPop, admixedOrigin=mate.poissonMating(pops[0][:,idxPop1[:NOFFSPRING]],
                                             pops[1][:,idxPop2[:NOFFSPRING]],
                                             snpPos, F_GM, NGENS, .5)

# #Save populations
mate.saveHaplotypes('ancestral_%s.%s.csv'%(POP1, CHR), subjects[0][idxPop1[NOFFSPRING:]], snpNames, 
                    snpPos, pops[0][:,idxPop1[NOFFSPRING:]])
mate.saveHaplotypes('ancestral_%s.%s.csv'%(POP2, CHR), subjects[1][idxPop2[NOFFSPRING:]], snpNames, 
                    snpPos,pops[1][:,idxPop2[NOFFSPRING:]])
mate.saveHaplotypes('admixed_%s_%s.%s.csv'%(POP1, POP2, CHR), ['ind']*NOFFSPRING, snpNames, snpPos, admixedPop)
mate.saveHaplotypes('admixed_origin_%s_%s.%s.csv'%(POP1, POP2, CHR), ['ind']*NOFFSPRING, snpNames, snpPos, admixedOrigin)

