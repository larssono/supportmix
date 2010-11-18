import simulateAdmixture as mate, numpy as np

DATA_DIR='../../../human_genome_data/data_hapmap3/'

POP_FILES=(DATA_DIR+'hapmap3_r2_b36_fwd.consensus.qc.poly.chr22_ceu.phased.gz',
           DATA_DIR+'data_hapmap3/hapmap3_r2_b36_fwd.consensus.qc.poly.chr22_yri.phased.gz',
           DATA_DIR+'human_genome_data/data_hapmap3/hapmap3_r2_b36_fwd.consensus.qc.poly.chr22_chd.unr.phased.gz')
F_GM=DATA_DIR+'genetic_map_chr22_b36.txt'
NGENS=8
NOFFSPRING=6

#Read populations
pops, nPops, subjects, nSNPs, snpPos, snpNames = mate.readFiles(POP_FILES)
idxPop1=np.random.permutation(nPops[0]) #Pick mating pairs
idxPop2=np.random.permutation(nPops[1])
idxPop3=np.random.permutation(nPops[2])
#Mix the first two populations
admixedPop, origin2=mate.poissonMating(pops[0][:,idxPop1[:NOFFSPRING]],
                                        pops[1][:,idxPop2[:NOFFSPRING]],
                                        snpPos, NGENS, .5)
# import pylab
# pylab.subplot(3,1,1)
# pylab.imshow(origin2.T, interpolation='nearest', cmap=pylab.cm.copper, vmin=0, vmax=3)
# pylab.ylabel('Sample ');pylab.yticks([]); pylab.xticks([]); pylab.axis('tight')

#Mix admixed with third population
admixedPop, admixedOrigin=mate.poissonMating(admixedPop,pops[2][:,idxPop3[:NOFFSPRING]],
                                        snpPos, NGENS, .7)

# pylab.subplot(3,1,2)
# pylab.imshow(admixedOrigin.copy().T, interpolation='nearest', cmap=pylab.cm.copper, vmin=0, vmax=2)
# pylab.ylabel('Sample ');pylab.yticks([]); pylab.xticks([]); pylab.axis('tight')

#Fix origin 
admixedOrigin[admixedOrigin==1]=2
idx=admixedOrigin==0
admixedOrigin[idx]=origin2[idx]

# pylab.subplot(3,1,3)
# pylab.imshow(admixedOrigin.T, interpolation='nearest', cmap=pylab.cm.copper, vmin=0, vmax=3)
# pylab.ylabel('Sample ');pylab.yticks([]); pylab.xticks([]); pylab.axis('tight')


#Save populations
mate.saveHaplotypes('ancestral_ceu.csv', subjects[0][idxPop1[NOFFSPRING:-10]], snpNames, 
                    snpPos, pops[0][:,idxPop1[NOFFSPRING:-10]])
mate.saveHaplotypes('ancestral_yri.csv', subjects[1][idxPop2[NOFFSPRING:-12]], snpNames, 
                    snpPos,pops[1][:,idxPop2[NOFFSPRING:-12]])
mate.saveHaplotypes('ancestral_chd.csv', subjects[2][idxPop3[NOFFSPRING:-10]], snpNames, 
                    snpPos,pops[2][:,idxPop3[NOFFSPRING:-10]])
mate.saveHaplotypes('admixed_ceu_yri_chd.csv', ['ind']*NOFFSPRING, snpNames, snpPos, admixedPop)
mate.saveHaplotypes('admixed_ceu_yri_chd_origin.csv', ['ind']*NOFFSPRING, snpNames, snpPos, admixedOrigin)

