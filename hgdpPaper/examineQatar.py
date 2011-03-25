###############################################
#Performs three analysis:
#  1. PCA of Qatar data and HGDP data together
#  3. Runs SupportMix on Qatar alone
###############################################
import sys; sys.path.append('../')
import fileReader, pylab, regionClassifier, time, glob, numpy as np
import string as str
from variables import *
from scipy.linalg import svd
from examineHGDP import classify

QATARFILES='data/CrystalQatar/phasedBeagle/qatarFlt.%(CHR)s.bgl.phased.gz'
HGDP_FILE='data/HGDP_raw_data/phasedBeagle/hgdp.%(CHR)s.bgl.phased.gz'
POPLABEL_FILE='data/HGDP_raw_data/samples_for_nstd27_5.csv'
QATARLABEL_FILE='data/CrystalQatar/haley_keepfile.txt'

#Create dictionary from name to population
popDict=dict([(l.split(',')[2], l.split(',')[7]) for l in open(POPLABEL_FILE)])
for ind, pop in (l.strip().split('\t')[1:] for l in open(QATARLABEL_FILE) ):
    popDict[ind]=pop
popDict.pop('Sample ID')

###############################################
#  Perform PCA
###############################################
#Read monlithic file containing all populations 
snpLabels=[]        #stores snp labels from in files
snpLocations=[]     #stores physical location from files
snpVals=[]
for i in range(1,23):
    CHR='chr%i' %i
    files=fileReader.concurrentFileReader(HGDP_FILE%locals(), QATARFILES%locals())
    subjects=files.next()
    for i, (snpName, snpLocation, snps) in enumerate(files):
        snpLabels.append(snpName)
        snpLocations.append(float(snpLocation))
        snpVals.append(fileReader.nucleotides2SNPs(sum(snps, [])))
subjects=np.hstack(subjects)
snpVals=np.asarray(snpVals)
popLabels=[popDict.get(s[:-2]) for s in subjects[::2]]
nSNPs, nSamples=snpVals.shape
#Normalize markers
snpVals=(snpVals-np.tile(snpVals.mean(1), (nSamples,1)).T)   #Mean center results
for i in range(nSNPs): snpVals[i,:]=snpVals[i,:]/np.sqrt(np.dot(snpVals[i,:],snpVals[i,:])) #Variance Scale

#Compute SVD and store results
U, S, Vt=svd(snpVals, 0)
S=S**2/np.sum(S**2)*100
np.savez(OUTPUT_PCA, popLabels=popLabels, Vt=Vt, S=S, subjects=subjects)  #Store results for plotting later


###################################################
#  Run SupportMix
###################################################
classifier=regionClassifier.SVMpymvpa(C)
NGENS=3
WINSIZE=100
C=1
print NGENS, WINSIZE, C

fpClass=open('data/qatarSupportMix/qatar.%(WINSIZE)i.admixedClass.csv'%locals(), 'w')
fpP=open('data/qatarSupportMix/qatar.%(WINSIZE)i.posterior.csv'%locals(), 'w')
fpPosition=open('data/qatarSupportMix/qatar.%(WINSIZE)i.position.csv'%locals(), 'w')
for i in range(1,23):
    t0=time.time()
    CHR='chr%i' %i
    admixedFile=QATARFILES%locals()
    fileNames=glob.glob('data/HGDP_raw_data/phasedBeagle/*/*%(CHR)s.bgl.phased.gz'%locals())
    fileNames.sort()
    smoother=regionClassifier.hmmFilter(geneticMapFile='data/hapmap2/genetic_map_%(CHR)s_b36.txt'%locals(),
                                               nGens=NGENS,nClasses=len(fileNames))
    pops=[file.split('/')[3] for file in fileNames]
    fileNames.append(admixedFile)
    admClassPre, admClass, p, subs, snpLocations, snpNames = classify(fileNames, smoother, classifier, WINSIZE)
    print '%s:%i:%i\t' %(CHR, (time.time()-t0)/60, (time.time()-t0)%60)
    #Save output
    for i in range(admClass.shape[0]): 
        fpClass.write('\t'.join(np.asarray(admClass[i,:], np.str))+'\n')
        fpP.write('\t'.join(np.asarray(p[i,:], np.str))+'\n')
        fpPosition.write('%s\t%i\t%i\t%s' %(CHR, snpLocations[i][0], snpLocations[i][-1], ','.join(snpNames[i]))+'\n')
np.save('data/qatarSupportMix/qatar.%(WINSIZE)i.populations'%locals(), pops)
np.save('data/qatarSupportMix/qatar.%(WINSIZE)i.subjects'%locals(), subs)
fpClass.close()
fpP.close()
fpPosition.close()
