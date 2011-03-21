import sys; sys.path.append('../')
import fileReader, pylab, runSVM, time, glob, numpy as np
import string as str
from variables import *
from scipy.linalg import svd
from examineHGDP import classify
#Performs three analysis:
#  1. PCA of Qatar data and HGDP data together
#  3. Runs SupportMix on Qatar alone


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
#np.savez(OUTPUT_PCA, popLabels=popLabels, Vt=Vt, S=S, subjects=subjects)  #Store results for plotting later


###################################################
#  Run SupportMix
###################################################
classifier=runSVM.regionClassifier.SVMpymvpa(C)
NGENS=3
WINSIZE=100
C=1
print NGENS, WINSIZE, C
for i in range(1,23):
    t0=time.time()
    CHR='chr%i' %i
    admixedFile=QATARFILES%locals()
    fileNames=glob.glob('data/HGDP_raw_data/phasedBeagle/*/*%(CHR)s.bgl.phased.gz'%locals())
    fileNames.sort()
    smoother=runSVM.regionClassifier.hmmFilter(geneticMapFile='data/hapmap2/genetic_map_%(CHR)s_b36.txt'%locals(),
                                               nGens=NGENS,nClasses=len(fileNames))
    pops=[file.split('/')[3] for file in fileNames]
    fileNames.append(admixedFile)
    ancSuccess, admClassPre, admClass, p, subs = classify(fileNames, smoother, WINSIZE, classifier)
    np.save('data/qatarSupportMix/qatar.%(CHR)s.%(WINSIZE)i.admixedClass'%locals(), admClass)
    np.save('data/qatarSupportMix/qatar.%(CHR)s.%(WINSIZE)i.posterior'%locals(), p)
    np.save('data/qatarSupportMix/qatar.%(CHR)s.%(WINSIZE)i.populations'%locals(), pops)
    np.save('data/qatarSupportMix/qatar.%(CHR)s.%(WINSIZE)i.subjects'%locals(), subs)
    print '%s:%i:%i\t' %(CHR, (time.time()-t0)/60, (time.time()-t0)%60)

