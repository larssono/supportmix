import sys; sys.path.append('..')
import gzip

import regionClassifier
import fileReader
import mvpa.suite as pymvpa
from scipy.linalg import svd
import numpy as np
import pylab

CHROM=1
DATADIR='/home/lom/current_projects/human_genome_data/HumanPhasedData/pygmy_bantu_tishkoff/'
LEMANDE =DATADIR+'bantu_lemande.chr%i.bgl.phased.gz'    
NGUMBA  =DATADIR+'bantu_ngumba.chr%i.bgl.phased.gz'     
TIKAR_S =DATADIR+'bantu_tikar_south.chr%i.bgl.phased.gz'
BAKA    =DATADIR+'pygmy_baka.chr%i.bgl.phased.gz'       
BAKOLA  =DATADIR+'pygmy_bakola.chr%i.bgl.phased.gz'     
BEDZAN  =DATADIR+'pygmy_bedzan.chr%i.bgl.phased.gz'     
MAPFILES=DATADIR+'../../hapmap2/genetic_map_chr%i_b36.txt'

def readFiles(fileNames):    
    snpNames=[]
    snpLocations=[]  #stores physical location from files
    vals=[]          #Stores Values of genotypes
    files=fileReader.concurrentFileReader(*fileNames, key=0)
    subjects=files.next()
    labels=np.asarray(sum([[i]*len(sub) for i, sub in enumerate(subjects)], []))
    for i, (snpName, snpLocation, snps) in enumerate(files):
        snpLocations.append(float(snpLocation))
        snpNames.append(snpName)
        vals.append(fileReader.nucleotides2Haplotypes(sum(snps, [])))
    vals=np.asarray(vals).T
    snpLocations=np.asarray(snpLocations)
    return subjects, labels, snpNames, snpLocations, vals

    
def plotPCA(vals, labels):
    """Calculate PCA and plot top PCs."""
    [u,s,vt]=svd(vals.T,0)
    pylab.subplot(2,2,1);pylab.scatter(vt[0,:], vt[1,:],20, labels, edgecolors='none', vmax=5);pylab.xticks([]);pylab.yticks([]);pylab.xlabel('PC1');pylab.ylabel('PC2'); #[pylab.text(vt[0,i], vt[1,i], str(i+1), fontsize=7) for i in range(116, 262)]
    pylab.subplot(2,2,2);pylab.scatter(vt[2,:], vt[1,:],20, labels, edgecolors='none', vmax=5);pylab.xticks([]);pylab.yticks([]);pylab.xlabel('PC3');pylab.ylabel('PC2'); #[pylab.text(vt[1,i], vt[2,i], str(i+1), fontsize=7) for i in range(116, 262)]
    pylab.subplot(2,2,3);pylab.scatter(vt[4,:], vt[3,:],20, labels, edgecolors='none', vmax=5);pylab.xticks([]);pylab.yticks([]);pylab.xlabel('PC5');pylab.ylabel('PC4'); #[pylab.text(vt[2,i], vt[3,i], str(i+1), fontsize=7) for i in range(116, 262)]
    pylab.subplot(2,2,4);pylab.scatter(vt[2,:], vt[3,:],20, labels, edgecolors='none', vmax=5);pylab.xticks([]);pylab.yticks([]);pylab.xlabel('PC3');pylab.ylabel('PC4'); #[pylab.text(vt[3,i], vt[4,i], str(i+1), fontsize=7) for i in range(116, 262)]
    pylab.subplots_adjust(right=.8)
    cax = pylab.axes([0.85, 0.1, 0.02, 0.8])
    pylab.colorbar(cax=cax)
    pylab.yticks([0,.2,.4,.6,.8,1], ['Lemande', 'Ngumba', 'Tikar S', 'Baka', 'Bakola', 'Bedzan'])
    pylab.suptitle('chromosome %i' %CHROM)


def selctAncestralCorrelation(winVals, idxBantu, idxPygmy,nAncestral=4):
    """Out of the nXm matrix winVals selects those labeled 2 that are
    most different by median of correlation  from those labeled 1 """
    c=np.corrcoef(winVals)
    bantuXpygmy=c[idxBantu,:][:,idxPygmy]  #Average correlation with Bantu
    idx=np.argsort(np.median(bantuXpygmy, 0))  #Average correlation with Bantu
    # labels=np.zeros(len(idxBantu)+len(idxPygmy))
    # labels[idxPygmy[idx[:nAncestral]]]=2
    # labels[idxPygmy[idx[nAncestral:]]]=5
    # pylab.figure(); plotPCA(winVals, labels)
    return idxPygmy[idx[:nAncestral]], idxPygmy[idx[nAncestral:]]

def selectAncestralL2Norm(winVals, idxBantu, idxPygmy, nAncestral):
    """Runs classification on whole dataset picks Pygmies with largest projection from classification plane"""
    labels=np.zeros(winVals.shape[0]); labels[idxPygmy]=1    
    ds=pymvpa.Dataset(winVals, sa={'targets':labels})
    classifier.train(ds)
    classifier.predict(ds)
    idx=np.argsort(-classifier.ca['estimates'].value[idxPygmy]) # Is this the L2Norm?
    #use classifier.model.get_sv(), classifier.model.get_sv_coef(), classifier.model.get_rho() to determine furthest sampels
    return idxPygmy[idx[:nAncestral]], idxPygmy[idx[nAncestral:]]

def selectAncestralRan(winVals, idxBantu, idxPygmy,nAncestral=4):
    idx=np.random.permutation(len(idxPygmy))
    return idxPygmy[idx[:nAncestral]], idxPygmy[idx[nAncestral:]]
 
def selctAncestralPCA(winVals, idxBantu, idxPygmy,nAncestral=4, pc=-1):
    """Using the top two PC from a PCA compute those Ancestral samples
    by distance from center of mass of other ancestral population most.
    Parameters:
        pc specifies the PC to use:  by default a combination of 1 and 2 (-1)
    """
    [u,s,vt]=svd(winVals.T, 0)
    if pc==-1:
        vb=vt[:2, idxBantu].mean(1)
        vp=vt[:2, idxPygmy].mean(1)
        vec=(vp-vb)/np.sqrt(np.sum((vp-vb)**2))
        cosTheta=(vp-vb)[0]/ np.sqrt(np.sum((vp-vb)**2))
        x=cosTheta*np.sum(vt[:2,idxPygmy]**2, 0)**.5
        idx=np.argsort(x)
    else:
        if vt[pc, idxBantu].mean>vt[pc, idxPygmy].mean():
            idx=np.argsort(vt[pc,idxPygmy])
        else:
            idx=np.argsort(-vt[pc,idxPygmy])
    return idxPygmy[idx[:nAncestral]], idxPygmy[idx[nAncestral:]]

def saveTPedFile(fileName, winSize, chroms, subjects, allSNPNames, allSNPLocations, allAdmixedClass,allP):
    """Given ancestry assignments for multiple chromosomes save a tped and tfam file
    Arguments:
    - `fileName`: Base name of filenames
    - `chroms`:   list of chromsome numbers
    - `subjects`: list of list of subjects
    - `snpNamesLocations`: list of all snp names in each chromsome
    - `snpLocations`: list of all snp locations in each chromsome
    - `ddmixedClass`: list of classifications, one item per chromsome
    - `p`: list of posterior probability for every window
    """
    nSubs=134#len(subjects)
    with open(fileName+'.tfam', 'w') as pedFile:
        for i, sub in enumerate(subjects[::2]):
            pedFile.write('%s %s 0 0 0 0\n' %(sub[:-2], sub[:-2]))
    tPedFp=gzip.open(fileName+".tped.gz","w")
    tPedProbFp=gzip.open(fileName+".Probs.tped.gz","w")
    for i, chrom in enumerate(chroms):
        snpNames=allSNPNames[i]
        snpLocations=allSNPLocations[i]
        admixedClass=allAdmixedClass[i]
        p=allP[i]
        for j, (rsId, rsPos) in enumerate(zip(snpNames, snpLocations)):
            tPedFp.write('%s %s %i %i' %(chrom, rsId, j/winSize, rsPos))
            tPedFp.write(' %i'*nSubs %tuple(admixedClass[j/winSize, :]+1))
            tPedFp.write('\n')
            tPedProbFp.write('%s %s %i %i' %(chrom, rsId, j/winSize, rsPos))
            tPedProbFp.write(' %0.3g'*nSubs %tuple(p[j/winSize, :]))
            tPedProbFp.write('\n')
    tPedFp.close()
    tPedProbFp.close()


###############################################################################################
## Main
###############################################################################################
chroms=range(1,23)
winSize=50; nAncestral=72; nGens=20; winStep=50
fileNames=[LEMANDE, NGUMBA, TIKAR_S, BAKOLA, BAKA, BEDZAN]
alphaFile='pygmy_bantu_alphas.txt'
# winSize=100; nAncestral=30; nGens=8; winStep=100
# fileNames=['tmp/ancestral_ceu.chr%i.csv','tmp/admixed_ceu_yri.chr%i.csv']    
# alphaFile='tmp/admixed_ceu_yri.alpha.chr22.csv'
# correctFile='tmp/admixed_ceu_yri_origin.chr22.csv'

classifier=pymvpa.LinearCSVMC(C=10)
cvte=pymvpa.CrossValidation(classifier, pymvpa.NGroupPartitioner(3, attr='runtype'))
alphas=[1-float(f.strip().split()[-1]) for f in open(alphaFile)]
allAdmixedClass=[]
allAdmixedP=[]
allSNPLocations=[]
allSNPNames=[]

#Run the ancestry assignemts
for CHROM in chroms:
    subjects, labels, snpNames, snpLocations, vals=readFiles([f%CHROM for f in fileNames])
    origLabels=labels.copy()
    if len(fileNames)==6:  #Real data
        idxBantu=pylab.find((labels==0)+(labels==1)+(labels==2))
        idxPygmy=pylab.find((labels==3)+(labels==4)+(labels==5))
    else: #Simulated data
        idxBantu=pylab.find(labels==0)
        idxPygmy=pylab.find(labels==1)
    labels[idxBantu]=0
    labels[idxPygmy]=1
    nSubs, nSNPs=vals.shape
    startPos=0; ancestralSuccess=[]; admixedClass=[]; ancestralSamples=np.ones((nAncestral, vals.shape[1]))
    while startPos < nSNPs:
        winVals=vals[:,startPos:startPos+winSize]
        idx=(np.abs(winVals.sum(0))!=nSubs)     #Filter any SNPs without information across samples
        winVals=winVals[:, idx]
        if winVals.shape[1]>3:  #Enough information to determine results
            ##Select Ancestral my based on...
            idxAncestralPygmy, idxAdmixedPygmy = selectAncestralL2Norm(winVals, idxBantu, idxPygmy, nAncestral)
            idxAncestral=list(idxBantu)+list(idxAncestralPygmy)
            ancestralSamples[:,startPos:startPos+winSize]=vals[idxAncestralPygmy,startPos:startPos+winSize]
            #Cross Validate ancestral Success rate
            runtype=np.zeros(len(idxAncestral)); runtype[0::3]=0;runtype[1::3]=1; runtype[2::3]=2 
            ds=pymvpa.Dataset(winVals[idxAncestral,:], sa={'targets':labels[idxAncestral], 'runtype':runtype, 'chunks':range(len(idxAncestral))})
            cv_result=1-np.mean(cvte(ds))
            ancestralSuccess.append(cv_result)
            ## Train on ancestral then classify admixed Pygmies
            classifier.train(ds)
            admixedClassWin=np.ones(len(idxPygmy))  #The ancestral "Pygmy" are classified as 1 i.e. Pygmy the rest are asigned below
            admixedClassWin[idxAdmixedPygmy-len(idxBantu)]=classifier.predict(winVals[idxAdmixedPygmy,:])
        else:  #Window is too small to determine anything just assign ancestry to 1
            admixedClassWin=np.ones(len(idxPygmy))
            ancestralSuccess.append(0.5)

        admixedClass.append(admixedClassWin)
        startPos+=winStep
    admixedClassPre=np.asarray(admixedClass)
    smoother=regionClassifier.hmmFilter(geneticMapFile=MAPFILES%CHROM,nGens=nGens,nClasses=2)
    admixedClass, p=smoother(snpLocations, ancestralSuccess, admixedClassPre)  
    #Store all results
    allAdmixedClass.append(admixedClass)
    allAdmixedP.append(p)
    allSNPLocations.append(snpLocations)
    allSNPNames.append(snpNames)
    #Print summary of results
    print 'Chromosome %i: Percent Pygmy: %0.2g%% (%0.2g) ' %(CHROM,(admixedClass==1).sum()/float(np.prod(admixedClass.shape))*100, float(nAncestral)/len(idxPygmy)*100)

#Plot results of multiple Chromsomes
pylab.figure(figsize=(19.5,5))
nPygmy=admixedClass.shape[1]
nWins=[x.shape[0] for x in allAdmixedClass]
xStart=0.03
for i, CHROM in enumerate(chroms):
    admixedClass=allAdmixedClass[i]
    p=allAdmixedP[i]
    nWin=nWins[i]
    snpLocations=allSNPLocations[i]
    xWidth=.84*(nWin/float(np.sum(nWins)))
    pylab.axes([xStart, 0.29, xWidth, .52])
    pylab.imshow(((admixedClass*2-1)*p).T, interpolation='nearest', cmap=pylab.cm.RdBu)
    pylab.axis('tight');     pylab.title('Chr %i' %CHROM, fontsize=8)
    xtickPos=np.asarray([2, round(nWin/2), nWin-2], np.int)
    pylab.axis('off'); 
    if xStart==0.03:
        pylab.axis('on'); pylab.xticks([]); 
        pylab.yticks([0, 24,58,83,108,121,134], ['|', 'Bakola', '|','Baka','|','Bedzan', '|'], fontsize=8, rotation=90)
    pylab.axes([xStart, 0.1, xWidth, .18])  #Average Pygmy Ancestry
    pylab.plot(admixedClass.mean(1))
    pylab.ylim([0.5, 1]); pylab.xlim(0, p.shape[0]); pylab.yticks([])
    pylab.xticks(xtickPos, np.asarray(np.round(snpLocations[::winStep][xtickPos]/1e6), np.int), fontsize=6)
    if xStart==0.03:
        pylab.yticks([.5,1], fontsize=8); pylab.ylabel('Percent Pygmy', fontsize=8);
    if xStart>0.45 and xStart<0.5:
        pylab.xlabel('Position [Mb]', fontsize=8)
    pylab.draw()
    xStart+=(xWidth+0.002)
#Plot average Pygmy Ancestry per sample
pylab.axes([0.92, 0.29, .07, .52]) 
pylab.plot(np.vstack(allAdmixedClass).mean(0), range(nPygmy)); pylab.ylim(nPygmy, 0); pylab.yticks([]); pylab.title('Percent Pygmy'); pylab.xlim(0,1); pylab.xticks([0,.5,1], fontsize=8)
pylab.plot(alphas, range(nPygmy)); pylab.ylim(nPygmy, 0); pylab.yticks([]); pylab.title('Percent Pygmy', fontsize=8); pylab.xlim(0,1); pylab.xticks([0,.5,1])
#Colorbar
pylab.axes([0.925, 0.1, .06, .18]);  pylab.xlim(-1, 1); pylab.ylim(-1, 1); pylab.axis('off')
pylab.colorbar(fraction=.5, aspect=4, orientation='horizontal', ticks=[-.99, 0, .99])
pylab.text(-1,-1.5, 'Bantu', fontsize=8)
pylab.text(1,-1.5, 'Pygmy', fontsize=8, horizontalalignment='right')

#Save output in tped plink files
saveTPedFile('pygmy_admixture_72b', winSize, chroms, np.hstack(subjects[3:]), allSNPNames, allSNPLocations, allAdmixedClass,allAdmixedP)


# If simulated data calculate Success Rate
if locals().get('correctFile'):
    correct=np.array([l.split()[2:] for l in fileReader.openfile(correctFile).readlines()[1:]], np.float)
    svmClass=np.repeat(admixedClassPre, winSize, 0)
    hmmClass=np.repeat(admixedClass, winSize, 0)
    svmSuccess=100-sum(abs(svmClass[:len(correct),:]-correct))/len(correct)*100
    hmmSuccess=100-sum(abs(hmmClass[:len(correct),:]-correct))/len(correct)*100
    print 'Correct %0.3g +/- %0.2g (%0.3g +/- %0.2g)' %(np.mean(hmmSuccess), np.std(hmmSuccess), np.mean(svmSuccess), np.std(svmSuccess))
    pylab.figure();
    pylab.subplot(2,1,1);
    pylab.imshow(((admixedClass*2-1)*p).T, interpolation='nearest', cmap=pylab.cm.copper, vmin=0, vmax=2)
    pylab.ylabel('Sample '); pylab.yticks([]); pylab.xticks([]); pylab.axis('tight'); pylab.title('Estimat')
    pylab.subplot(2,1,2);
    pylab.imshow(correct[:, :].T, interpolation='nearest', cmap=pylab.cm.copper, vmin=0, vmax=2)
    pylab.ylabel('Sample');pylab.yticks([]); pylab.xticks([])
    pylab.xlabel('Position along %s' %CHROM);  pylab.axis('tight'); pylab.title('True ancestry')

pylab.show()
