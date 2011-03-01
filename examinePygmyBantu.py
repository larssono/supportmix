from runSVM import *
import mvpa.suite as pymvpa
from scipy.linalg import svd

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




###############################################################################################
## Main
###############################################################################################
winSize=200; nAncestral=60; nGens=100; winStep=200
classifier=pymvpa.LinearCSVMC(C=10)
cvte=pymvpa.CrossValidatedTransferError(pymvpa.TransferError(classifier),
                                        splitter=pymvpa.OddEvenSplitter(npertarget='equal'))
AllAncestralSuccess=[]
AllAdmixedClass=[]

for CHROM in range(3,4):
    fileNames=[LEMANDE%CHROM, NGUMBA%CHROM, TIKAR_S%CHROM, BAKOLA%CHROM, BAKA%CHROM, BEDZAN%CHROM]
    # fileNames=['tmp/pygmybantu_experimentation_20101014/ancestral_yoruba.chr%i.csv' %CHROM,
    #           'tmp/pygmybantu_experimentation_20101014/admixed_yoruba_karitiana.chr%i.csv.gz' %CHROM]    
    subjects, labels, snpNames, snpLocations, vals=readFiles(fileNames)
    origLabels=labels.copy()
    idxBantu=pylab.find((labels==0)+(labels==1)+(labels==2))
    idxPygmy=pylab.find((labels==3)+(labels==4)+(labels==5))
    #idxBantu=pylab.find(labels==0)
    #idxPygmy=pylab.find(labels==1)
    labels[idxBantu]=0
    labels[idxPygmy]=1
    nSubs, nSNPs=vals.shape
    startPos=0; ancestralSuccess=[]; admixedClass=[]; ancestralSamples=np.ones((nAncestral, vals.shape[1]))
    while startPos < nSNPs:
        winVals=vals[:,startPos:startPos+winSize]
        idx=(np.abs(winVals.sum(0))!=nSubs)     #Filter any SNPs without information across samples
        winVals=winVals[:, idx]

        ##Select Ancestral Pygmy based on cross validation
        #idxAncestralPygmy, idxAdmixedPygmy = selctAncestralCorrelation(winVals, idxBantu, idxPygmy, nAncestral)
        #idxAncestralPygmy, idxAdmixedPygmy = selctAncestralPCA(winVals, idxBantu, idxPygmy, nAncestral, -1)
        #idxAncestralPygmy, idxAdmixedPygmy = selectAncestralRan(winVals, idxBantu, idxPygmy, nAncestral)
        idxAncestralPygmy, idxAdmixedPygmy = selectAncestralL2Norm(winVals, idxBantu, idxPygmy, nAncestral)

        idxAncestral=list(idxBantu)+list(idxAncestralPygmy)
        ancestralSamples[:,startPos:startPos+winSize]=vals[idxAncestralPygmy,startPos:startPos+winSize]
        #Cross Validate ancestral Success rate
        ds=pymvpa.Dataset(winVals[idxAncestral,:], sa={'targets':labels[idxAncestral], 'chunks':range(len(idxAncestral))})
        cv_result=1-np.mean(cvte(ds))
        ancestralSuccess.append(cv_result)
        # print startPos, winVals.shape, cv_result

        ## Train on ancestral then classify admixed Pygmies
        classifier.train(ds)
        admixedClassWin=np.ones(len(idxPygmy))  #The ancestral "Pygmy" are classified as 1 i.e. Pygmy the rest are asigned below
        admixedClassWin[idxAdmixedPygmy-len(idxBantu)]=classifier.predict(winVals[idxAdmixedPygmy,:])

        #Play with the number of admixed on each side of plane
        # estimates=classifier.ca['estimates'].value

        admixedClass.append(admixedClassWin)
        startPos+=winStep
    admixedClassPre=np.asarray(admixedClass)
    #admixedClass=admixedClassPre
    smoother=regionClassifier.hmmFilter(geneticMapFile=MAPFILES%CHROM,nGens=nGens,nClasses=2)
    admixedClass, p=smoother(snpLocations, ancestralSuccess, admixedClassPre)  #TODO fix so this work

    print 'Percent Pygmy: %0.2g%% (%0.2g)' %((admixedClass==1).sum()/float(np.prod(admixedClass.shape))*100, float(nAncestral)/len(idxPygmy)*100)
    print 'Percent Success: %0.2g %%' %(np.mean(ancestralSuccess)*100)

    # sizeBP, sizeCM=winSizeBPandCM(snpLocations, winSize, CHROM)
    pylab.figure()
    pylab.axes([0.1, 0.4, .8, .5])
    #p[p<=np.sort(p.flatten())[-p.size*.30]]=0
    pylab.imshow(((admixedClass*2-1)*p).T, interpolation='nearest', cmap=pylab.cm.RdBu)
    print pylab.xticks()
    pylab.axis('tight'); xticks=range(pylab.xticks()[0][1], pylab.xticks()[0][-2], 10)# xticks=range(#xticks=np.asarray(pylab.xticks()[0][1:-1], np.int)
    pylab.xticks(xticks, np.asarray(snpLocations[::winStep][xticks], np.int))
    pylab.ylabel('Pygmy Samples')
    pylab.axes([0.1, 0.1, .8, .2])
    #pylab.plot(ancestralSuccess, color='blue' )
    pylab.plot(p.mean(1))
    pylab.plot(admixedClass.mean(1))
    pylab.ylim([.5, 1])
    pylab.ylabel('<Posterior p>')
    pylab.xlabel('Chromosome %i' %CHROM)
    pylab.xlim(0, p.shape[0])
    #pylab.xticks(xticks, np.asarray(snpLocations[::winStep][xticks], np.int))


# pylab.figure();plotPCA(np.vstack((vals, ancestralSamples)), np.hstack((labels, [4]*nAncestral)))
# pylab.show()

#####Simulated Results
# correctFile='tmp/pygmybantu_experimentation_20101014/admixed_origin_yoruba_karitiana.chr%i.csv.gz'%CHROM
# correct=np.array([l.split()[2:] for l in fileReader.openfile(correctFile).readlines()[1:]], np.float)
# svmClass=np.repeat(admixedClassPre, winSize, 0)
# hmmClass=np.repeat(admixedClass, winSize, 0)
# svmSuccess=100-sum(abs(svmClass[:len(correct),:]-correct))/len(correct)*100
# hmmSuccess=100-sum(abs(hmmClass[:len(correct),:]-correct))/len(correct)*100
# print 'Correct %0.3g +/- %0.2g (%0.3g +/- %0.2g)' %(np.mean(hmmSuccess), np.std(hmmSuccess), np.mean(svmSuccess), np.std(svmSuccess))
# pylab.figure();
# pylab.subplot(2,1,1);
# pylab.imshow(((admixedClass*2-1)*p).T, interpolation='nearest', cmap=pylab.cm.copper, vmin=0, vmax=2)
# pylab.ylabel('Sample '); pylab.yticks([]); pylab.xticks([]); pylab.axis('tight'); pylab.title('Estimat')
# pylab.subplot(2,1,2);
# pylab.imshow(correct[:, :].T, interpolation='nearest', cmap=pylab.cm.copper, vmin=0, vmax=2)
# pylab.ylabel('Sample');pylab.yticks([]); pylab.xticks([])
# pylab.xlabel('Position along %s' %CHROM);  pylab.axis('tight'); pylab.title('True ancestry')
# ###########################
