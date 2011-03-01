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
    pylab.subplot(2,2,2);pylab.scatter(vt[1,:], vt[2,:],20, labels, edgecolors='none', vmax=5);pylab.xticks([]);pylab.yticks([]);pylab.xlabel('PC2');pylab.ylabel('PC3'); #[pylab.text(vt[1,i], vt[2,i], str(i+1), fontsize=7) for i in range(116, 262)]
    pylab.subplot(2,2,3);pylab.scatter(vt[2,:], vt[3,:],20, labels, edgecolors='none', vmax=5);pylab.xticks([]);pylab.yticks([]);pylab.xlabel('PC3');pylab.ylabel('PC4'); #[pylab.text(vt[2,i], vt[3,i], str(i+1), fontsize=7) for i in range(116, 262)]
    pylab.subplot(2,2,4);pylab.scatter(vt[3,:], vt[4,:],20, labels, edgecolors='none', vmax=5);pylab.xticks([]);pylab.yticks([]);pylab.xlabel('PC4');pylab.ylabel('PC5'); #[pylab.text(vt[3,i], vt[4,i], str(i+1), fontsize=7) for i in range(116, 262)]
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
    return idxPygmy[idx[:nAncestral]], idxPygmy[idx[nAncestral:]]

def selectAncestralL2Norm(winVals, idxBantu, idxPygmy, nAncestral):
    """Runs classification on whole dataset picks Pygmies with largest projection from classification plane"""
    ds=pymvpa.Dataset(winVals)
    ds.sa['targets']=labels
    ds.sa['chunks']=range(len(labels))
    classifier.train(ds)
    #use classifier.model.get_sv(), classifier.model.get_sv_coef(), classifier.model.get_rho() to determine furthest sampels
    return None

def selectAncestralExhaustive(winVals, idxmBantu, idxPygmy):
    """Exhaustively attempts to find 4 ancestral Pygmies based on cross validation
    of ancestral population.  !EXTREMELY INEFICIENT
    """
    #TODO!  fix this to be flexible
    cvte=pymvpa.CrossValidatedTransferError(pymvpa.TransferError(classifier),
                                            splitter=pymvpa.CustomSplitter([(None, [1,119]),(None, [2,118]),(None, [3,117]),(None, [4,116]) ]))
    for i in range(0, len(idxPygmy)):
        for j in range(i+1, len(idxPygmy),2):
            for k in range(j+1, len(idxPygmy),2):
                for l in range(k+1, len(idxPygmy),2):
                    idxAncestral=list(idxBantu) + list(idxPygmy[[i,j,k,l]])
                    ds=pymvpa.Dataset(winVals[idxAncestral,:])
                    ds.sa['targets']=labels[idxAncestral]
                    ds.sa['chunks']=range(len(idxAncestral))
                    cv_result=1-np.mean(cvte(ds))
                    if cv_result>0.75:
                        print idxPygmy[[i,j,k,l]], cv_result
    #TODO return when finnished

###############################################################################################
## Main
###############################################################################################
winSize=200; nAncestral=20; nGens=100; winStep=200
classifier=pymvpa.LinearCSVMC(C=10)
cvte=pymvpa.CrossValidatedTransferError(pymvpa.TransferError(classifier),
                                        splitter=pymvpa.OddEvenSplitter(npertarget='equal'))
AllAncestralSuccess=[]
AllAdmixedClass=[]

for CHROM in range(1,23):
    fileNames=[LEMANDE%CHROM, NGUMBA%CHROM, TIKAR_S%CHROM, BAKOLA%CHROM, BAKA%CHROM, BEDZAN%CHROM]
    subjects, labels, snpNames, snpLocations, vals=readFiles(fileNames)
    idxBantu=pylab.find((labels==0)+(labels==1)+(labels==2))
    idxPygmy=pylab.find((labels==3)+(labels==4)+(labels==5))
    labels[idxBantu]=0
    labels[idxPygmy]=1
    nSubs, nSNPs=vals.shape
    startPos=0; ancestralSuccess=[]; admixedClass=[]; ancestralSamples=np.ones((nAncestral, vals.shape[1]))
    while startPos < nSNPs:
        winVals=vals[:,startPos:startPos+winSize]
        idx=(np.abs(winVals.sum(0))!=nSubs)     #Filter any SNPs without information across samples
        winVals=winVals[:, idx]

        ##Select Ancestral Pygmy based on cross validation
        idxAncestralPygmy, idxAdmixedPygmy = selctAncestralCorrelation(winVals, idxBantu, idxPygmy, nAncestral)
        idxAncestral=list(idxBantu)+list(idxAncestralPygmy)
        ancestralSamples[:,startPos:startPos+winSize]=vals[idxAncestralPygmy,startPos:startPos+winSize]
        #Cross Validate ancestral Success rate
        ds=pymvpa.Dataset(winVals[idxAncestral,:])
        ds.sa['targets']=labels[idxAncestral]
        ds.sa['chunks']=range(len(idxAncestral))
        cv_result=1-np.mean(cvte(ds))
        ancestralSuccess.append(cv_result)
        print startPos, winVals.shape, cv_result

        ## Train on ancestral then classify admixed Pygmies
        classifier.train(ds)
        admixedClassWin=np.ones(len(idxPygmy))  #The ancestral "Pygmy" are classified as 1 i.e. Pygmy the rest are asigned below
        admixedClassWin[idxAdmixedPygmy-len(idxBantu)]=classifier.predict(winVals[idxAdmixedPygmy,:])
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
    pylab.imshow(((admixedClass*2-1)*p).T, interpolation='nearest', cmap=pylab.cm.RdBu)
    pylab.axis('tight'); xticks=np.asarray(pylab.xticks()[0][1:-1], np.int)
    pylab.xticks(xticks, np.asarray(snpLocations[::winStep][xticks], np.int))
    pylab.ylabel('Pygmy Samples')
    pylab.axes([0.1, 0.1, .8, .2])
    pylab.plot(ancestralSuccess, color='blue' )
    pylab.plot(p.mean(1))
    pylab.ylim([.5, 1])
    pylab.ylabel('Classifacation Rate')
    pylab.xlabel('Chromosome %i' %CHROM)
    pylab.xlim(0, p.shape[0])
    pylab.xticks(xticks, np.asarray(snpLocations[::winStep][xticks], np.int))


#pylab.figure();plotPCA(np.vstack((vals, ancestralSamples)), np.hstack((labels, [4]*nAncestral)))
pylab.show()

