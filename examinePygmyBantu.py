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
winSize=3000; nAncestral=3; nGens=10
classifier=pymvpa.LinearCSVMC(C=10)
cvte=pymvpa.CrossValidatedTransferError(pymvpa.TransferError(classifier),
                                        splitter=pymvpa.OddEvenSplitter(npertarget='equal'))
AllAncestralSuccess=[]
AllAdmixedClass=[]

# ######## Simulated Data
# fileNames=['pygmybantu_experimentation/ancestral_yoruba.chr1.csv',
#            'pygmybantu_experimentation/admixed_yoruba_karitiana.chr1.csv']
# subjects, labels, snpNames, snpLocations, vals=readFiles(fileNames)
# idxBantu=pylab.find(labels==0)
# idxPygmy=pylab.find(labels==1)
# nSubs, nSNPs=vals.shape
# ############################


for pygmyPop in [BEDZAN, BAKOLA, BAKA]:
    admixedFp=open(pygmyPop.split('/')[-1].split('.')[0]+'.output.csv', 'w')
    for CHROM in range(1,23):
        print 'Chrom: %i Population: %s' %(CHROM, pygmyPop.split('/')[-1].split('.')[0])
        fileNames=[LEMANDE%CHROM, NGUMBA%CHROM, TIKAR_S%CHROM, pygmyPop%CHROM]
        subjects, labels, snpNames, snpLocations, vals=readFiles(fileNames)
        if  CHROM==1:
            admixedFp.write('Chrom\t Position\t AncestralSuccess\t')
            admixedFp.write('\t'.join(subjects[-1])+'\n')
        #plotPCA(vals, labels); pylab.show()
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
            # #Evaluate choice of ancestral Pygmies with PCA
            # [u,s,vt]=svd(winVals.T,0)
            # pylab.subplot(4,3,i+1)
            # pylab.scatter(vt[1,:], vt[2,:],20, labels, edgecolors='none', vmax=3);pylab.xticks([]);pylab.yticks([])
            # pylab.scatter(vt[1,idxAncestralPygmy], vt[2,idxAncestralPygmy],20, 'green', edgecolors='none', vmax=3);pylab.xticks([]);pylab.yticks([])
            # i+=1
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
            startPos+=winSize
        admixedClassPre=np.asarray(admixedClass)
        smoother=regionClassifier.hmmFilter(geneticMapFile=MAPFILES%CHROM,nGens=nGens,nClasses=2)
        #admixedClass=smoother(snpLocations, np.ones_like(ancestralSuccess)*np.mean(ancestralSuccess), admixedClassPre)  #TODO fix so this works
        admixedClass=smoother(snpLocations, ancestralSuccess, admixedClassPre)  #TODO fix so this works
        #admixedClass=admixedClassPre
        
        #Save output
        idx=np.asarray(ancestralSuccess)<=0.5
        admixedClass[idx, :]=-1
        j=-1
        for i in range(len(snpLocations)):
            if i%winSize==0: j+=1
            admixedFp.write('%s\t%s\t%s\t' %(CHROM, snpLocations[i], ancestralSuccess[j])) #             
            admixedFp.write('\t'.join(np.asarray(admixedClass[j,:], np.str_)))
            admixedFp.write('\n')
                       
        sizeBP, sizeCM=winSizeBPandCM(snpLocations, winSize, CHROM)
        # pylab.figure()
        # pylab.axes([0.1, 0.1, .8, .2])
        # pylab.plot(ancestralSuccess, color='blue' )
        # pylab.ylim([.5, 1])
        # pylab.ylabel('Success Rate besed Cross Valditation')
        # pylab.title('Chromosome %i' %CHROM)
        # pylab.axes([0.1, 0.4, .8, .5])
        # pylab.imshow(np.asarray(admixedClass).T, interpolation='nearest')
        # pylab.axis('tight')
        # pylab.figure();plotPCA(np.vstack((vals, ancestralSamples)), np.hstack((labels, [4]*nAncestral)))
        # pylab.show()
    admixedFp.close()


# #####Simulated Results
# correctFile='pygmybantu_experimentation/admixed_origin_yoruba_karitiana.chr1.csv'
# correct=np.array([l.split()[2:] for l in fileReader.openfile(correctFile).readlines()[1:]], np.float)
# #Compare and find successRate
# #admixedClass=admixedClassPre
# svmClass=np.repeat(admixedClassPre, winSize, 0)
# hmmClass=np.repeat(admixedClass, winSize, 0)
# svmSuccess=100-sum(abs(svmClass[:len(correct),:]-correct))/len(correct)*100
# hmmSuccess=100-sum(abs(hmmClass[:len(correct),:]-correct))/len(correct)*100
# print 'Correct %0.3g +/- %0.2g (%0.3g +/- %0.2g)' %(np.mean(hmmSuccess), np.std(hmmSuccess), 
#                                                     np.mean(svmSuccess), np.std(svmSuccess))
# pylab.figure();
# pylab.subplot(2,1,1);
# pylab.imshow(admixedClass.T, interpolation='nearest', 
#              cmap=pylab.cm.copper, vmin=0, vmax=2)
# pylab.ylabel('Sample '); pylab.yticks([]); pylab.xticks([]); pylab.axis('tight'); pylab.title('Estimated Ancestry')
# pylab.subplot(2,1,2);
# pylab.imshow(correct[:, :].T, interpolation='nearest', cmap=pylab.cm.copper, vmin=0, vmax=2)
# pylab.ylabel('Sample');pylab.yticks([]); pylab.xticks([])
# pylab.xlabel('Position along %s' %CHROM);  pylab.axis('tight'); pylab.title('True ancestry')
# ###########################

