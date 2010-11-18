import PyML, pylab
import fileReader, sys, gzip
import numpy as np
from  hmm import hmm

WIN_SIZE=100
SVM=PyML.SVM(C=10000, optimizer='mysmo')
nGens=1
DOSNPS=False

DOTHREE=False
GENETIC_MAP_FILE='data_simulated/genetic_map_chr22_b36.txt'

ANC_1_FILE='data_simulated/ancestral1.csv'
ANC_2_FILE='data_simulated/ancestral2.csv'
#ANC_2_FILE='data_hapmap3/hapmap3_r2_b36_fwd.consensus.qc.poly.chr22_mkk.phased.gz'
#ANC_2_FILE='data_hapmap3/hapmap3_r2_b36_fwd.consensus.qc.poly.chr22_lwk.unr.phased.gz'
ADMIXED_FILE='data_simulated/admixed.csv'
ADMIXED_CORRECT_FILE='data_simulated/admixed_origin.csv'

# ANC_1_FILE='data_simulated/three/ancestral1.csv'
# ANC_2_FILE='data_simulated/three/ancestral2.csv'
# ANC_3_FILE='data_simulated/three/ancestral2.csv'
# ADMIXED_FILE='data_simulated/three/admixed.csv'
# ADMIXED_CORRECT_FILE='data_simulated/three/admixed_origin.csv'

#------------------------------------------------------------
# hmm stuff
#------------------------------------------------------------
class geneticMap(object):
    """keeps track of genetic Map locations and returns closests genetic map location 
    given a snp location. """

    def __init__(self,file ):
        """ """
        fp=open(file)
        self.m=np.asarray([np.asarray(l.split())[[0,2]] for l in fp.readlines()[1:]], np.float)

    def pos2gm(self, pos):
        """Converts position in bp to position in centiMorgans"""
        m=self.m
        i=m[:,0].searchsorted(pos)
        try:
            if m[i,0]==pos or i==0:
                return m[i,1]
            elif i==0:
                return m[0,1]
            else:  #linear interpolation
                return (m[i,1]-m[i-1,1])/(m[i,0]-m[i-1,0])*(pos-m[i-1,0]) + m[i-1,1]
        except IndexError:
            if i==len(m):
                return m[-1,1]
            else:
                raise IndexError

def runHMM(snpLocations, successRate, admixedClass):
    """Smoothes out admixedClass using HMM"""
    gm=geneticMap(GENETIC_MAP_FILE)
    mapLocations=map(gm.pos2gm, snpLocations)
    #determine transition matrices
    a=[]; b=[]
    oldPos=0
    for i in range(0,len(mapLocations), WIN_SIZE):
        newPos=np.mean(mapLocations[i:i+WIN_SIZE])
        dM=-(newPos - oldPos)/100*nGens
        e=np.exp(dM)  #TODO this will not work for snps with three classes
        oldPos=newPos
        a.append([[e, 1-e],[1-e, e]])#TODO this will not work for snps with three classes
    a=np.asarray(a)
    for s in successRate:
        b.append([[s, 1-s],[1-s, s]])#TODO this will not work for snps with three classes
    b=np.asarray(b)
    model=hmm(a, b)

    results=[]
    for i in range(admixedClass.shape[1]):
        model.forward_backward(admixedClass[:,i])
        maxIdx=model.pstate.argsort(1)[:,-1]
        results.append(maxIdx)
    return np.array(results).T

#------------------------------------------------------------
# Classifier stuff
#------------------------------------------------------------
def createSVMData(vals, labels):
    """Given array and labels creates PyML with normalization and  Kernel"""
    dataPyML=PyML.VectorDataSet(vals, L=labels)
    dataPyML.normalize(1)
    dataPyML.attachKernel('polynomial', degree = 1)
    return dataPyML

def stepSVM(valsTrain, labelsTrain, valsTest):
    """Does cv on ancestral population then trains on ancestral
    population followed by testing on admixed population.
    """
    #test ancestral populations

    if DOTHREE:
        import PyML.classifiers
        classifier=PyML.classifiers.multi.OneAgainstRest(SVM)
    else:
        classifier=SVM
    haplotypeData=createSVMData(valsTrain, labelsTrain)
    results=classifier.cv(haplotypeData, 3);
    ancestralSuccess = results.getBalancedSuccessRate()
    #train on ancestral populations
    classifier.train(haplotypeData);
    #classify admixed population
    testData=createSVMData(valsTest,[1]*valsTest.shape[0])
    admixedClass = [classifier.classify(testData,i)[0] for i in range(valsTest.shape[0])]
    return ancestralSuccess, admixedClass

def readCorrect():
    """Read the file of correct classifications and returns
    the classes or the genotyped class.

    For the latter the conversions are made using the values in
    adjacent columns and converting: 0,0->0, 0,1->1, 1,1->2, 1,2->3,
    2,0->4, 2,2->5 etc."""
    fileData=np.array([l.split()[2:] for l in fileReader.openfile(ADMIXED_CORRECT_FILE).readlines()[1:]], np.float)
    #if DOSNPS:
    #    nGroups=fileData.max()
    # correct=fileReader.nucluotides2SNPs(correct) if DOSNPS else fileReader.nucluotides2Haplotypes(correct)
    return fileData


if __name__ == '__main__':
    if DOTHREE:
        files=fileReader.concurrentFileReader(ANC_1_FILE, ANC_2_FILE, ANC_3_FILE, ADMIXED_FILE)
    else:
        files=fileReader.concurrentFileReader(ANC_1_FILE, ANC_2_FILE, ADMIXED_FILE)
    subjects=files.next()
    nTrain=np.sum(map(len, subjects[:-1]))  #Number of samples in training set
    nTest=len(subjects[-1]);
    labelsTrain =sum([[i]*len(sub) for i, sub in enumerate(subjects[:-1])],[])   #Arbitrary labels for ancestral populations
    if DOSNPS:
        nTrain/=2; nTest/=2
        labelsTrain=labelsTrain[::2]

    snpLocations=[]     #stores physical location from files
    ancestralSuccess=[] #stores success of ancestral classification
    admixedClass=[]     #stores classification of test Subjects
    vals=np.zeros((nTrain+nTest, WIN_SIZE))  #temporary storage of output
    while True: #for j in range(100): #To go through all positions in file
        for i, (snpName, snpLocation, snps) in enumerate(files):
            snpLocations.append(float(snpLocation))
            vals[:,i] = fileReader.nucleotides2SNPs(sum(snps, [])) if DOSNPS else fileReader.nucleotides2Haplotypes(sum(snps, []))
            if i==vals.shape[1]-1:
                break
        ancestral, admixed=stepSVM(vals[:nTrain,:i], labelsTrain, vals[-nTest:, :i])
        ancestralSuccess.append(ancestral)
        admixedClass.append(admixed)
        if i<WIN_SIZE-1:
            break
    admixedClassPre=np.array(admixedClass)
    admixedClass=runHMM(snpLocations, ancestralSuccess, admixedClassPre)

    #Compare and find successRate
    correct=readCorrect()
    svmClass=np.repeat(admixedClassPre, WIN_SIZE, 0)
    hmmClass=np.repeat(admixedClass, WIN_SIZE, 0)

    svmSuccess=100-sum(abs(svmClass[:len(correct),:]-correct))/len(correct)*100
    hmmSuccess=100-sum(abs(hmmClass[:len(correct),:]-correct))/len(correct)*100
    print np.mean(svmSuccess),'+/-', np.std(svmSuccess), np.mean(hmmSuccess),'+/-',np.std(hmmSuccess)
    print np.mean(svmSuccess),'+/-', np.std(svmSuccess)

    pylab.figure()
    pylab.clf()
    pylab.subplot(4,1,1);
    pylab.errorbar(range(len(ancestralSuccess)), map(np.mean, ancestralSuccess), map(np.std, ancestralSuccess))
    pylab.xlim([0, len(ancestralSuccess)]); pylab.xticks([])
    pylab.ylim([0,1]); pylab.ylabel('Success rate')
    pylab.title('Ancestral Population')


    pylab.subplot(4,1,2);
    pylab.imshow(admixedClassPre.T, interpolation='nearest', cmap=pylab.cm.copper, vmin=0, vmax=2)
    pylab.ylabel('Sample '); pylab.yticks([]); pylab.xticks([])
    pylab.axis('tight')

    pylab.subplot(4,1,3);
    pylab.imshow(admixedClass.T, interpolation='nearest', cmap=pylab.cm.copper, vmin=0, vmax=2)
    pylab.ylabel('Sample ');pylab.yticks([]); pylab.xticks([])
    pylab.axis('tight')

    pylab.subplot(4,1,4);
    pylab.imshow(correct[:, :].T, interpolation='nearest', cmap=pylab.cm.copper, vmin=0, vmax=2)
    pylab.ylabel('Sample ');pylab.yticks([]); pylab.xticks([])
    pylab.xlabel('Position along Chromosome 22')
    pylab.axis('tight')
    pylab.show()
    
