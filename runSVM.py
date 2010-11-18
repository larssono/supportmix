import PyML, pylab
import numpy as np
import fileReader
import sys
from  hmm import hmm
import gzip

WIN_SIZE=100
SVM=PyML.SVM(C=10000, optimizer='mysmo')
nGens=1

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

def runHMM(mapLocations, successRate, admixedClass):
    """Smoothes out admixedClass using HMM"""

    #determine transition matrices
    a=[]; b=[]
    oldPos=0
    for i in range(0,len(mapLocations), WIN_SIZE):
        newPos=np.mean(mapLocations[i:i+WIN_SIZE])
        dM=-(newPos - oldPos)/100*nGens
        e=np.exp(dM)
        oldPos=newPos
        a.append([[e, 1-e],[1-e, e]])
    a=np.asarray(a)
    for s in successRate:
        b.append([[s, 1-s],[1-s, s]])
    b=np.asarray(b)
    model=hmm(a, b)

    results=[]
    for i in range(admixedClass.shape[1]):
        model.forward_backward(admixedClass[:,i])
        maxIdx=model.pstate.argsort(1)[:,-1]
        results.append(maxIdx)
    return np.array(results).T


def createSVMData(vals, labels):
    """Given array and labels creates PyML with normalization and  Kernel"""
    dataPyML=PyML.VectorDataSet(vals, L=labels)
    dataPyML.normalize(1)
    dataPyML.attachKernel('polynomial', degree = 1)
    return dataPyML

def stepSVM(valsTrain, labelsTrain, valsTest, labelsTest):
    """Does cv on ancestral population then trains on ancestral
    population followed by testing on admixed population.
    """
    #test ancestral populations
    haplotypeData=createSVMData(valsTrain, labelsTrain)
    results=SVM.cv(haplotypeData, 3);
    ancestralSuccess = results.getBalancedSuccessRate()
    #train on ancestral populations
    SVM.train(haplotypeData);
    #classify admixed population
    testData=createSVMData(valsTest, labelsTest)
    admixedClass = [SVM.classify(testData,i)[0] for i in range(valsTest.shape[0])]
    return ancestralSuccess, admixedClass
    

if __name__ == '__main__':
    ceuFile='data_simulated/ancestral_ceu_chr22.csv.gz'
    yriFile='data_simulated/ancestral_yri_chr22.csv.gz'
    aswFile='data_simulated/admixed_asw_chr22_8gens.csv.gz'
    aswCorrectFile='data_simulated/admixed_asw_chr22_8gens_origin.csv.gz'
    #mapping reference genome
    ceuFile='data_hapmap3/hapmap3_r2_b36_fwd.consensus.qc.poly.chr1_ceu.phased.gz'
    yriFile='data_hapmap3/hapmap3_r2_b36_fwd.consensus.qc.poly.chr1_yri.phased.gz'
    aswFile='data_hg18/hg18_chr1_hapmapIII.txt.gz'
    #read Genetic Map
    gm=geneticMap('data_simulated/genetic_map_chr1_b36.txt')
    #Cs=[]
    #Cs_svmSuccessMEAN=[]
    #Cs_svmSuccessSTD=[]
    #Cs_hmmSuccessMEAN=[]
    #Cs_hmmSuccessSTD=[]
    
    #for penalty in [1, 10, 100, 1000, 10000, 100000]:
    #    SVM=PyML.SVM(C=penalty, optimizer='mysmo')

    #Set up storage variables
    files=fileReader.concurrentFileReader(ceuFile, yriFile, aswFile)
    subjects=files.next()
    nTrain, nTest=len(subjects[0])+len(subjects[1]), len(subjects[2])  #Not Flexible
    labelsTrain,labelsTest =  ['ceu']*len(subjects[0]) + ['yri']*len(subjects[1]), ['yri']*len(subjects[2])
    idxTrain, idxTest=range(nTrain), range(nTrain, nTrain+nTest)
    snpLabels=[]        #stores snp labels from in files
    snpLocations=[]     #stores physical location from files
    ancestralSuccess=[] #stores success of ancestral classification
    admixedClass=[]     #stores classification of test Subjects
    mapLocations=[]
    vals=np.zeros((nTrain+nTest, WIN_SIZE))  #temporary storage of output

    while True: #for j in range(100): #To go through all positions in file
        for i, (snpName, snpLocation, snps) in enumerate(files):#range(WIN_SIZE):
            snpLabels.append(snpName)
            snpLocations.append(float(snpLocation))
            mapLocations.append(gm.pos2gm(snpLocations[-1]))
            vals[:,i]=fileReader.nucleotides2Haplotypes(sum(snps, []))
            #vals[:,i]=fileReader.nucleotides2SNPs(sum(snps, []))
            if i==vals.shape[1]-1:
                break
        ancestral, admixed=stepSVM(vals[idxTrain,:i], labelsTrain, 
                                   vals[idxTest, :i], labelsTest)
        ancestralSuccess.append(ancestral)
        admixedClass.append(admixed)
        if i<WIN_SIZE-1:
            break
    admixedClassPre=np.array(admixedClass)
    admixedClass=runHMM(mapLocations, ancestralSuccess, admixedClassPre)

    #Compare and find successRate
    #correct=np.array([l.split()[2:] for l in gzip.open(aswCorrectFile).readlines()[1:]], np.float)
    svmClass=np.repeat(admixedClassPre, WIN_SIZE, 0)
    hmmClass=np.repeat(admixedClass, WIN_SIZE, 0)

    #svmSuccess=100-sum(abs(svmClass[:len(correct),:]-correct))/len(correct)*100
    #hmmSuccess=100-sum(abs(hmmClass[:len(correct),:]-correct))/len(correct)*100

    #Visualize and collect output data
    #Cs.append(penalty)
    #Cs_svmSuccessMEAN.append(np.mean(svmSuccess))
    #Cs_svmSuccessSTD.append(np.std(svmSuccess))
    #Cs_hmmSuccessMEAN.append(np.mean(hmmSuccess))
    #Cs_hmmSuccessSTD.append(np.std(hmmSuccess))
    #print np.mean(svmSuccess), np.std(svmSuccess), np.mean(hmmSuccess), np.std(hmmSuccess)

    pylab.figure()
    pylab.clf()
    pylab.subplot(4,1,1);
    pylab.errorbar(range(len(ancestralSuccess)), map(np.mean, ancestralSuccess), map(np.std, ancestralSuccess))
    pylab.xlim([0, len(ancestralSuccess)]); pylab.xticks([])
    pylab.ylim([0,1]); pylab.ylabel('Success rate')
    pylab.title('Ancestral Population')


    pylab.subplot(4,1,2);
    pylab.imshow(admixedClassPre.T, interpolation='nearest')
    pylab.ylabel('Sample '); pylab.yticks([]); pylab.xticks([])
    pylab.axis('tight')

    pylab.subplot(4,1,3);
    pylab.imshow(admixedClass.T, interpolation='nearest')
    pylab.ylabel('Sample ');pylab.yticks([]); pylab.xticks([])
    pylab.axis('tight')

    #pylab.subplot(4,1,4);
    #pylab.imshow(correct[:, :].T, interpolation='nearest')
    #pylab.ylabel('Sample ');pylab.yticks([]); pylab.xticks([])
    #pylab.xlabel('Position along Chromosome 22')
    #pylab.axis('tight')

#     pylab.figure()
#     h1,h2,h3= pylab.errorbar(Cs, Cs_svmSuccessMEAN,  Cs_svmSuccessSTD)
#     h2,h3,h4=pylab.errorbar(Cs, Cs_hmmSuccessMEAN,  Cs_hmmSuccessSTD)
#     pylab.ylim([90, 110])
#     pylab.xlabel('Cs')
#     pylab.ylabel('Success Rate')
#     pylab.legend((h1, h2),['SVM', 'SVM and HMM'], 4)
#     ax=pylab.gca()
#     ax.set_xscale('log')
    pylab.show()
    
