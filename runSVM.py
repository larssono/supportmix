import PyML, pylab
import numpy as np
import fileReader
import sys

WIN_SIZE=100
SVM=PyML.SVM(C=10, optimizer='mysmo')

def runHMM(snpLocations, admixedClass):
    """Smoothes out admixedClass using HMM"""
    import hmm
    states=[0,1]
    start_probability = {1: 0.75, 0: 0.25}
    transition_probability = {0 : {0: 0.7, 1: 0.3},
                              1 : {0: 0.3, 1: 0.7}}
    emission_probability = {0 : {0: 0.8, 1: 0.2},
                            1 : {0: 0.2, 1: 0.8}}
    results=[]
    for i in range(admixedClass.shape[1]):
        p1, states, p2 = hmm.forward_viterbi(admixedClass[:,i],
                                             states, start_probability, 
                                             transition_probability,
                                             emission_probability)
        results.append(states)
    return np.array(results).T


def createSVMData(vals, labels):
    """Given array and labels creates PyML with normalization and  Kernel"""
    dataPyML=PyML.VectorDataSet(vals, L=labels)
    dataPyML.normalize(1)
    #dataPyML.attachKernel('polynomial', degree = 1)
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
    

def trainAndTestSVM(fileNames):
    """Iterates through files reading WIN_SIZE number of snps then trains and tests
    """
    files=fileReader.concurrentFileReader(fileNames['ceu'], fileNames['yri'], fileNames['asw'])
    subjects=files.next()
    nTrain, nTest=len(subjects[0])+len(subjects[1]), len(subjects[2])  #Not Flexible
    labelsTrain,labelsTest =  ['ceu']*len(subjects[0]) + ['yri']*len(subjects[1]), ['yri']*len(subjects[2])
    idxTrain, idxTest=range(nTrain), range(nTrain, nTrain+nTest)
    snpLabels=[]        #stores snp labels from in files
    snpLocations=[]     #stores physical location from files
    ancestralSuccess=[] #stores success of ancestral classification
    admixedClass=[]     #stores classification of test Subjects
    vals=np.zeros((nTrain+nTest, WIN_SIZE))
    print vals.shape, nTrain, nTest

    #while True:  #To go through all positions in file
    for j in range(200):
        for i, (snpName, snpLocation, snps) in enumerate(files):#range(WIN_SIZE):
            snpLabels.append(snpName)
            snpLocations.append(snpLocation)
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
    admixedClass=np.array(admixedClass)
    admixedClass=runHMM(snpLocations, admixedClass)
    return snpLabels, snpLocations, ancestralSuccess, admixedClass


if __name__ == '__main__':
    fileNames={'ceu': 'data_hapmap3/hapmap3_r2_b36_fwd.consensus.qc.poly.chr22_ceu.phased.gz',
           'yri': 'data_hapmap3/hapmap3_r2_b36_fwd.consensus.qc.poly.chr22_yri.phased.gz',
           'asw': 'data_abra_test/LD_0.80.chr22.filtered.txt.gz'}
           #'asw': 'data_hapmap3/hapmap3_r2_b36_fwd.consensus.qc.poly.chr22_asw.phased.gz'}
    snpLabels, snpLocations, ancestralSuccess, admixedClass =  trainAndTestSVM(fileNames)

    pylab.clf()
    pylab.subplot(2,1,1);
    pylab.errorbar(range(len(ancestralSuccess)), map(np.mean, ancestralSuccess), map(np.std, ancestralSuccess))
    pylab.title('Success of SVM on ancestral population')

    pylab.subplot(2,1,2);
    pylab.imshow(admixedClass.T, interpolation='nearest')
    pylab.ylabel('Sample ');pylab.title('Genotype Origin')
    pylab.axis('tight')
    pos, vals=pylab.xticks()

    pylab.figure()
    pylab.hist(admixedClass.sum(0), 10) #float(NSTEPS), 10)
    pylab.title('Distribution of Percent African in Population')
    pylab.show()
