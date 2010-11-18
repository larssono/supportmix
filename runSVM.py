import fileReader, regionClassifier, sys, gzip,  pylab, numpy as np

WIN_SIZE=150

FILES=['data_simulated/yri_mkk/ancestral_yri.chr22.csv',
       'data_simulated/yri_mkk/ancestral_mkk.chr22.csv',
       'data_simulated/yri_mkk/admixed_yri_mkk.chr22.csv']
ADMIXED_CORRECT_FILE='data_simulated/yri_mkk/admixed_origin_yri_mkk.chr22.csv'
# FILES=['data_simulated/ceu_yri_chd/ancestral_ceu.chr22.csv',
#        'data_simulated/ceu_yri_chd/ancestral_yri.chr22.csv',
#        'data_simulated/ceu_yri_chd/ancestral_chd.chr22.csv',
#        'data_simulated/ceu_yri_chd/admixed_ceu_yri_chd.chr22.csv']
# ADMIXED_CORRECT_FILE='data_simulated/ceu_yri_chd/admixed_origin_ceu_yri_chd.chr22.csv'

MAPFILE='../../human_genome_data/genetic_map_chr22_b36.txt'

def readCorrect():
    """Read the file of correct classifications and returns
    the classes or the genotyped class."""
    fileData=np.array([l.split()[2:] for l in fileReader.openfile(ADMIXED_CORRECT_FILE).readlines()[1:]], np.float)
    return fileData

if __name__ == '__main__':
    snpLocations=[]     #stores physical location from files
    ancestralSuccess=[] #stores success of ancestral classification
    admixedClass=[]     #stores classification of test Subjects

    #classifier=regionClassifier.SVMpyml2(C=100)
    classifier=regionClassifier.SVMpymvpa(C=100)

    filter=regionClassifier.hmmFilter3(geneticMapFile=MAPFILE,nGens=1)
    files=fileReader.concurrentFileReader(*FILES)
    subjects=files.next()
    nTrain=np.sum(map(len, subjects[:-1]))  #Number of samples in training set
    nTest=len(subjects[-1]);
    labelsTrain =sum([[i]*len(sub) for i, sub in enumerate(subjects[:-1])],[])   #Arbitrary labels for ancestral populations

    vals=np.zeros((nTrain+nTest, WIN_SIZE))  #temporary storage of output
    while True: #for j in range(100): #To go through all positions in file
        for i, (snpName, snpLocation, snps) in enumerate(files):
            snpLocations.append(float(snpLocation))
            vals[:,i] = fileReader.nucleotides2Haplotypes(sum(snps, []))
            if i==vals.shape[1]-1:
                break
        ancestral, admixed=classifier(vals[:nTrain,:i], labelsTrain, vals[-nTest:, :i])
        ancestralSuccess.append(ancestral)
        admixedClass.append(admixed)
        if i<WIN_SIZE-1:
            break
    admixedClassPre=np.array(admixedClass)
    admixedClass=filter(snpLocations, ancestralSuccess, admixedClassPre)

    #Compare and find successRate
    correct=readCorrect()
    svmClass=np.repeat(admixedClassPre, WIN_SIZE, 0)
    hmmClass=np.repeat(admixedClass, WIN_SIZE, 0)

    svmSuccess=100-sum(abs(svmClass[:len(correct),:]-correct))/len(correct)*100
    hmmSuccess=100-sum(abs(hmmClass[:len(correct),:]-correct))/len(correct)*100
    print np.mean(svmSuccess),'+/-', np.std(svmSuccess), np.mean(hmmSuccess),'+/-',np.std(hmmSuccess)

    pylab.figure(1); pylab.clf()
    pylab.subplot(4,1,1);
    pylab.errorbar(range(len(ancestralSuccess)), map(np.mean, ancestralSuccess), map(np.std, ancestralSuccess))
    pylab.xlim([0, len(ancestralSuccess)]); pylab.xticks([])
    pylab.ylim([0,1]); pylab.ylabel('Success rate')
    pylab.title('Ancestral Population')

    pylab.subplot(4,1,2);
    pylab.imshow(admixedClassPre.T, interpolation='nearest', cmap=pylab.cm.copper, vmin=0, vmax=2)
    pylab.ylabel('Sample '); pylab.yticks([]); pylab.xticks([]); pylab.axis('tight')

    pylab.subplot(4,1,3);
    pylab.imshow(admixedClass.T, interpolation='nearest', cmap=pylab.cm.copper, vmin=0, vmax=2)
    pylab.ylabel('Sample ');pylab.yticks([]); pylab.xticks([]); pylab.axis('tight')

    pylab.subplot(4,1,4);
    pylab.imshow(correct[:, :].T, interpolation='nearest', cmap=pylab.cm.copper, vmin=0, vmax=2)
    pylab.ylabel('Sample ');pylab.yticks([]); pylab.xticks([])
    pylab.xlabel('Position along Chromosome 1');  pylab.axis('tight')
    pylab.show()
    
