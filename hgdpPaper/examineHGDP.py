import sys; sys.path.append('../')
import re, os, glob, cPickle, regionClassifier, numpy as np, fileReader, time, runLamp
from variables import *

def success(originFile, admixedClassPre, admixedClass, winSize=WINSIZE):
    correct=np.array([l.split()[2:] for l in fileReader.openfile(originFile).readlines()[1:]], np.float)
    #Compare and find successRate
    svmClass=np.repeat(admixedClassPre, winSize, 0)[:len(correct),:]
    hmmClass=np.repeat(admixedClass, winSize, 0)[:len(correct),:]
    svmSuccess=100*(svmClass==correct).sum(0)/float(len(correct))
    hmmSuccess=100*(hmmClass==correct).sum(0)/float(len(correct))
    return np.mean(hmmSuccess), np.std(hmmSuccess), np.mean(svmSuccess), np.std(svmSuccess)

def readFst():
    fst={}
    fp=fileReader.openfile(FILEFST+'.gz')
    for l in fp:
        pop1, pop2, val=l.strip().split('\t') 
        fst.setdefault(pop1, {})[pop2]=float(val)
        fst.setdefault(pop2, {})[pop1]=float(val)
    fp.close()
    return fst
    
def classify(fileNames, smoother, classifier=regionClassifier.SVMpymvpa(C), win_size=100):
    """Deconvolves ancestry in last file based on ancestral
    populations in first files.
    Arguments:
    - `fileNames`: list of fileNames
    - `win_size`: number of snps in each window (defualt=100)
    Returns:
    - `ancestralSuccess`: success of cross validation in ancestral populations
    - `admixedClassPre`: classification of admixed samples before hmm filter
    - `admixedClass`:    classification of admixed samples after hmm filter
    """
    snpLocations=[]     #stores physical location from files
    snpNames=[]     #stores physical location from files
    ancestralSuccess=[] #stores success of ancestral classification
    admixedClass=[]     #stores classification of test Subjects
    files=fileReader.concurrentFileReader(*fileNames, key=1)

    subjects=files.next()[0]
    nTrain=np.sum(map(len, subjects[:-1]))  #Number of samples in training set
    nTest=len(subjects[-1]);
    labelsTrain =sum([[i]*len(sub) for i, sub in enumerate(subjects[:-1])],[])
    vals=np.zeros((nTrain+nTest, win_size))  #temporary storage of output
    while True: 
        rsIds=[]
        pos=[]
        for i, ([snpName, snpLocation], snps) in enumerate(files):
            pos.append(float(snpLocation))
            rsIds.append(snpName)
            vals[:,i] = fileReader.nucleotides2Haplotypes(sum(snps, []))
            if i==vals.shape[1]-1:
                break
        snpLocations.append(pos)
        snpNames.append(rsIds)
        #print  len(snpLocations), snpLocations[-1][0], '->', snpLocations[-1][-1], snpNames[-1][0], '->', snpNames[-1][-1]
        ancestral, admixed=classifier(vals[:nTrain,:i+1], labelsTrain, vals[-nTest:, :i+1])
        ancestralSuccess.append(ancestral)
        admixedClass.append(admixed)
        if i<win_size-1:
            break
    admixedClassPre=np.array(admixedClass)
    admixedClass, p=smoother(np.hstack(snpLocations), ancestralSuccess, admixedClassPre)
    return  admixedClassPre, admixedClass, p, subjects[-1], snpLocations, snpNames


if __name__ == '__main__':
    ####################################
    #Find all ancestral populations and set up "global variables"
    ####################################
    ancestralFiles=glob.glob(FILEANCESTRAL %'*'+'.gz')
    pops=[re.findall(FILEANCESTRAL%'(.*)'+'.gz', fileName)[0] for fileName in ancestralFiles]
    pops.sort()
    fst=readFst()
    classifier=regionClassifier.SVMpymvpa(C)
    ####################################
    #Run all 2-way
    ####################################
    twoPopResults=results(WINSIZE, NGENS)
    smoother=regionClassifier.hmmFilter(geneticMapFile=GM_FILE,nGens=NGENS,nClasses=2)
    print NGENS, WINSIZE, C
    for pop1 in pops:
        for pop2 in pops:
            admixedFile=FILE2ADMPOPS%(pop1, pop2) +'.gz'
            if not  os.path.isfile(admixedFile): continue  #If file does not exist
            t0=time.time()
            originFile=FILE2ADMPOPSORIGIN%(pop1, pop2) + '.gz'
            fileNames=[FILEANCESTRAL%pop1+'.gz', FILEANCESTRAL%pop2+'.gz', admixedFile]
            admClassPre, admClass, p, subs, snpLocations, snpNames = classify(fileNames, smoother, classifier, WINSIZE)
            hmmSuccess, hmmStd, svmSuccess, svmStd =success(originFile, admClassPre, admClass)
            twoPopResults.append([pop1, pop2], fst[pop1][pop2], [hmmSuccess, hmmStd], p, admClass)
            print '%i:%i\t%s-%s\t%0.3g\t%0.3g' %((time.time()-t0)/60, (time.time()-t0)%60, pop1, pop2, fst[pop1][pop2], hmmSuccess)
    with open(OUTPUT_TWO_POP_SVM,'w') as fp: cPickle.dump(twoPopResults, fp)

    ####################################
    #Run Two way LAMP
    ####################################
    print '\nLAMP populations'
    twoPopResults=results(1, 1)
    for pop1, pop2 in POPS: 
        t0=time.time()
        admixedFile=FILE2ADMPOPS%(pop1, pop2) +'.gz'
        if not  os.path.isfile(admixedFile): continue  #If file does not exist
        ancestralFile1=FILEANCESTRAL%pop1+'.gz'
        ancestralFile2=FILEANCESTRAL%pop2+'.gz'
        originFile=FILE2ADMPOPSORIGIN%(pop1, pop2) + '.gz'
        runLamp.convertFiles(ancestralFile1, ancestralFile2, admixedFile)
        proc=runLamp.subprocess.Popen('lamp config.txt', shell=True, stdout=runLamp.subprocess.PIPE, stderr=runLamp.subprocess.STDOUT)
        tmp=proc.stdout.readlines()
        succMean, succStd, ancestry, correct= runLamp.readResults(originFile)
        runLamp.cleanUp()
        twoPopResults.append([pop1, pop2], fst[pop1][pop2], [succMean, succStd], [], ancestry)
        print '%i:%i\t%s-%s\t%0.3g\t%0.3g' %((time.time()-t0)/60, (time.time()-t0)%60, pop1, pop2, fst[pop1][pop2], succMean)
    with open(OUTPUT_TWO_POP_LAMP,'w') as fp: cPickle.dump(twoPopResults, fp)

    ####################################
    #Run all 3-way admixture
    ####################################
    print '\nThree populations'
    threePopResults=[results(WINSIZE, NGENS), results(WINSIZE, NGENS)]
    smoother=regionClassifier.hmmFilter(geneticMapFile=GM_FILE,nGens=NGENS,nClasses=3)
    for pop1 in pops:
        popNames=[['yoruba', 'french', pop1],['yoruba', 'bedouin', pop1]]
        admixedFile=[FILE3ADMPOPSCEUYRI%pop1+'.gz', FILE3ADMPOPSYORBED%pop1+'.gz']
        originFiles=[FILE3ADMPOPSCEUYRIORIGIN%pop1+'.gz', FILE3ADMPOPSYORBEDORIGIN%pop1+'.gz']
        for i, (anc1, anc2, anc3) in enumerate(popNames):
            if not  os.path.isfile(admixedFile[i]): continue  #If file does not exist
            t0=time.time()
            fileNames=[FILEANCESTRAL%anc1+'.gz', FILEANCESTRAL%anc2+'.gz', FILEANCESTRAL%anc3+'.gz', admixedFile[i]]
            admClassPre, admClass, p, subs, snpLocations, snpNames = classify(fileNames, smoother, classifier, WINSIZE)
            hmmSuccess, hmmStd, svmSuccess, svmStd =success(originFiles[i], admClassPre, admClass)
            currFst=min(fst[anc3][anc1], fst[anc3][anc2])
            threePopResults[i].append([anc1,anc2,anc3], currFst, [hmmSuccess, hmmStd], p, admClass)
            print '%i:%i\t%s-%s-%s\t%0.3g\t%0.3g' %((time.time()-t0)/60, (time.time()-t0)%60, anc1, anc2, anc3, currFst, hmmSuccess)
    with open(OUTPUT_THREE_AFRIC_SVM,'w') as fp: cPickle.dump(threePopResults[0], fp)
    with open(OUTPUT_THREE_ASIAN_SVM,'w') as fp: cPickle.dump(threePopResults[1], fp)

    ####################################
    #Run all Variable generations
    ####################################
    print '\nChanges in generations'
    twoPopResults=results(WINSIZE, NGENS)
    for pop1, pop2 in POPS:
        for nGens in SAMPLEGENERATIONS:
            t0=time.time()
            admixedFile=FILEGENSADMIX%(pop1, pop2, nGens)+'.gz'
            if not  os.path.isfile(admixedFile): continue  #If file does not exist
            fileNames=[FILEANCESTRAL%pop1+'.gz', FILEANCESTRAL%pop2+'.gz', admixedFile]
            originFile=FILEGENSADMIXORIGIN%(pop1, pop2, nGens)+'.gz'
            smoother=regionClassifier.hmmFilter(geneticMapFile=GM_FILE,nGens=nGens,nClasses=2)
            admClassPre, admClass, p, subs, snpLocations, snpNames = classify(fileNames, smoother, classifier, WINSIZE)
            hmmSuccess, hmmStd, svmSuccess, svmStd =success(originFile, admClassPre, admClass, WINSIZE)
            twoPopResults.append([pop1, pop2], nGens, [hmmSuccess, hmmStd], p, admClass)
            print '%i:%i\t%s-%s\t%0.3g\t%0.3g' %((time.time()-t0)/60, (time.time()-t0)%60, pop1, pop2, nGens, hmmSuccess)
    with open(OUTPUT_TWO_POP_SVM_GENS,'w') as fp: cPickle.dump(twoPopResults, fp)

    ####################################
    #Run all Variable alpha
    ####################################
    print '\nChanges in alpha'
    twoPopResults=results(WINSIZE, NGENS)
    smoother=regionClassifier.hmmFilter(geneticMapFile=GM_FILE,nGens=NGENS,nClasses=2)
    for pop1, pop2 in POPS:
        for a in ALPHAS:
            t0=time.time()
            admixedFile=FILEALPHAADMIX%(pop1, pop2, NGENS,a)+'.gz'
            if not  os.path.isfile(admixedFile): continue  #If file does not exist
            fileNames=[FILEANCESTRAL%pop1+'.gz', FILEANCESTRAL%pop2+'.gz', admixedFile]
            originFile=FILEALPHAADMIXORIGIN%(pop1, pop2, NGENS,a)+'.gz'
            admClassPre, admClass, p, subs, snpLocations, snpNames = classify(fileNames, smoother, classifier, WINSIZE)
            hmmSuccess, hmmStd, svmSuccess, svmStd =success(originFile, admClassPre, admClass)
            twoPopResults.append([pop1, pop2], a, [hmmSuccess, hmmStd], p, admClass)
            print '%i:%i\t%s-%s\t%0.3g\t%0.3g' %((time.time()-t0)/60, (time.time()-t0)%60, pop1, pop2, a, hmmSuccess)
    with open(OUTPUT_TWO_POP_SVM_ALPHA,'w') as fp: cPickle.dump(twoPopResults, fp)

    ####################################
    #Experiment with wrong generations
    ####################################
    print '\nChanges delta Generations'
    twoPopResults=results(WINSIZE, NGENS)
    for pop1, pop2 in POPS:
        for nGens in [.25, .5, 1.,5.,10.,20.,100.]:
            t0=time.time()
            admixedFile=FILEGENSADMIX%(pop1, pop2, 5)+'.gz'
            if not  os.path.isfile(admixedFile): continue  #If file does not exist
            fileNames=[FILEANCESTRAL%pop1+'.gz', FILEANCESTRAL%pop2+'.gz', admixedFile]
            originFile=FILEGENSADMIXORIGIN%(pop1, pop2, 5)+'.gz'
            smoother=regionClassifier.hmmFilter(geneticMapFile=GM_FILE,nGens=nGens,nClasses=2)
            admClassPre, admClass, p, subs, snpLocations, snpNames = classify(fileNames, smoother, classifier, WINSIZE)
            hmmSuccess, hmmStd, svmSuccess, svmStd =success(originFile, admClassPre, admClass)
            twoPopResults.append([pop1, pop2], 5/nGens, [hmmSuccess, hmmStd], p, admClass)
            print '%i:%i\t%s-%s\t%0.3g\t%0.3g' %((time.time()-t0)/60, (time.time()-t0)%60, pop1, pop2, 5/nGens, hmmSuccess)
    with open(OUTPUT_TWO_POP_SVM_DELTA_GENS,'w') as fp: cPickle.dump(twoPopResults, fp)
            

    ####################################
    #Experiment with window size
    ####################################
    print '\nChanges in window size with 5 generations'
    twoPopResults=results(WINSIZE, NGENS)
    smoother=regionClassifier.hmmFilter(geneticMapFile=GM_FILE,nGens=NGENS,nClasses=2)
    for pop1, pop2 in POPS:
        for winSize in [10, 50, 100, 200, 500, 1000, 2000, 5000]:
            t0=time.time()
            admixedFile=FILEGENSADMIX%(pop1, pop2, 5)+'.gz'
            if not  os.path.isfile(admixedFile): continue  #If file does not exist
            fileNames=[FILEANCESTRAL%pop1+'.gz', FILEANCESTRAL%pop2+'.gz', admixedFile]
            originFile=FILEGENSADMIXORIGIN%(pop1, pop2, 5)+'.gz'
            admClassPre, admClass, p, subs, snpLocations, snpNames = classify(fileNames, smoother, classifier, winSize)
            hmmSuccess, hmmStd, svmSuccess, svmStd =success(originFile, admClassPre, admClass, winSize)
            twoPopResults.append([pop1, pop2], winSize, [hmmSuccess, hmmStd], p, admClass)
            print '%i:%i\t%s-%s\t%0.3g\t%0.3g' %((time.time()-t0)/60, (time.time()-t0)%60, pop1, pop2, winSize, hmmSuccess)
    with open(OUTPUT_TWO_POP_SVM_WIN,'w') as fp: cPickle.dump(twoPopResults, fp)


    ####################################
    #Run Yoruba-Bedouin-Brahui with all populations
    ####################################
    classifier=regionClassifier.SVMpymvpa(C)
    C=1; NGENS=5;  WINSIZE=200; 
    fpClass=open('data/simulatedQatar.admixedClass.csv', 'w')
    fpP=open('data/simulatedQatar.posterior.csv', 'w')
    t0=time.time()
    admixedFile='data/hgdp3/admixed_hgdp_yoruba_bedouin_brahui.chr1.csv.gz'
    fileNames=glob.glob('data/hgdp_ancestral/ancestral_hgdp*')
    fileNames.sort()
    smoother=regionClassifier.hmmFilter('data/hapmap2/genetic_map_chr1_b36.txt', NGENS, len(fileNames))
    pops=[file.split('/')[2].replace('ancestral_hgdp_', '').replace('.chr1.csv.gz', '') for file in fileNames]
    fileNames.append(admixedFile)
    admClassPre, admClass, p, subs, snpLocations, snpNames = classify(fileNames, smoother, classifier, WINSIZE)
    print '%i:%i\t' %((time.time()-t0)/60, (time.time()-t0)%60)   
    #Save output
    for i in range(admClass.shape[0]): 
        fpClass.write('\t'.join(np.asarray(admClass[i,:], np.str))+'\n')
        fpP.write('\t'.join(np.asarray(p[i,:], np.str))+'\n')
    np.save('data/simulatedQatar.populations', pops)
    fpClass.close()
    fpP.close()
