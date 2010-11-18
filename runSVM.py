import fileReader, regionClassifier, sys, gzip,  pylab, numpy as np
from optparse import OptionParser

MAPFILES='/'.join(__file__.split('/')[:-1])+'/../../human_genome_data/hapmap2/genetic_map_%s_b36.txt'
USAGE="""%prog [options] ancestralFile1 ancestralFile2 [...] admixedFile

where the ancestralFiles and admixedFiles are contain phased samples
in tab delimited format snps in rows and samples in columns.  One row
of column headers and two columns of header information (i.e. rsId and
position).  The files have to be from the same chromosome indicated by
including chr[0-9]* in the name.
"""

def winSizeBPandCM(snpLocations, winSize, chrom):
    chrom='chr%i'%chrom
    gm=regionClassifier.geneticMap(MAPFILES%chrom)
    winStarts=[0]; winStarts.extend(snpLocations[::winSize])
    gmPos=map(gm.pos2gm, winStarts)
    return np.diff(winStarts), np.diff(gmPos)

def determineChromosome(fileNames):
    """Estimates the chromosome name by searching for "chr" in input
    file names.
    Arguments:
    - `fileNames`:List of fileNames 
    """
    import re 
    p = re.compile('chr\d*', re.IGNORECASE)
    try:
        found=[p.search(name).group() for name in fileNames]
        if not np.all(found[0]==np.asarray(found)): raise Error
        return found[0] 
    except:
        print "ERROR: names of input files should contain the same chr[0-9]*"
        sys.exit()

def runSVM(fileNames, nGens=6, svmC=100, win_size=100):
    """Deconvolves ancestry in last file based on ancestral
    populations in first files.
    Arguments:
    - `fileNames`: list of fileNames
    - `nGens`: Number of generations since admixture used in hmm filter (default=6)
    - `svmC`: missclassification penalty term in svm (default=100)
    - `win_size`: number of snps in each window (defualt=100)

    Returns:
    - `snpLocations`: list of positions in bp
    - `ancestralSuccess`: success of cross validation in ancestral populations
    - `admixedClassPre`: classification of admixed samples before hmm filter
    - `admixedClass`:    classification of admixed samples after hmm filter
    """
    files=fileReader.concurrentFileReader(*fileNames, key=1)
    chrom=determineChromosome(fileNames)
    snpLocations=[]     #stores physical location from files
    ancestralSuccess=[] #stores success of ancestral classification
    admixedClass=[]     #stores classification of test Subjects
    classifier=regionClassifier.SVMpymvpa(C=svmC)
    smoother=regionClassifier.hmmFilter(geneticMapFile=MAPFILES%chrom,nGens=nGens,nClasses=len(fileNames)-1)

    subjects=files.next()
    nTrain=np.sum(map(len, subjects[:-1]))  #Number of samples in training set
    nTest=len(subjects[-1]);
    labelsTrain =sum([[i]*len(sub) for i, sub in enumerate(subjects[:-1])],[])
    vals=np.zeros((nTrain+nTest, win_size))  #temporary storage of output
    while True: #for j in range(100): #To go through all positions in file
        for i, (snpName, snpLocation, snps) in enumerate(files):
            snpLocations.append(float(snpLocation))
            vals[:,i] = fileReader.nucleotides2Haplotypes(sum(snps, []))
            if i==vals.shape[1]-1:
                break
        ancestral, admixed=classifier(vals[:nTrain,:i], labelsTrain, vals[-nTest:, :i])
        ancestralSuccess.append(ancestral)
        admixedClass.append(admixed)
        if i<win_size-1:
            break
    admixedClassPre=np.array(admixedClass)
    admixedClass=smoother(snpLocations, ancestralSuccess, admixedClassPre)

    return snpLocations, ancestralSuccess, admixedClassPre, admixedClass


if __name__ == '__main__':
    parser = OptionParser(usage=USAGE)
    #I want to add stuff describing the arguments
    parser.add_option("-f", "--fileCorrect", type='string', dest="correctFile", 
                      help="FILE contains correct classifications in tab delimited format with SNPs in rows and samples in columns (first two columns contain rsID and position. ", metavar="FILE")
    parser.add_option("-w", "--window", type='int', dest="win", default=10, 
                      help="Number of SNPs in each window (default 10)", metavar="N")
    parser.add_option("-g", "--generations", type='float', dest="nGens", default=6,
                      help="Number of generations sinces admixture used in hmm. (default 6)", metavar="N")
    parser.add_option("-s", "--save", type='str', dest="saveFile", default=None,
                      help="Destination file to save output", metavar="N")
    (options, args) = parser.parse_args()

    fileNames=args
    snpLocations, ancestralSuccess, admixedClassPre, admixedClass = runSVM(fileNames, 
                                                                           nGens=options.nGens, 
                                                                           win_size=options.win)
    for i, pop in enumerate(fileNames[:-1]):
        pop=pop.split('/')[-1]
        print 'Pop %s: %2.2g' %(pop, 100*np.sum(admixedClass.flatten()==i)/float(np.prod(admixedClass.shape)))

    chrom=determineChromosome(fileNames)
    pylab.subplot(2,1,1);
    pylab.errorbar(range(len(ancestralSuccess)), map(np.mean, ancestralSuccess), map(np.std, ancestralSuccess))
    pylab.xlim([0, len(ancestralSuccess)]); pylab.xticks([])
    pylab.ylim([0,1]); pylab.ylabel('Success rate')
    pylab.title('Ancestral Population Success Rate')
    pylab.subplot(2,1,2);
    pylab.imshow(admixedClass.T, interpolation='nearest', 
                 cmap=pylab.cm.copper, vmin=0, vmax=2)
    pylab.ylabel('Sample ');pylab.yticks([]); pylab.xticks([]); pylab.axis('tight'); 
    pylab.xlabel('Position along %s' %chrom)
    pylab.title('Admixed Classification')


    if options.saveFile:
        fp=open(options.saveFile, 'w')  #Save output to file
        starts=snpLocations[::options.win]
        fp.write('Chrom\t Position\t AncestralSuccess\t AdmixedClass\n')
        for i in range(len(starts)):
            fp.write('%s\t%s\t%s\t' %(chrom, starts[i], ancestralSuccess[i]))
            fp.write('\t'.join(np.asarray(admixedClass[i,:], np.str_)))
            fp.write('\n')
        fp.close()


    if options.correctFile:
        correct=np.array([l.split()[2:] for l in fileReader.openfile(options.correctFile).readlines()[1:]], np.float)
        #Compare and find successRate
        svmClass=np.repeat(admixedClassPre, options.win, 0)
        hmmClass=np.repeat(admixedClass, options.win, 0)
        svmSuccess=100-sum(abs(svmClass[:len(correct),:]-correct))/len(correct)*100
        hmmSuccess=100-sum(abs(hmmClass[:len(correct),:]-correct))/len(correct)*100
        print 'Correct %0.3g +/- %0.2g (%0.3g +/- %0.2g)' %(np.mean(hmmSuccess), np.std(hmmSuccess), 
                                                            np.mean(svmSuccess), np.std(svmSuccess))
        pylab.figure();
        pylab.subplot(2,1,1);
        pylab.imshow(admixedClass.T, interpolation='nearest', 
                     cmap=pylab.cm.copper, vmin=0, vmax=2)
        pylab.ylabel('Sample '); pylab.yticks([]); pylab.xticks([]); pylab.axis('tight'); pylab.title('Estimated Ancestry')
        pylab.subplot(2,1,2);
        pylab.imshow(correct[:, :].T, interpolation='nearest', cmap=pylab.cm.copper, vmin=0, vmax=2)
        pylab.ylabel('Sample');pylab.yticks([]); pylab.xticks([])
        pylab.xlabel('Position along %s' %chrom);  pylab.axis('tight'); pylab.title('True ancestry')
        pylab.show()
    pylab.show()
