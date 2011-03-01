import sys; sys.path.append('../')
import glob, runSVM, numpy as np, fileReader, time
NGENS=1
NWIN=100
ADMIXEDNAMES='hgdp/admixed_hgdp_%s_%s.chr1.csv.gz'
CHROM=1
C=1000

def success(originFile, admixedClassPre, admixedClass):
    correct=np.array([l.split()[2:] for l in fileReader.openfile(originFile).readlines()[1:]], np.float)
    #Compare and find successRate
    svmClass=np.repeat(admixedClassPre, WINSIZE, 0)[:len(correct),:]
    hmmClass=np.repeat(admixedClass, WINSIZE, 0)[:len(correct),:]
    svmSuccess=100*(svmClass==correct).sum(0)/float(len(correct))
    hmmSuccess=100*(hmmClass==correct).sum(0)/float(len(correct))
    return np.mean(hmmSuccess), np.std(hmmSuccess), np.mean(svmSuccess), np.std(svmSuccess)

summaryOutFile=open('summary_WIN%i_GENS%i_C%i.txt' %(NWIN, NGENS, C), 'w')
summaryOutFile.write('pop1\tpop2\twin[bp]\twin[cm]\tsuccess[svm]\tstd[svm]\tsuccess[hmm]\tstd[hmm]\n')

admixedFiles=glob.glob('hgdp/admixed_hgdp*')
ancestralFiles=glob.glob('hgdp/ancestral_hgdp*')
t0=time.time()
for ancestral1 in ancestralFiles:
    for ancestral2 in ancestralFiles:
        pop1=ancestral1.replace('hgdp/ancestral_hgdp_', '').replace('.chr1.csv.gz', '')
        pop2=ancestral2.replace('hgdp/ancestral_hgdp_', '').replace('.chr1.csv.gz', '')
        admixedFile=ADMIXEDNAMES%(pop1, pop2)
        if admixedFile in admixedFiles:
            originFile=admixedFile.replace('_hgdp_', '_origin_hgdp_')
            fileNames=[ancestral1, ancestral2, admixedFile]
            outFile='tmp/%s_%s_win%i_gens%i_c%i.csv' %(pop1, pop2, NWIN, NGENS, C)
            snpLocations, ancestralSuccess, admixedClassPre, admixedClass = runSVM.runSVM(fileNames, 
                                                                                          nGens=NGENS, 
                                                                                          win_size=NWIN,
                                                                                          svmC=C)
            hmmSuccess, hmmStd, svmSuccess, svmStd =success(originFile, admixedClassPre, admixedClass)
            winSizeBP, winSizeCM = runSVM.winSizeBPandCM(snpLocations, NWIN, CHROM)
            #output summary data
            print '%i:%i' %((time.time()-t0)/60, (time.time()-t0)%60),
            print '%s(%2.2g) + %s(%2.2g)' %(pop1, 100*np.sum(admixedClass.flatten()==0)/float(np.prod(admixedClass.shape)), pop2, 100*np.sum(admixedClass.flatten()==1)/float(np.prod(admixedClass.shape))),
            print 'Correct: %0.3g+/-%0.2g (%0.3g+/-%0.2g)' %(hmmSuccess, hmmStd, svmSuccess,svmStd)
            summaryOutFile.write('%s\t%s\t%0.4g\t%0.4g\t%0.4g\t%0.4g\t%0.4g\t%0.4g\n' %(pop1, pop2, 
                                                                                       winSizeBP.mean(), 
                                                                                       winSizeCM.mean(), 
                                                                                       svmSuccess, 
                                                                                       svmStd, 
                                                                                       hmmSuccess, 
                                                                                       hmmStd ))
            summaryOutFile.flush()
            #Output detail results
            fp=open(outFile, 'w')  #Save output to file
            starts=snpLocations[::NWIN]
            fp.write('Chrom\t Position\t AncestralSuccess\t WinSize[BP]\t WinSize[CM]\t AdmixedClass\n')
            for i in range(len(starts)):
                fp.write('%s\t%s\t%s\t%i\t%g\t' %(CHROM, starts[i], ancestralSuccess[i], winSizeBP[i], winSizeCM[i]))
                fp.write('\t'.join(np.asarray(admixedClass[i,:], np.str_)))
                fp.write('\n')
            fp.close()
summaryOutFile.close()
