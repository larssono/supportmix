#Test different genetic maps
from variables import *
from examineHGDP import *

WINSIZE=100
twoPopResults=results(WINSIZE, NGENS)
smoother=regionClassifier.hmmFilter(WINSIZE,nGens=NGENS,nClasses=2)
classifier=regionClassifier.SVMpymvpa(C)
for pop1, pop2 in POPS:
    admixedFile=FILE2ADMPOPS%(pop1, pop2) +'.gz'
    t0=time.time()
    originFile=FILE2ADMPOPSORIGIN%(pop1, pop2) + '.gz'
    fileNames=[FILEANCESTRAL%pop1+'.gz', FILEANCESTRAL%pop2+'.gz', admixedFile]
    admClassPre, admClass, p, subs, snpLocations, snpNames = classify(fileNames, smoother, classifier, WINSIZE, CHR, mapFile='data/hapmap2/genetic_map_%(CHR)s_b36.txt')
    hmmSuccess, hmmStd, svmSuccess, svmStd =success(originFile, admClassPre, admClass, WINSIZE)
    print '%i:%i\t%s-%s\t%0.3g' %((time.time()-t0)/60, (time.time()-t0)%60, pop1, pop2, hmmSuccess)

    admClassPre, admClass, p, subs, snpLocations, snpNames = classify(fileNames, smoother, classifier, WINSIZE, CHR, mapFile='data/genetic_map_const_chr1_.csv')
    hmmSuccess, hmmStd, svmSuccess, svmStd =success(originFile, admClassPre, admClass, WINSIZE)
    print '%i:%i\t%s-%s\t%0.3g' %((time.time()-t0)/60, (time.time()-t0)%60, pop1, pop2, hmmSuccess)


# french-bedouin	77.0   	76.1
# bedouin-yoruba	97.7 	97.7
# han-bedouin	97.8    97.7
# french-yoruba	98.5 	98.3
# han-yoruba	99.1    98.9
# papuan-yoruba	99.0	99.0
# papuan-karitiana	98.1	98.2
