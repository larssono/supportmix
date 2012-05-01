import runLamp
import numpy as np
import regionClassifier
from variables import *
from examineHGDP import classify, success

ancestralFile1='ancestral_hgdp_french.chr1_short.csv.gz'
ancestralFile2='ancestral_hgdp_yoruba.chr1_short.csv.gz'
admixedFile='admixed_hgdp_french_yoruba.chr1_short.csv.gz'
originFile='admixed_origin_hgdp_french_yoruba.chr1_short.csv.gz'
fileNames=[ancestralFile1, ancestralFile2, admixedFile]


#Run SupportMix on two shortened populations
classifier=regionClassifier.SVMpymvpa(C)
smoother=regionClassifier.hmmFilter(WINSIZE,nGens=NGENS,nClasses=2)
admClassPre, admClass, p, subs, snpLocations, snpNames = classify(fileNames, smoother, classifier, WINSIZE, CHR)
hmmSuccess, hmmStd, svmSuccess, svmStd =success(originFile, admClassPre, admClass)
print 'SupportMix:  %0.3g+/-%0.3g' %(hmmSuccess, hmmStd)


#Run LAMP on two shortened populations
runLamp.convertFiles(ancestralFile1, ancestralFile2, admixedFile)
proc=runLamp.subprocess.Popen('lamp config.txt', shell=True, stdout=runLamp.subprocess.PIPE, stderr=runLamp.subprocess.STDOUT)
tmp=proc.stdout.readlines()
succMean, succStd, ancestry, correct= runLamp.readResults(originFile)
runLamp.cleanUp()
print 'LAMP:  %0.3g+/-%0.3g' %(succMean, succStd)



