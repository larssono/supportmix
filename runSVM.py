import PyML, pylab
import numpy as np
import gzip

def readSNPs(file):
    """Reads first and second column of file storing SNP names and position in bp
    """
    fp=gzip.open(file); fp.readline()
    snps = np.asarray([l.split()[0:2] for l in fp], dtype=np.str_)
    snpNames=snps[:,0]
    snpLocations=snps[:,1].astype(int)
    return snpNames, snpLocations

    

def readChrFile(file):
    """Reads one hapmap3 phased chromosome file
    Returns subjectList, SNPs
    """
    fp=gzip.open(file)
    subjectList=np.asarray(fp.readline().strip().split()[2:], dtype=np.str_)
    snps=[line.strip().split()[2:] for line in fp]
    return subjectList, snps, 

def convertSNPs(snps):
    """Recieves a list of lists of nucleotiedes and converts to numpy array of 0,1
    """
    nSubj=len(snps[0])
    out=np.zeros((len(snps), nSubj), dtype=int)

    for i, snp in enumerate(snps):
        genotypes=list(set(snp))
        if snp.count(genotypes[0]) >nSubj/2.:
            majorAllele=genotypes[0]
        else:
            majorAllele=genotypes[1]
    
        out[i,:]=np.asarray(snp)==majorAllele
    return out
    

if __name__ == '__main__':
    STEPSIZE=200
    NSTEPS=582
    asw='hapmap3_r2_b36_fwd.consensus.qc.poly.chr1_asw.phased.gz'
    ceu='hapmap3_r2_b36_fwd.consensus.qc.poly.chr1_ceu.phased.gz'
    yri='hapmap3_r2_b36_fwd.consensus.qc.poly.chr1_yri.phased.gz'
    LABELS=['asw','ceu', 'yri']
    classes=[]
    allSubjects=[]
    allSNPs=[]
    for i,file in  enumerate([asw, ceu, yri]):
        subj, snps=readChrFile(file)
        allSubjects.extend(subj)
        classes.extend([LABELS[i]]*len(subj))
        if len(allSNPs)==0:
            allSNPs=snps
        else:
            for i, snp in enumerate(snps):
                allSNPs[i].extend(snp)
    data=convertSNPs(allSNPs)
    data[data==0] = -1
    snpNames, snpLocations = readSNPs(asw)
    classes=np.asarray(classes)
    nSNPs, nSubjects=data.shape
    
    idxNASW = classes!='asw'
    idxASW =  classes=='asw'
    
    s=PyML.SVM(C=100)
    
    ancestralSuccess=[]
    admixedYRI=[]
    for i in range(NSTEPS):
        #kernel=np.corrcoef(data[i*STEPSIZE:(i+1)*STEPSIZE, idxNASW].T)
        haplotypeData=PyML.VectorDataSet(data[i*STEPSIZE:(i+1)*STEPSIZE, idxNASW].T, L=list(classes[idxNASW]))
        haplotypeData.normalize(1)
        #haplotypeData.attachKernel('polynomial', degree = 1, normalization = 'cosine')

        result=s.cv(haplotypeData, 5)
        ancestralSuccess.append(result.getBalancedSuccessRate())
        
        s.train(haplotypeData)
        testData=PyML.VectorDataSet(data[i*STEPSIZE:(i+1)*STEPSIZE, idxASW].T, L=['ceu']*sum(idxASW))
        testData.normalize(1)
        #testData.attachKernel('polynomial', degree = 1, normalization = 'cosine')

        result = s.test(testData)
        print result
        admixedYRI.append(result.getPredictedClass())

    admixedYRI=np.asarray(admixedYRI)

    pylab.subplot(2,1,1);
    pylab.plot(ancestralSuccess); pylab.xlim([0,NSTEPS])
    pos, vals=pylab.xticks()    
    pylab.xticks(pos[1:-1], (snpLocations[(pos[1:-1]*STEPSIZE).astype(int)]/1e6).round())
    pylab.title('Success of SVM on ancestral population')
    pylab.ylim([0, 1])
    pylab.subplot(2,1,2);
    pylab.imshow(admixedYRI.T, interpolation='nearest')
    pylab.ylabel('Sample ');pylab.title('Genotype Origin')
    pylab.axis('tight')
    pos, vals=pylab.xticks()    
    pylab.xticks(pos[1:-1], (snpLocations[(pos[1:-1]*STEPSIZE).astype(int)]/1e6).round())
    pylab.xlabel('Mb along chromosome 1')

    pylab.figure()
    pylab.hist(admixedYRI.sum(0)/float(NSTEPS), 5)
    pylab.title('Distribution of Percent African in Population')
    
