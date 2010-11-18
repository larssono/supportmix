"""Region Classifier is a module with different classifiers for
genotype and haplotype classifiers.  Primarily based on support vector
machines (SVM)."""
import PyML, numpy as np
from  hmm import hmm
import mvpa.suite as pymvpa

class regionClassifier:
    def __init__(self):
        """Abstract constructor"""
        abstract()

    def __call__(self):
        """Trains classifier on valsTrains and lavelsTrain then
        tests the valsTest returning their classifications.
        """
        abstract()
    
class SVMpyml2(regionClassifier):
    def __init__(self, C=100):
        """SVM classifier for to classes using PyML
        Arguments:
        - `C`: Penalty term for missclassified samples in SVM
        """
        self.svm=PyML.SVM(C=C, optimizer='mysmo')

    def _createSVMData(self, vals, labels):
        """Given array and labels creates PyML with normalization and  Kernel"""
        dataPyML=PyML.VectorDataSet(vals, L=labels)
        dataPyML.normalize(1)
        #dataPyML.attachKernel('polynomial', degree = 1)
        dataPyML.attachKernel('linear')
        return dataPyML


    def __call__(self, valsTrain, labelsTrain, valsTest, doAncestralCV=True):
        """Does cv on ancestral population then trains on ancestral
        population followed by testing on admixed population.
        
        Arguments:
        - `valsTrain`: numpy array (nSamplesxnFeatures) of trainig samples 
        - `labelsTrain`: list of nSamples labels
        - `valsTest`:  numpy array of (nSamples2xnFeatures) of test samples
        """
        haplotypeData=self._createSVMData(valsTrain, labelsTrain)
        #train on ancestral populations
        self.svm.train(haplotypeData);
        #classify admixed population
        testData=self._createSVMData(valsTest,[1]*valsTest.shape[0])
        admixedClass = [self.svm.classify(testData,i)[0] for i in range(valsTest.shape[0])]
        if doAncestrallCV:
            results=self.svm.cv(haplotypeData, 3);
            ancestralSuccess = results.getBalancedSuccessRate()
            return ancestralSuccess, admixedClass
        return  admixedClass

class SVMpymvpa(regionClassifier):

    def __init__(self, C=100):
        """SVM classifier for to classes using PyML
        Arguments:
        - `C`: Penalty term for missclassified samples in SVM
        """
        #self.classifier = pymvpa.kNN(k=1, dfx=pymvpa.one_minus_correlation, voting='majority')
        self.classifier = pymvpa.LinearCSVMC(C=10)


    def __call__(self, valsTrain, labelsTrain, valsTest, doAncestralCV=True):
        """Does cv on ancestral population then trains on ancestral
        population followed by testing on admixed population.
        
        Arguments:
        - `valsTrain`: numpy array (nSamplesxnFeatures) of trainig samples 
        - `labelsTrain`: list of nSamples labels
        - `valsTest`:  numpy array of (nSamples2xnFeatures) of test samples
        """
        #Create and normalize data
        ds=pymvpa.Dataset(valsTrain)
        ds.sa['targets']=labelsTrain
        ds.sa['runtype']=np.random.randint(0, 2, valsTrain.shape[0])
        #pymvpa.zscore(ds, param_est=('targets', [1])) #Normalize somehow
        
        #Train on ancestral
        self.classifier.train(ds)
        admixedClass=self.classifier.predict(valsTest)
        #Cross Validated ancestral population
        if doAncestralCV:
            terr=pymvpa.TransferError(self.classifier)  #tracks error based on two datasets
            hspl = pymvpa.HalfSplitter(attr='runtype')
            cvte = pymvpa.CrossValidatedTransferError(terr, splitter=hspl)
            cv_results=cvte(ds)
            ancestralSuccess=1-np.mean(cv_results)
            return ancestralSuccess, admixedClass
        return admixedClass



#------------------------------------------------------------
# Post classfication Filtering
#------------------------------------------------------------
class globalFilter(object):
    """Given a region of classifiers filters based on some method"""

    def __init__(self):
        """abstract method"""
        abstract()

    def __call__(self):
        abstract()
        
class hmmFilter2(globalFilter):
    """Uses hmm and transition probabilites to filter previously
    classified regions """

    def __init__(self, geneticMapFile, nGens):
        """Constructor
        Arguments:
        - `geneticMapFile`: file containing mapping from physical distance to genetic distance
        - `nGens`: number of generations since admixture
        """
        self.gm=geneticMap(geneticMapFile)
        self.nGens=nGens

    def __call__(self,snpLocations, successRate, admixedClass):
        """Filters transitions based on hmm model 
        Arguments:
        - `snpLocations`: Locations of all the SNPs classified
        - `successRate`:  Probabilities of succesfully classifying each snp
        - `admixedClass`: classification made
        """
        mapLocations=map(self.gm.pos2gm, snpLocations)
        win_size=int(np.ceil(len(mapLocations)/float(admixedClass.shape[0]))) #TODO verify that this always works
        #determine transition matrices
        a=[]; b=[]
        oldPos=0
        for i in range(0,len(mapLocations), win_size):
            newPos=np.mean(mapLocations[i:i+win_size])
            dM=-(newPos - oldPos)/100*self.nGens
            e=np.exp(dM)  
            oldPos=newPos
            a.append([[e, 1-e],[1-e, e]])#TODO this will not work for snps with three classes
        a=np.asarray(a)
        for s in successRate:
            b.append([[s, 1-s],[1-s, s]])#TODO this will not work for snps with three classes
        b=np.asarray(b)
        model=hmm(a, b)

        #Go through and calculate hmm values
        results=[]
        for i in range(admixedClass.shape[1]):
            model.forward_backward(admixedClass[:,i])
            maxIdx=model.pstate.argsort(1)[:,-1]
            results.append(maxIdx)
        return np.array(results).T

class hmmFilter3(globalFilter):
    """Uses hmm and transition probabilites to filter previously
    classified regions """

    def __init__(self, geneticMapFile, nGens):
        """Constructor
        Arguments:
        - `geneticMapFile`: file containing mapping from physical distance to genetic distance
        - `nGens`: number of generations since admixture
        """
        self.gm=geneticMap(geneticMapFile)
        self.nGens=nGens

    def __call__(self,snpLocations, successRate, admixedClass):
        """Filters transitions based on hmm model 
        Arguments:
        - `snpLocations`: Locations of all the SNPs classified
        - `successRate`:  Probabilities of succesfully classifying each snp
        - `admixedClass`: classification made
        """
        mapLocations=map(self.gm.pos2gm, snpLocations)
        win_size=int(np.ceil(len(mapLocations)/float(admixedClass.shape[0]))) #TODO verify that this always works
        #determine transition matrices
        a=[]; b=[]
        oldPos=0
        for i in range(0,len(mapLocations), win_size):
            newPos=np.mean(mapLocations[i:i+win_size])
            dM=-(newPos - oldPos)/100*self.nGens
            e=np.exp(dM)  
            oldPos=newPos
            a.append([[e, (1-e)/2, (1-e)/2],
                      [(1-e)/2, e, (1-e)/2],
                      [(1-e)/2, (1-e)/2, e]])
        a=np.asarray(a)
        for s in successRate:
            b.append([[s,       (1-s)/2, (1-s)/2],
                      [(1-s)/2, s,       (1-s)/2],
                      [(1-s)/2, (1-s)/2, s]])#TODO this will not work for snps with three classes
        b=np.asarray(b)
        model=hmm(a, b)

        #Go through and calculate hmm values
        results=[]
        for i in range(admixedClass.shape[1]):
            model.forward_backward(admixedClass[:,i])
            maxIdx=model.pstate.argsort(1)[:,-1]
            results.append(maxIdx)
        return np.array(results).T
        

#--------------------------Helper functions--------------------------------------
def abstract():
    import inspect
    caller = inspect.getouterframes(inspect.currentframe())[1][3]
    raise NotImplementedError(caller + ' must be implemented in subclass')

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
