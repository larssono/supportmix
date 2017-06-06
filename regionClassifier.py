"""Region Classifier is a module with different classifiers for
genotype and haplotype classifiers.  Primarily based on support vector
machines (SVM)."""
import numpy as np
from  hmm import hmm
import mvpa2.suite as pymvpa


class regionClassifier:
    def __init__(self):
        """Abstract constructor"""
        abstract()

    def __call__(self):
        """Trains classifier on valsTrains and lavelsTrain then
        tests the valsTest returning their classifications.
        """
        abstract()
    
class SVMpymvpa(regionClassifier):

    def __init__(self, C=100):
        """SVM classifier for to classes using PyML
        Arguments:
        - `C`: Penalty term for missclassified samples in SVM
        """
        #self.classifier = pymvpa2.kNN(k=1, dfx=pymvpa.one_minus_correlation, voting='majority')
        self.classifier = pymvpa2.LinearCSVMC(C=C)


    def __call__(self, valsTrain, labelsTrain, valsTest, doAncestralCV=True):
        """Trains on ancestral population followed by testing on
        admixed population.  Optionally does cross validation on
        ancestral population.
        
        Arguments:
        - `valsTrain`: numpy array (nSamplesxnFeatures) of training samples 
        - `labelsTrain`: list of nSamples labels
        - `valsTest`:  numpy array of (nSamples2xnFeatures) of test samples
        """
        #Create and normalize data
        ds=pymvpa2.Dataset(valsTrain)
        ds.sa['targets']=labelsTrain
        runtype=np.zeros(valsTrain.shape[0]); runtype[0::3]=0;runtype[1::3]=1; runtype[2::3]=2 
        ds.sa['runtype']=runtype
        try:     #Train on ancestral
            self.classifier.train(ds)
            admixedClass=self.classifier.predict(valsTest)
        except pymvpa2.DegenerateInputError:  #The valsTrain is to small to contain information
            print "WARNING: Window is degenerate; guessing ancestry"
            admixedClass=np.zeros(valsTest.shape[0])  #Just assign ancestry to first pop
            if doAncestralCV:
                return 1./len(np.unique(labelsTrain)), admixedClass  #Assign success to create equal 
            return admixedClass
        if doAncestralCV:          #Cross Validated ancestral population
            hspl=pymvpa2.NGroupPartitioner(3, attr='runtype')
            # cvte = pymvpa2.CrossValidation(self.classifier, hspl)
            cvte = pymvpa2.CrossValidation(self.classifier, hspl, enable_ca='stats')
            cv_results=cvte(ds)
            return cvte.ca.stats.matrix, admixedClass 
            # ancestralSuccess=1-np.mean(cv_results)
            # return ancestralSuccess, admixedClass
        return admixedClass



#--------------------- Post classification Filtering ---------------------------
class globalFilter(object):
    """Given a region of classifiers filters based on some method"""
    def __init__(self):
        """abstract method"""
        abstract()

    def __call__(self):
        abstract()
        
class hmmFilter(globalFilter):
    """Uses hmm and transition probabilities to filter previously
    classified regions """

    def __init__(self, winSize, nGens, nClasses):
        """Constructor
        Arguments:
        - `winSize` - number of snps in each window
        - `nGens`: number of generations since admixture
        - `nClasses`: number of output classifications
        """
        self.winSize=winSize
        self.nGens=nGens
        self.nClasses=nClasses

    def __call__(self,mapLocations, successRate, admixedClass):
        """Filters transitions based on hmm model 
        Arguments:
        - `mapLocations`: Locations of all the SNPs classified in [cM]
        - `successRate`:  Probabilities of successfully classifying each snp
        - `admixedClass`: classification made
        """
        win_size=self.winSize  #int(np.ceil(len(mapLocations)/float(admixedClass.shape[0])))
        #determine transition matrices
        a=[]; b=[]
        oldPos=0
        for i in range(0,len(mapLocations), win_size):
            newPos=np.mean(mapLocations[i:i+win_size])
            dM=-(newPos - oldPos)/100*self.nGens
            e=np.exp(dM)  
            oldPos=newPos
            x=np.empty((self.nClasses, self.nClasses)) #Create state transitions
            for j in range(self.nClasses):
                x[j,:]=(1.-e)/(self.nClasses-1)
                x[j,j]=e
            a.append(x)

        a=np.asarray(a) 
        for s in successRate:
            # x=np.empty((self.nClasses, self.nClasses)) #Create output transitions
            # for i in range(self.nClasses):
            #     x[i,:]=(1.-s)/(self.nClasses-1)
            #     x[i,i]=s
            x=(s/s.sum(0).astype(np.float)).T #new1
            ##x=(s.T/s.sum(1).astype(float)).T  #new2
            b.append(x)
        b=np.asarray(b)
        model=hmm(a, b)

        #Go through and calculate hmm values
        results=[]
        posterior=[]
        for i in range(admixedClass.shape[1]):
            model.forward_backward(admixedClass[:,i])
            maxIdx=model.pstate.argsort(1)[:,-1]
            results.append(maxIdx)
            p=[np.asarray(model.pstate)[k][j] for (k,j) in enumerate(maxIdx)]
            posterior.append(p)
        return np.array(results).T, np.asarray(posterior).T

        
#--------------------------Helper functions--------------------------------------
def abstract():
    import inspect
    caller = inspect.getouterframes(inspect.currentframe())[1][3]
    raise NotImplementedError(caller + ' must be implemented in subclass')



if __name__ == '__main__':
    hm=hmmFilter('../../human_genome_data/data_hapmap3/genetic_map_chr22_b36.txt', 10, 3)
    snpLocations=[14431347,16211813,16989647,17960185,19135369,20429295,21053116,21773909,22522081,23704434,24530858,25026783,25592569,25979024,26384814,27576605,28538608,29525941,30726524,31319141,31887585,32430944,32932547,33457979,33936496,34908093,35571213,36162392,37337120,38357467,39774632,40982261,41820819,42478803,43037663,43632337,44287817,45168693,45771245,46338815,46775561,47215631,47719425,48069691,48554269]
    successRate=np.asarray([.3,]*len(snpLocations))
    admixedClass=np.asarray([[0,]*20+[1,]*25]).T
    x = hm(snpLocations, successRate, admixedClass)
    print x
