CHR='chr1'
GM_FILE='data/hapmap2/genetic_map_%(CHR)s_b36.txt'%locals()
FILEFST='data/hgdp_fst.%(CHR)s.csv'%locals()
HGDPCOORDSFILE='data/rosenbergEtAl2005.coordinates.txt' #from http://rosenberglab.bioinformatics.med.umich.edu/diversity.html
FILEANCESTRAL='data/hgdp_ancestral/ancestral_hgdp_%s.'+CHR+'.csv'
FILE2ADMPOPS='data/hgdp2/admixed_hgdp_%s_%s.'+CHR+'.csv'
FILE2ADMPOPSORIGIN='data/hgdp2/admixed_origin_hgdp_%s_%s.'+CHR+'.csv'
FILE3ADMPOPSCEUYRI='data/hgdp3/admixed_hgdp_yoruba_french_%s.'+CHR+'.csv'
FILE3ADMPOPSCEUYRIORIGIN='data/hgdp3/admixed_hgdp_origin_yoruba_french_%s.'+CHR+'.csv'
FILE3ADMPOPSHANJPN='data/hgdp3/admixed_hgdp_han_japanese_%s.'+CHR+'.csv'
FILE3ADMPOPSHANJPNORIGIN='data/hgdp3/admixed_hgdp_origin_han_japanese_%s.'+CHR+'.csv'
FILEGENSADMIX='data/hgdp_generations/admixed_hgdp_%s_%s_%ggens.'+CHR+'.csv'
FILEGENSADMIXORIGIN='data/hgdp_generations/admixed_origin_hgdp_%s_%s_%ggens.'+CHR+'.csv'
FILEALPHAADMIX='data/hgdp_alpha/admixed_hgdp_%s_%s_%ggens_%galpha.'+CHR+'.csv'
FILEALPHAADMIXORIGIN='data/hgdp_alpha/admixed_origin_hgdp_%s_%s_%ggens_%galpha.'+CHR+'.csv'
SAMPLEGENERATIONS=[1,5,10,50,100,200]
ALPHAS=[0.1, 0.2, 0.3, 0.4, 0.5]
POPS=['yoruba', 'french', 'han', 'bedouin']
NOFFSPRING=2
NGENS=5

#Parameters for Using and Storing SupportMix
C=1
WINSIZE=200
NGENS=5
OUTPUT_TWO_POP_SVM='data/twoPopResults.P'
OUTPUT_TWO_POP_LAMP='data/twoPopResultsLamp.P'
OUTPUT_THREE_ASIAN_SVM='data/threePopResultsHanJpn.P'
OUTPUT_THREE_AFRIC_SVM='data/threePopResultsCeuYri.P'
OUTPUT_TWO_POP_SVM_GENS='data/twoPopResultsGens.P'
OUTPUT_TWO_POP_SVM_ALPHA='data/twoPopResultsAlpha.P'
OUTPUT_TWO_POP_SVM_DELTA_GENS='data/twoPopResultsDeltaGens.P'
OUTPUT_TWO_POP_SVM_WIN='data/twoPopResultsWins.P'
OUTPUT_PCA='data/pcaResults.npz'

POPCOLORS={'Qatar1': [255,178,0,255], #Qatar 
'Qatar2':            [255,142,0,255], 
'Qatar3':            [255,106,0,255], 
'Palestinian':       [64,0,255,255],  #Arab like populations
'Bedouin':           [32,41,255,255],
'Druze':             [0,81,255,255],
'Makrani':           [255,0,0,255],  #Persian type
'Sindhi':            [255,0,12,255],
'Balochi':           [255,0,24,255],
'Brahui':            [255,0,36,255],
'Hazara':            [255,0,48,255],
'Pathan':            [255,0,60,255],
'Kalash':            [255,0,72,255],
'Burusho':           [255,0,84,255],  #N Eastern Asia
'Xibo':              [255,0,90,255], 
'Uygur':             [255,0,100,255], 
'Mozabite':          [0,183,255,255], #North Africa
'Mandenka':          [0,138,32,255], #SubSaharan Africa  
'Yoruba':            [0,167,56,255],                      
'Biaka Pygmies':     [0,196,80,255],                     
'Mbuti Pygmies':     [0,225,104,255],                     
'Bantu N.E.':        [0,255,128,255],                      
'Bantu S.W. Ovambo': [185,255,0,255], #S Africa
'San':               [199,255,0,255], 
'Bantu S.W. Herero': [213,255,0,255],
'Bantu S.E. S.Sotho':[227,255,0,255], 
'Bantu S.E. Tswana': [241,255,0,255],
'Bantu S.E. Zulu':   [255,255,0,255], 
'Yakut':             [255,0,136,255], #Asian
'Oroqen':            [253,0,144,255],
'Daur':              [251,0,152,255],
'Hezhen':            [249,0,160,255],
'Mongola':           [247,0,168,255],
'Japanese':          [245,0,176,255], 
'Tu':                [243,0,184,255],
'Han':               [241,0,192,255],
'Tujia':             [239,0,199,255],
'She':               [237,0,207,255],
'Miaozu':            [235,0,215,255],
'Yizu':              [233,0,223,255],
'Naxi':              [231,0,231,255],
'Lahu':              [229,0,239,255],
'Dai':               [227,0,247,255],
'Cambodians':        [225,0,255,255],
'Adygei':[0,0,0,255], 'Tuscan':[0,0,0,255], 'French':[0,0,0,255],  'French Basque':[0,0,0,255], 'Sardinian':[0,0,0,255], 'North Italian':[0,0,0,255], 'Orcadian':[0,0,0,255], 'Russian':[0,0,0,255],  #European
'Pima':[0,0,0,255],  'Maya':[0,0,0,255], 'Colombians':[0,0,0,255], 'Karitiana':[0,0,0,255], 'Surui':[0,0,0,255],  #Native American
'Papuan':[0,0,0,255], 'NAN Melanesian':[0,0,0,255]} #Oceania

for key,val in POPCOLORS.items():  #Fix colors to match fileNames
    POPCOLORS[key.lower().replace(' ', '_')]=val


class results(object):
    """Keeps track of multiple results. """
    def __init__(self, win, gens):
        """Assigns the results to variables
        Arguments:
        - `win`:       The window sized used to determine the admixture
        - `gens`:      Number of generations since admixture used for HMM model
        - `ancestral`: List of ancestral population names
        - `fst`:       List of fst values between all the populations
        - `success`:   List of mean success rates
        - `posterior`: List of matrices of subjects x positions of posterior probability
        - `ancestral`: List of matrices of subjects x positions of ancestral assignment
        """
        self.win = win
        self.gens = gens
        self.files = list()
        self.fst = list()
        self.success = list()
        self.posterior = list()
        self.ancestral = list()
        
    def append(self, files, fst, success, posterior, ancestral):
        """Appends the results of one supportMix run    
        Arguments:
        - `files`:
        - `fst`:
        - `success`:
        - `posterior`:
        - `ancestral`:  """
        self.files.append(files)
        self.fst.append(fst)
        self.success.append(success)
        self.posterior.append(posterior)
        self.ancestral.append(ancestral)
        
    def __str__(self):
        str='Window: %i\nGenerations:%i\n\n' %(self.win, self.gens)
        for i in range(len(self.files)):
            str+='%s\t%g\t%3.2g+/-%3.2g\n' %('-'.join(self.files[i]), self.fst[i], self.success[i][0], self.success[i][1])
        return str
