from variables import *
import cPickle, pylab, numpy as np
from mpl_toolkits.basemap import Basemap
import string, sys, fileReader

FILETYPE='eps'

popLegend=['Palestinian','Bedouin','Druze','Makrani','Sindhi','Balochi','Brahui','Hazara','Pathan','Kalash','Burusho','Mozabite','Mandenka','Yoruba','Biaka Pygmies','Mbuti Pygmies',
'Bantu N.E.','Bantu S.W. Ovambo','San','Bantu S.W. Herero','Bantu S.E. S.Sotho', 'Bantu S.E. Tswana','Bantu S.E. Zulu','Yakut',
'Oroqen','Daur','Hezhen','Mongola','Japanese','Tu','Han','Tujia','She','Miaozu','Yizu','Naxi','Lahu','Dai','Cambodians', 'Xibo', 'Uygur']


if __name__ == '__main__':
    params = {'axes.labelsize': 8,
              'text.fontsize': 8,
              'legend.fontsize': 7,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'text.usetex': False}
    pylab.rcParams.update(params)
    #Sort order of Qatari's to match structure results from Haley's paper?
    qatarOrder=[l.strip().split()[0] for l in open('data/CrystalQatar/haley_keepfile.txt')]

    ################################
    #Read all qatar SupportMix results
    ################################
    qatarPops=np.load('data/qatarSupportMix/qatar.100.populations.npy')
    colors=[POPCOLORS[label] for label in qatarPops]
    subs=np.load('data/qatarSupportMix/qatar.100.subjects.npy')
    subs=np.asarray([sub[:-2] for sub in subs])
    idx=np.hstack([np.nonzero(subs==rightOrd)[0] for rightOrd in qatarOrder])  #Get index of sorted subjects
    qatarAncestry=np.loadtxt('data/qatarSupportMix/qatar.100.admixedClass.csv')[:,idx]
    qatarAncestryP=np.loadtxt('data/qatarSupportMix/qatar.100.posterior.csv')[:,idx]
    qatarPositions=np.loadtxt('data/qatarSupportMix/qatar.100.position.csv', dtype=np.str)
    qatarColors=np.zeros((qatarAncestry.shape[0],qatarAncestry.shape[1],4))
    nQatarWins, nQatarSubs=qatarAncestry.shape
    for i in range(nQatarWins):
        for j in range(nQatarSubs):
            qatarColors[i,j,:]=colors[int(qatarAncestry[i,j])]
    qatarColors=qatarColors/255.

    ################################
    # Figure 1: PCA and Maps
    ################################
    pylab.figure(figsize=(7.08, 7.03))
    ######### PCA plot ################
    ax=pylab.axes([0.093, 0.51, 0.45, 0.45])
    pca=np.load(OUTPUT_PCA)
    Vt=pca['Vt']; S=pca['S']; popLabels=pca['popLabels']; subjects=pca['subjects']
    colors=np.asarray([POPCOLORS[l] for l in popLabels ])/255.
    idxQatar1=(popLabels=='Qatar1')
    idxQatar2=(popLabels=='Qatar2')
    idxQatar3=(popLabels=='Qatar3')
    idxNonQatar=np.logical_not(idxQatar1+idxQatar2+idxQatar3)
    pylab.scatter(-Vt[0,idxNonQatar], -Vt[1,idxNonQatar], s=15, c=colors[idxNonQatar,:], linewidths=0)
    pylab.scatter(-Vt[0,idxQatar1], -Vt[1, idxQatar1], s=15, c=colors[idxQatar1,:], linewidths=.1, marker='s')  #Qatar 1
    pylab.scatter(-Vt[0,idxQatar2], -Vt[1, idxQatar2], s=15, c=colors[idxQatar2,:], linewidths=.1, marker='v')
    pylab.scatter(-Vt[0,idxQatar3], -Vt[1, idxQatar3], s=15, c=colors[idxQatar3,:], linewidths=.1, marker='d')
    pylab.xlabel('Principal component 1 (%0.2g%%)' %S[0])
    pylab.ylabel('Principal component 2 (%0.2g%%)' %S[1])
    ax.annotate('q1', xy=(-Vt[0,981], -Vt[1, 981] ),  xycoords='data', xytext=(-30, 10), textcoords='offset points',arrowprops=dict(facecolor='black', shrink=0.1, width=1, headwidth=4), horizontalalignment='center', verticalalignment='top' )
    ax.annotate('q2', xy=(-Vt[0,925], -Vt[1, 925] ),  xycoords='data', xytext=(-15, -20), textcoords='offset points',arrowprops=dict(facecolor='black', shrink=0.1, width=1, headwidth=4), horizontalalignment='center', verticalalignment='top' )
    ax.annotate('q3', xy=(-Vt[0,1029], -Vt[1, 1029] ),  xycoords='data', xytext=(-30, 10), textcoords='offset points',arrowprops=dict(facecolor='black', shrink=0.1, width=1, headwidth=4), horizontalalignment='center', verticalalignment='top' )
    pylab.text(0.03, 1.02, 'A', transform = ax.transAxes, horizontalalignment='right', fontsize=12)
    ######### Map plot ################
    ax=pylab.axes([0.54, 0.51, 0.45, 0.45])
    map = Basemap(llcrnrlon=-180,llcrnrlat=-70,urcrnrlon=180,urcrnrlat=70, projection='merc',lat_ts=20, resolution = 'l',area_thresh = 100000. )
    map.drawcoastlines(linewidth=0.25)
    map.drawcountries(linewidth=0.05)
    map.fillcontinents(color='white',lake_color=[.9, .9, .9])
    pops=np.asarray([l.strip().split('\t')[:4] for l in open(HGDPCOORDSFILE).readlines()[1:]])
    lats=np.asarray(pops[:,1], np.float); lons=np.asarray(pops[:,3], np.float)
    pops=pops[:,0]; 
    # compute the native map projection coordinates for cities.
    xc,yc = map(lons,lats)
    colors=np.asarray([POPCOLORS[l] for l in pops])/255.
    # plot filled circles at the locations of the cities.
    pylab.scatter(xc,yc,s=40, zorder=10, c=colors, linewidth=0)
    xQatar, yQatar=map(51.53333, 25.28667)
    pylab.scatter(xQatar,yQatar,s=40, zorder=10, c=np.asarray([POPCOLORS['Qatar1']])/255., linewidth=0)
    pylab.arrow(xQatar, yQatar, 3e6,-3e6, zorder=11, lw=.5)
    pylab.text(xQatar+3e6, yQatar-3.1e6, 'Qatar', horizontalalignment='center', verticalalignment='top', fontsize=8)
    # for name,xpt,ypt in zip(pops,xc,yc):
    #     pylab.text(xpt+50000,ypt+50000,name,fontsize=9)
    pylab.axis([17e6, 30e6, 6.5e6, 20e6])
    pylab.axis('off')
    pylab.subplots_adjust(left=.08, bottom=.02, right=.98, top=.98, hspace=.01)
    pylab.text(0.03, 1.02, 'B', transform = ax.transAxes, horizontalalignment='right', fontsize=12)
    ############# Legend ###################
    pylab.axes([.07, 0.34, .91, .11])
    x,y=0,0
    X=(x for x in [5.8, 11.4, 18.4, 25.2, 31.4, 37.2])
    for pop in popLegend:
        pylab.scatter(x,y, s=20, c=np.asarray(POPCOLORS[pop])/255., linewidth=0)
        pylab.text(x+.5, y, pop, fontsize=6, verticalalignment='center')
        y+=5;
        if y==30:
            x=X.next(); y=0
    pylab.scatter(x,y, s=20, c=[0,0,0], linewidth=0)
    pylab.text(x+.5, y, 'Other (Europe/Americas/Oceania)', fontsize=6, verticalalignment='center')
    pylab.xlim(-.25, 48)
    pylab.ylim(27, -2)
    pylab.axis('off')
    ######### SupportMix Plots ################
    startPos=0.09; hSpace=0.01
    for i in range(3,0, -1):
        CHR='chr%i' %i
        chrIdx=qatarPositions[:,0]==CHR
        height=sum(chrIdx)/float(nQatarWins)
        ax=pylab.axes([.07, startPos, .91, .8*height])
        startPos+=.8*height+hSpace
        pylab.imshow(qatarColors[chrIdx,:,:], interpolation='nearest')
        pylab.axis('tight'); pylab.draw(); 
        pylab.xticks([]); pylab.yticks([])
        pylab.ylabel(CHR, fontsize=6, rotation='horizontal', horizontalalignment='right')
    pylab.text(-20, -20, 'Positions across chromosomes', rotation='vertical', fontsize=8, verticalalignment='top')
    pylab.text(0.0, 1.03, 'C', transform = ax.transAxes, horizontalalignment='right', fontsize=12)
    ########## Arrows to sample individuals ######
    ax=pylab.axes([.07, .03, .91, .06])
    pylab.axis([0, 312, 0, -10]); pylab.axis('off')
    pylab.text(156,0,  'Qatari individuals', horizontalalignment='center', fontsize=8)
    for number, pos in [('q1',100),('q2',210),('q3',284)]:
        ax.annotate(number, xy=(pos, -9.5),  xycoords='data',
                    xytext=(pos, -5), textcoords='data',
                    arrowprops=dict(facecolor='black', shrink=0.1, width=1, headwidth=4),
                    horizontalalignment='center', verticalalignment='top')
    pylab.savefig('fig1.'+FILETYPE,format=FILETYPE) 

    ################################
    # Figure 2 - average ancestry
    ################################
    pylab.figure(figsize=(7.08,2.5))
    q1=qatarAncestry[:,0:206]
    q2=qatarAncestry[:,206:278]
    q3=qatarAncestry[:,278:312]
    selectedPop=['bedouin','palestinian','druze', 'makrani', 'sindhi', 'balochi', 'brahui','hazara','pathan', 'kalash','burusho',  'mozabite','mandenka', 'yoruba','bantu_n.e.', 'Other']
    selectedLab=['Bedouin','Palestinian','Druze', 'Makrani', 'Sindhi', 'Balochi', 'Brahui','Hazara','Pathan', 'Kalash','Burusho',  'Mozabite','Mandenka', 'Yoruba','Bantu N.E.', 'Other']
    for i, qPop in enumerate([q1,q2,q3]):
        ax=pylab.subplot(1,3,i+1)
        totPercent=0
        for pos, pop in enumerate(selectedPop[:-1]):
            popNumber=pylab.find(qatarPops==pop)
            percent=np.sum(qPop==popNumber)/float(qPop.size)*100
            totPercent+=percent
            pylab.bar(pos, percent, color=np.asarray(POPCOLORS[pop])/255.)
        pylab.bar(pos+1, 100-totPercent, color='k')
        pylab.ylim(0, 30); pylab.xlim(-0.2, 16.2)
        if i==0: 
            pylab.ylabel('Percent of Loci')
            ax.annotate('63%', xy=(.8, 30),  xycoords='data',
                        xytext=(4, 25), textcoords='data', arrowprops=dict(facecolor='black', shrink=0.1, width=1, headwidth=4),
                        horizontalalignment='center', verticalalignment='top', fontsize=8 )
        pylab.xticks(np.arange(len(selectedPop))+.45, selectedLab, rotation=90)
        pylab.text(0.03, 1.02, ['A', 'B', 'C'][i], transform = ax.transAxes, horizontalalignment='right', fontsize=12)
    pylab.subplots_adjust(left=.07, bottom=.28, right=.98, top=.9, hspace=.01)
    pylab.savefig('fig2.'+FILETYPE,format=FILETYPE) 
    ################################
    # Summary Output
    ################################
    summary=[(pop, np.sum(qatarAncestry==i)/float(qatarAncestry.size)*100) for (i,pop) in enumerate(qatarPops)]
    summary.sort(lambda s1,s2: int(np.sign(s2[1]-s1[1])))
    #for s in summary:
    #    print '%15s\t%0.2g' %s
    fp=open('supplemental_table.tex', 'w')
    populationGroups={'Middle East': ['bedouin','palestinian','druze','mozabite', 'All'],
                      'Persian': ['makrani','sindhi','balochi','brahui','hazara','pathan','kalash','burusho', 'All'],
                      'European': ['french','russian','north_italian','french_basque','adygei','orcadian','sardinian','tuscan', 'All'],
                      'sub-Saharan Africa': ['mandenka','yoruba','biaka_pygmies','mbuti_pygmies','bantu_n.e.', 'All'],
                      'South African': ['bantu_s.e._s.sotho','bantu_s.e._tswana','bantu_s.e._zulu', 'bantu_s.w._herero','bantu_s.w._ovambo', 'san', 'All'],
                      'Asian': ['cambodians', 'dai', 'daur', 'han',  'hezhen', 'japanese','lahu', 'miaozu', 'mongola',  'naxi', 'oroqen', 'she','tu', 'tujia', 'uygur', 'xibo', 'yakut', 'yizu', 'All'],
                      'Other': ['colombians', 'karitiana', 'maya', 'nan_melanesian', 'papuan', 'pima',  'surui', 'All']}
    for key, popNames in populationGroups.items():
        fp.write(key)
        tot=[np.zeros(q1.shape[1]),np.zeros(q2.shape[1]),np.zeros(q3.shape[1])]
        for pop in popNames:
            show=False
            str=' \t&%s\t' %pop.replace('_', ' ').capitalize()
            for i, q in enumerate([q1,q2,q3]):
                if pop=='All':
                    percent=tot[i]
                    show=True
                else:
                    popNumber=pylab.find(qatarPops==pop)
                    percent=np.sum(q==popNumber, 0)/(nQatarWins/100.)
                    tot[i]+=percent
                if np.mean(percent)>0.1:
                    show=True
                str+='& $%3.2g\pm %3.2g$ \t' %(np.mean(percent), np.std(percent))
                #fp.write('%0.2g+/-%0.2g\t %0.2g-%0.2g' %(np.mean(percent), np.std(percent), np.min(percent), np.max(percent)))
            if show:
                if pop=='All':
                    fp.write( '\cline{2-5} '+ str+'\\\\ \hline \hline \n')
                else:
                    fp.write(str+'\\\\ \n ')
        

    ################################
    # Figure 3: Simulation results
    ################################
    pylab.figure(figsize=(7.08,4.5))
    ######### HGDP SupportMix v Lamp ################
    ax=pylab.subplot(2,2,1)
    with open(OUTPUT_TWO_POP_SVM, 'r') as fp:  svm2=cPickle.load(fp)
    with open(OUTPUT_TWO_POP_LAMP, 'r') as fp: lamp=cPickle.load(fp)
    fst=np.asarray(svm2.fst)
    success=np.asarray(svm2.success)[:,0]
    pylab.plot(fst, success, '.b')
    pylab.plot(lamp.fst, np.asarray(lamp.success)[:,0], 'ro')
    for files in lamp.files:
        i=pylab.find([files==s for s in svm2.files])
        pylab.scatter(fst[i], success[i], s=30, facecolors='none', edgecolors='r', linewidths=1, zorder=10)
    pylab.axis([0, 0.25, 50, 100])
    pylab.xlabel('Fst')
    pylab.ylabel('Correctly classified loci [%]')
    pylab.text(0.03, 1.02, 'A', transform = ax.transAxes, horizontalalignment='right', fontsize=12)
    ######### HGDP Generations ################
    ax=pylab.subplot(2,2,2)    
    with open(OUTPUT_TWO_POP_SVM_GENS,'r') as fp: gens=cPickle.load(fp)
    success=np.asarray(gens.success)
    pylab.semilogx(gens.fst[:6], success.reshape((7,6,2))[:,:,0].T, '-o', linewidth=.5, markersize=5)
    pylab.xticks(gens.fst[:6], gens.fst[:6]); pylab.yticks(range(50,101,10), [])
    pylab.xlim(-4,210)
    pylab.xlabel('Generations since admixture')
    pylab.text(0.03, 1.02, 'B', transform = ax.transAxes, horizontalalignment='right', fontsize=12)
    ######### HGDP Wrong g ################
    ax=pylab.subplot(2,2,3)
    with open(OUTPUT_TWO_POP_SVM_DELTA_GENS,'r') as fp: g=cPickle.load(fp)
    success=np.asarray(g.success)
    for i in range(7):
        pylab.semilogx(g.fst[:7], success[i*7:(i+1)*7,0], '-o', linewidth=.5, markersize=5)
    pylab.xticks(g.fst[:7], g.fst[:7])
    pylab.xlabel('g/g\'')
    pylab.ylabel('Correctly classified loci [%]')
    pylab.axis([0.045, 21, 50, 100])
    pylab.text(0.03, 1.02, 'C', transform = ax.transAxes, horizontalalignment='right', fontsize=12)
    ######### HGDP WinSize ################
    ax=pylab.subplot(2,2,4)
    with open(OUTPUT_TWO_POP_SVM_WIN,'r') as fp: g=cPickle.load(fp)
    success=np.asarray(g.success)
    for i in range(7):
        pylab.semilogx(g.fst[:8], success[i*8:(i+1)*8,0], '-o', linewidth=.5, markersize=5)
    pylab.xticks(g.fst[:8], g.fst[:8]); pylab.yticks(range(50,101,10), [])
    pylab.xlabel('Window size [# SNPs]')
    pylab.axis([9, 5050, 50, 100])
    pylab.legend(['-'.join([l.capitalize() for l in g.files[i]]) for i in [0,8,16,24,32,40,48]], 4, ncol=2)
    pylab.subplots_adjust(left=.078, bottom=.09, right=.97, top=.93, hspace=.19, wspace=.1)
    pylab.text(0.03, 1.02, 'D', transform = ax.transAxes, horizontalalignment='right', fontsize=12)
    pylab.savefig('fig3.'+FILETYPE,format=FILETYPE) 

    ################################
    # Supplemental Figure 1 - Structure+all chroms
    ################################
    pylab.figure(figsize=(8.8, 10.8))
    ####### SupportMix Plots ############x##    
    startPos=.95
    hSpace=0.001
    for i in range(1,23):
        CHR='chr%i' %i
        chrIdx=qatarPositions[:,0]==CHR
        height=.77*sum(chrIdx)/nQatarWins
        pylab.axes([.1, startPos-height, .8, height])
        startPos-=height+hSpace
        pylab.imshow(qatarColors[chrIdx,:,:], interpolation='nearest')
        pylab.axis('tight'); pylab.draw(); 
        pylab.xticks([]); pylab.yticks([])
        pylab.ylabel(CHR, rotation='horizontal')
    pylab.text(-25, -150, 'SupportMix locus-specific ancestry assignments across chromosomes', rotation='vertical', fontsize=10)

    ####### STRUCTURE PLOT ############x##
    ax=pylab.axes([.1,.05, .8, .095])
    names=np.asarray([l.strip().split() for l in open('data/qatar_admix_q')])
    alphas=names[:,1:].astype(np.float)
    names=names[:,0]
    idx=[list(names).index(rightOrd) for rightOrd in qatarOrder]
    pylab.bar(np.arange(156), alphas[idx,0], width=1, color='b', linewidth=0)
    pylab.bar(np.arange(156), alphas[idx,1], width=1, bottom=alphas[idx,0], linewidth=0,color='g')
    pylab.bar(np.arange(156), alphas[idx,2], width=1, bottom=alphas[idx,:2].sum(1), linewidth=0, color='r')
    pylab.axis([-0.2, 156.2, 0, 1]); 
    pylab.yticks(np.linspace(0,1,6)), pylab.xticks([])
    pylab.text(-12.5, -.1, 'STRUCTURE ancestry', rotation='vertical', fontsize=8)
    pylab.xticks([-.2, 51.5, 103.01, 119.5, 139, 147.5,  156], ['|','Qatar 1','|', 'Qatar 2', '|', 'Qatar 3', '|'])
    pylab.xlabel('Individual genomes', fontsize=8)

    pylab.savefig('supplemental_fig1.'+FILETYPE,format=FILETYPE) 


    ################################
    # Supplemental figure 2
    ################################
    pylab.figure(figsize=(7.08,3.8))
    ######### HGDP SupportMix three ################
    ax=pylab.subplot(1,2,1)
    svm3=[]
    with open(OUTPUT_THREE_AFRIC_SVM, 'r') as fp:  svm3.append(cPickle.load(fp))
    with open(OUTPUT_THREE_ASIAN_SVM, 'r') as fp:  svm3.append(cPickle.load(fp))
    pylab.plot(svm3[0].fst, np.asarray(svm3[0].success)[:,0], '.g')
    pylab.plot(svm3[1].fst, np.asarray(svm3[1].success)[:,0], '.r')
    pylab.axis([0, 0.08, 50, 100])
    pylab.xlabel('Fst')
    pylab.ylabel('Correctly classified loci [%]')
    pylab.legend(['Yoruba-French-other', 'Yoruba-Bedouin-other'], 4)
    pylab.text(0.0, 1.02, 'A', transform = ax.transAxes, horizontalalignment='right', fontsize=12)
    ######### HGDP SupportMix alpha ################
    ax=pylab.subplot(1,2,2)
    with open(OUTPUT_TWO_POP_SVM_ALPHA,'r') as fp: alphas=cPickle.load(fp)
    success=np.asarray(alphas.success)
    pylab.plot(alphas.fst[:5], success.reshape((7,5,2))[:,:,0].T, '-o')    
    for i in range(6):
        pylab.errorbar(alphas.fst[:5], success[i*5:(i+1)*5,0], success[i*5:(i+1)*5,1], fmt=None, ecolor='k')
    pylab.xticks(alphas.fst[:5], alphas.fst[:5])
    pylab.legend(['-'.join(alphas.files[i]) for i in [0,5,10,15,20,25,30]], 3)
    pylab.ylim(50, 100); pylab.yticks(range(50,101, 10), [])
    pylab.xlabel('Ancestry fraction')
    pylab.xlim(0.08, .52)
    pylab.text(0.0, 1.02, 'B', transform = ax.transAxes, horizontalalignment='right', fontsize=12)
    pylab.subplots_adjust(left=.078, bottom=.1, right=.97, top=.92, wspace=.1)
    pylab.savefig('supplemental_fig2.'+FILETYPE,format=FILETYPE) 
    


    ################################
    # Supplemental Figure 3 - Simulated Qatari
    ################################
    simQatarPops=np.load('data/simulatedQatar.populations.npy')
    colors=[POPCOLORS[label] for label in simQatarPops]
    simQatarAncestry=np.loadtxt('data/simulatedQatar.admixedClass.csv')
    simQatarAncestryP=np.loadtxt('data/simulatedQatar.posterior.csv')
    nsimQatarWins, nsimQatarSubs=simQatarAncestry.shape
    simQatarColors=np.zeros((nsimQatarSubs,nsimQatarWins,4))
    for i in range(nsimQatarWins):
        for j in range(nsimQatarSubs):
            simQatarColors[j,i,:]=colors[int(simQatarAncestry[i,j])]
    simQatarColors=simQatarColors/255.


    simQatarCorrect=np.asarray([l.strip().split('\t')[2:] for l in fileReader.openfile('data/hgdp3/admixed_hgdp_origin_yoruba_bedouin_brahui.chr1.csv.gz').readlines()[1:]], np.int)
    simQatarCorrect[simQatarCorrect==0]=np.nonzero(simQatarPops==['yoruba'])[0][0]
    simQatarCorrect[simQatarCorrect==1]=np.nonzero(simQatarPops==['bedouin'])[0][0]
    simQatarCorrect[simQatarCorrect==2]=np.nonzero(simQatarPops==['brahui'])[0][0]
    nsimQatarWins, nsimQatarSubs=simQatarCorrect.shape
    simQatarCorrectColors=np.zeros((nsimQatarSubs,nsimQatarWins,4))
    for i in range(nsimQatarWins):
        for j in range(nsimQatarSubs):
            simQatarCorrectColors[j,i,:]=colors[int(simQatarCorrect[i,j])]
    simQatarCorrectColors=simQatarCorrectColors/255.

    comparison=np.repeat(simQatarAncestry, 400, 0)[:nsimQatarWins,:]
    success=100*(comparison==simQatarCorrect).sum(0)/float(nsimQatarWins)
    print 'Correct %0.3g +/- %0.2g' %(np.mean(success), np.std(success))

    pylab.figure(figsize=(8.8, 2.8))
    pylab.axes([.1, .52, .75, .38 ])
    pylab.imshow(simQatarColors, interpolation='nearest')
    pylab.axis('tight'); pylab.draw(); 
    pylab.xticks([]); pylab.yticks([]);
    pylab.ylabel('Estimated ancestry')
    pylab.axes([.1, .1, .75, .38 ])
    pylab.imshow(simQatarCorrectColors, interpolation='nearest')
    pylab.axis('tight'); pylab.draw(); 
    pylab.xticks([]); pylab.yticks([])
    pylab.ylabel('Correct ancestry');pylab.xlabel('Position along chr1')
    pylab.text(-2100, 2.5, 'Four haploid genomes', horizontalalignment='center', rotation='vertical')
    pylab.axes([.86, 0.1, .15, .8])
    x,y=0,0
    popLegend=['Adygei', 'French', 'Yoruba', 'Bantu N.E.','Mbuti Pygmies','Biaka Pygmies', 'Mandenka', 
               'Brahui', 'Pathan','Burusho', 'Hazara', 'Makrani','Mozabite','Druze','Palestinian', 
               'Bedouin']
    for pop in popLegend:
        pylab.scatter(x,y, s=20, c=np.asarray(POPCOLORS[pop])/255., linewidth=0)
        pylab.text(x+.5, y, pop, fontsize=8, verticalalignment='center')
        y+=1;
    pylab.axis([-0.5, 8, -.3, 15.3])
    pylab.axis('off')



    pylab.savefig('supplemental_fig3.'+FILETYPE,format=FILETYPE) 

    pylab.show()
