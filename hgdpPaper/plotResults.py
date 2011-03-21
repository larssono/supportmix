from variables import *
import cPickle, pylab, numpy as np
from mpl_toolkits.basemap import Basemap
import string

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
    # Figure 1: PCA and Maps
    ################################
    pylab.figure(figsize=(7.08, 7.08))
    ######### PCA plot ################
    ax=pylab.axes([0.09, 0.51, 0.45, 0.45])
    pca=np.load(OUTPUT_PCA)
    Vt=pca['Vt']; S=pca['S']; popLabels=pca['popLabels']; subjects=pca['subjects']
    colors=np.asarray([POPCOLORS[l] for l in popLabels ])/255.
    idxQatar1=(popLabels=='Qatar1')
    idxQatar2=(popLabels=='Qatar2')
    idxQatar3=(popLabels=='Qatar3')
    idxNonQatar=np.logical_not(idxQatar1+idxQatar2+idxQatar3)
    pylab.scatter(-Vt[0,idxNonQatar], -Vt[1,idxNonQatar], s=15, c=colors[idxNonQatar,:], linewidths=0)
    pylab.scatter(-Vt[0,idxQatar1], -Vt[1, idxQatar1], s=15, c=colors[idxQatar1,:], linewidths=0, marker='s')  #Qatar 1
    pylab.scatter(-Vt[0,idxQatar2], -Vt[1, idxQatar2], s=15, c=colors[idxQatar2,:], linewidths=0, marker='v')
    pylab.scatter(-Vt[0,idxQatar3], -Vt[1, idxQatar3], s=15, c=colors[idxQatar3,:], linewidths=0, marker='d')
    pylab.xlabel('PC 1 (%0.2g%%)' %S[0])
    pylab.ylabel('PC 2 (%0.2g%%)' %S[1])
    ax.annotate('a', xy=(-Vt[0,981], -Vt[1, 981] ),  xycoords='data', xytext=(-30, 10), textcoords='offset points',arrowprops=dict(facecolor='black', shrink=0.1, width=1, headwidth=4), horizontalalignment='center', verticalalignment='top' )
    ax.annotate('b', xy=(-Vt[0,925], -Vt[1, 925] ),  xycoords='data', xytext=(-15, -20), textcoords='offset points',arrowprops=dict(facecolor='black', shrink=0.1, width=1, headwidth=4), horizontalalignment='center', verticalalignment='top' )
    ax.annotate('c', xy=(-Vt[0,1029], -Vt[1, 1029] ),  xycoords='data', xytext=(-30, 10), textcoords='offset points',arrowprops=dict(facecolor='black', shrink=0.1, width=1, headwidth=4), horizontalalignment='center', verticalalignment='top' )
    ######### Map plot ################
    ax=pylab.axes([0.54, 0.51, 0.45, 0.45])
    map = Basemap(llcrnrlon=-180,llcrnrlat=-70,urcrnrlon=180,urcrnrlat=70, projection='merc',lat_ts=20, resolution = 'l',area_thresh = 100000. )
    map.drawcoastlines(linewidth=0.25)
    map.drawcountries(linewidth=0.05)
    map.fillcontinents(color='white',lake_color=[.9, .9, .9])
    #map.drawmapboundary(fill_color='aqua')
    #map.drawmeridians(np.arange(0,360,30))
    #map.drawparallels(np.arange(-90,90,30))
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
        admClass=np.load('data/qatarSupportMix/qatar.%(CHR)s.100.admixedClass.npy'%locals())
        p=np.load('data/qatarSupportMix/qatar.%(CHR)s.100.posterior.npy'%locals())
        pops=np.load('data/qatarSupportMix/qatar.%(CHR)s.100.populations.npy'%locals())
        subs=np.load('data/qatarSupportMix/qatar.%(CHR)s.100.subjects.npy'%locals())
        subs=np.asarray([sub[:-2] for sub in subs])
        idx=np.hstack([np.nonzero(subs==rightOrd)[0] for rightOrd in qatarOrder])  #Get index of sorted subjects
        colors=[POPCOLORS[label] for label in pops]
        vals=np.zeros((admClass.shape[0],admClass.shape[1],4))
        for i in range(admClass.shape[0]):
            for j in range(admClass.shape[1]):
                vals[i,j,:]=colors[admClass[i,j]]
        vals=vals/255.
        ax=pylab.axes([.07, startPos, .91, .8*admClass.shape[0]/602.])
        startPos+=.8*admClass.shape[0]/602.+hSpace
        pylab.imshow(vals[:,idx,:], interpolation='nearest')
        pylab.axis('tight'); pylab.draw(); 
        pylab.yticks([0,40], ['0', '250'], fontsize=6); 
        pylab.xticks([])
        pylab.ylabel(CHR+' pos[Mb]', fontsize=6)
    ########## Arrows to sample individuals ######
    ax=pylab.axes([.07, .03, .91, .06])
    pylab.axis([0, 312, 0, -10]); pylab.axis('off')
    pylab.text(156,0,  'Qatari individuals', horizontalalignment='center', fontsize=8)

    for number, pos in [('a',100),('b',210),('c',284)]:
        ax.annotate(number, xy=(pos, -9.5),  xycoords='data',
                    xytext=(pos, -5), textcoords='data',
                    arrowprops=dict(facecolor='black', shrink=0.1, width=1, headwidth=4),
                    horizontalalignment='center', verticalalignment='top')

    pylab.savefig('fig1.'+FILETYPE,format=FILETYPE) 

    ################################
    # Figure 2 - average ancestry
    ################################
    pylab.figure(figsize=(7.08,2.5))
    qatar=[]; oldSubs=[]; oldPops=[]
    for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,21,22]:
        CHR='chr%i' %i
        admClass=np.load('data/qatarSupportMix/qatar.%(CHR)s.100.admixedClass.npy'%locals())
        p=np.load('data/qatarSupportMix/qatar.%(CHR)s.100.posterior.npy'%locals())
        pops=np.load('data/qatarSupportMix/qatar.%(CHR)s.100.populations.npy'%locals())
        subs=np.load('data/qatarSupportMix/qatar.%(CHR)s.100.subjects.npy'%locals())
        subs=np.asarray([sub[:-2] for sub in subs])
        if oldSubs==[]:
            oldSubs=subs
            oldPops=pops
        if np.any(pops!=oldPops) or np.any(oldSubs!=subs):
            print 'Shiiite:  order of subjects or populations changed'
            break
        qatar.append(admClass)
    qatar=np.vstack(qatar)
    idx=np.hstack([np.nonzero(subs==rightOrd)[0] for rightOrd in qatarOrder])  #Get index of sorted subjects
    qatar=qatar[:, idx]
    q1=qatar[:,0:206]
    q2=qatar[:,206:278]
    q3=qatar[:,278:312]
    selectedPop=['bedouin','palestinian','druze', 'makrani', 'sindhi', 'balochi', 'brahui','hazara','pathan', 'kalash','burusho',  'mozabite','mandenka', 'yoruba','bantu_n.e.', 'Other']
    for i, qPop in enumerate([q1,q2,q3]):
        ax=pylab.subplot(1,3,i+1)
        totPercent=0
        for pos, pop in enumerate(selectedPop[:-1]):
            popNumber=pylab.find(pops==pop)
            percent=np.sum(qPop==popNumber)/float(qPop.size)*100
            totPercent+=percent
            pylab.bar(pos, percent, color=np.asarray(POPCOLORS[pop])/255.)
        pylab.bar(pos+1, 100-totPercent, color='k')
        pylab.ylim(0, 30)
        if i==0: 
            pylab.ylabel('Percent of Loci')
            ax.annotate('62.6%', xy=(.8, 30),  xycoords='data',
                        xytext=(4, 25), textcoords='data', arrowprops=dict(facecolor='black', shrink=0.1, width=1, headwidth=4),
                        horizontalalignment='center', verticalalignment='top', fontsize=8 )

        pylab.text(1, 65, ['a', 'b', 'c'][i])
        pylab.xticks(np.arange(len(selectedPop))+.45, selectedPop, rotation=90)
    pylab.subplots_adjust(left=.07, bottom=.28, right=.98, top=.95, hspace=.01)
    summary=[(pop, np.sum(qatar==i)/float(qatar.size)*100) for (i,pop) in enumerate(pops)]
    summary.sort(lambda s1,s2: int(np.sign(s2[1]-s1[1])))
    for s in summary:
        print '%15s\t%0.2g' %s
    for q in [q1,q2,q3]:
        print 'Subpopulation\n------------------'
        tot=np.zeros(q[1,:].shape)
        for pop in ['bedouin','palestinian','druze','mozabite']:
            popNumber=pylab.find(pops==pop)
            percent=np.sum(q==popNumber, 0)/7.1
            tot+=percent
            print '%15s\t%0.2g+/-%0.2g\t %0.2g-%0.2g' %(pop, np.mean(percent), np.std(percent), np.min(percent), np.max(percent))
        print '%5s\t%0.2g+/-%0.2g\t %0.2g-%0.2g' %('Arab', np.mean(tot), np.std(tot), np.min(tot), np.max(tot))
        tot=np.zeros(q[1,:].shape)
        for pop in ['makrani','sindhi','balochi','brahui','hazara','pathan','kalash','burusho']:
            popNumber=pylab.find(pops==pop)
            percent=np.sum(q==popNumber, 0)/7.1
            tot+=percent
            print '%15s\t%0.2g+/-%0.2g\t %0.2g-%0.2g' %(pop, np.mean(percent), np.std(percent), np.min(percent), np.max(percent))
        print '%5s\t%0.2g+/-%0.2g\t %0.2g-%0.2g' %('Persian', np.mean(tot), np.std(tot), np.min(tot), np.max(tot))
        tot=np.zeros(q[1,:].shape)
        for pop in ['french','russian','north_italian','french_basque','adygei','orcadian','sardinian','tuscan']:
            popNumber=pylab.find(pops==pop)
            percent=np.sum(q==popNumber, 0)/7.1
            tot+=percent
        print '%5s\t%0.2g+/-%0.2g\t %0.2g-%0.2g' %('European', np.mean(tot), np.std(tot), np.min(tot), np.max(tot))
        tot=np.zeros(q[1,:].shape)
        for pop in ['mandenka','yoruba','biaka_pygmies','mbuti_pygmies','bantu_n.e.']:
            popNumber=pylab.find(pops==pop)
            percent=np.sum(q==popNumber, 0)/7.1
            tot+=percent
            print '%15s\t%0.2g+/-%0.2g\t %0.2g-%0.2g' %(pop, np.mean(percent), np.std(percent), np.min(percent), np.max(percent))
        print '%5s\t%0.2g+/-%0.2g\t %0.2g-%0.2g' %('African', np.mean(tot), np.std(tot), np.min(tot), np.max(tot))
        tot=np.zeros(q[1,:].shape)
        for pop in ['bantu_s.e._s.sotho','bantu_s.e._tswana','bantu_s.e._zulu', 'bantu_s.w._herero','bantu_s.w._ovambo', 'san',]:
            popNumber=pylab.find(pops==pop)
            percent=np.sum(q==popNumber, 0)/7.1
            tot+=percent
        print '%5s\t%0.2g+/-%0.2g\t %0.2g-%0.2g' %('S. African', np.mean(tot), np.std(tot), np.min(tot), np.max(tot))
    pylab.savefig('fig2.'+FILETYPE,format=FILETYPE) 

    ################################
    # Figure 3: Simulation results
    ################################
    pylab.figure(figsize=(7.08,4.5))
    ######### HGDP SupportMix v Lamp ################
    pylab.subplot(2,2,1)
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
    ######### HGDP Generations ################
    pylab.subplot(2,2,2)    
    with open(OUTPUT_TWO_POP_SVM_GENS,'r') as fp: gens=cPickle.load(fp)
    success=np.asarray(gens.success)
    pylab.semilogx(gens.fst[:6], success.reshape((6,6,2))[:,:,0].T, '-o', linewidth=.5, markersize=5)
    pylab.xticks(gens.fst[:6], gens.fst[:6])
    pylab.xlim(-4,210)
    pylab.xlabel('Generations since admixture')
    ######### HGDP Wrong g ################
    pylab.subplot(2,2,3)
    with open(OUTPUT_TWO_POP_SVM_DELTA_GENS,'r') as fp: g=cPickle.load(fp)
    success=np.asarray(g.success)
    for i in range(6):
        pylab.semilogx(g.fst[:7], success[i*7:(i+1)*7,0], '-o', linewidth=.5, markersize=5)
    pylab.xticks(g.fst[:7], g.fst[:7])
    pylab.xlabel('g/g\'')
    pylab.ylabel('Correctly classified loci [%]')
    pylab.axis([0.045, 21, 50, 100])
    ######### HGDP WinSize ################
    pylab.subplot(2,2,4)
    with open(OUTPUT_TWO_POP_SVM_WIN,'r') as fp: g=cPickle.load(fp)
    success=np.asarray(g.success)
    for i in range(6):
        pylab.semilogx(g.fst[:8], success[i*8:(i+1)*8,0], '-o', linewidth=.5, markersize=5)
    pylab.xticks(g.fst[:8], g.fst[:8])
    pylab.xlabel('Window size')
    pylab.axis([9, 5050, 50, 100])
    pylab.legend(['-'.join(g.files[i]) for i in [0,8,16,24,32,40]], 4, ncol=2)
    pylab.subplots_adjust(left=.078, bottom=.08, right=.97, top=.97, hspace=.18, wspace=.1)
    pylab.savefig('fig3.'+FILETYPE,format=FILETYPE) 

    ################################
    # Supplemental Figure 1 - Structure+all chroms
    ################################
    pylab.figure(figsize=(8.8, 10.8))
    ####### STRUCTURE PLOT ############x##
    pylab.axes([.1,.85, .8, .1])
    names=np.asarray([l.strip().split() for l in open('data/qatar_admix_q')])
    alphas=names[:,1:].astype(np.float)
    names=names[:,0]
    idx=[list(names).index(rightOrd) for rightOrd in qatarOrder]
    pylab.bar(np.arange(156), alphas[idx,0], width=1, color='b', linewidth=0)
    pylab.bar(np.arange(156), alphas[idx,1], width=1, bottom=alphas[idx,0], linewidth=0,color='g')
    pylab.bar(np.arange(156), alphas[idx,2], width=1, bottom=alphas[idx,:2].sum(1), linewidth=0, color='r')
    pylab.axis([-0.2, 156.2, 0, 1]); 
    pylab.yticks(np.linspace(0,1,6)), pylab.xticks([])
    pylab.ylabel('STRUCTURE ancestry', rotation='vertical', fontsize=7)
    ####### SupportMix Plots ############x##    
    startPos=.05
    hSpace=0.001
    for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,21,22]:
        CHR='chr%i' %i
        admClass=np.load('data/qatarSupportMix/qatar.%(CHR)s.100.admixedClass.npy'%locals())
        p=np.load('data/qatarSupportMix/qatar.%(CHR)s.100.posterior.npy'%locals())
        pops=np.load('data/qatarSupportMix/qatar.%(CHR)s.100.populations.npy'%locals())
        subs=np.load('data/qatarSupportMix/qatar.%(CHR)s.100.subjects.npy'%locals())
        subs=np.asarray([sub[:-2] for sub in subs])
        idx=np.hstack([np.nonzero(subs==rightOrd)[0] for rightOrd in qatarOrder])  #Get index of sorted subjects
        colors=[POPCOLORS[label] for label in pops]
        vals=np.zeros((admClass.shape[0],admClass.shape[1],4))
        for i in range(admClass.shape[0]):
            for j in range(admClass.shape[1]):
                vals[i,j,:]=colors[admClass[i,j]]
        vals=vals/255.
        #print startPos
        pylab.axes([.1, startPos, .8, .77*admClass.shape[0]/710.])
        startPos+=.77*admClass.shape[0]/710.+hSpace
        vals[:,:,3]=p**.2
        pylab.imshow(vals[:,idx,:], interpolation='nearest')
        pylab.axis('tight'); pylab.draw(); 
        pylab.xticks([]); pylab.yticks([])
        pylab.ylabel(CHR, rotation='horizontal')
        if CHR=='chr1':
            pylab.xlabel('Qatari individuals')
    pylab.text(-28, 450, 'SupportMix ancestry assignment across chromosomes', rotation='vertical', fontsize=8)
    pylab.savefig('supplemental_fig1.'+FILETYPE,format=FILETYPE) 


    # ################################
    # # Supplemental figure 2
    # ################################
    pylab.figure(figsize=(7.08,4.5))
    ######### HGDP SupportMix three ################
    pylab.subplot(1,2,1)
    svm3=[]
    with open(OUTPUT_THREE_AFRIC_SVM, 'r') as fp:  svm3.append(cPickle.load(fp))
    with open(OUTPUT_THREE_ASIAN_SVM, 'r') as fp:  svm3.append(cPickle.load(fp))
    pylab.plot(svm3[0].fst, np.asarray(svm3[0].success)[:,0], '.g')
    pylab.plot(svm3[1].fst, np.asarray(svm3[1].success)[:,0], '.r')
    pylab.axis([0, 0.10, 50, 100])
    pylab.xlabel('Fst')
    pylab.ylabel('Correctly classified loci [%]')
    ######### HGDP SupportMix alpha ################
    pylab.subplot(1,2,2)
    with open(OUTPUT_TWO_POP_SVM_ALPHA,'r') as fp: alphas=cPickle.load(fp)
    success=np.asarray(alphas.success)
    pylab.plot(alphas.fst[:5], success.reshape((6,5,2))[:,:,0].T, '-o')    
    for i in range(6):
        pylab.errorbar(alphas.fst[:5], success[i*5:(i+1)*5,0], success[i*5:(i+1)*5,1], fmt=None, ecolor='k')
    pylab.xticks(alphas.fst[:5], alphas.fst[:5])
    pylab.legend(['-'.join(alphas.files[i]) for i in [0,5,10,15,20,25]], 3)
    pylab.ylim(50, 100)
    pylab.xlabel('Ancestry fraction')
    pylab.xlim(0.08, .52)
    pylab.subplots_adjust(left=.078, bottom=.09, right=.97, top=.97, wspace=.1)
    pylab.savefig('supplemental_fig2.'+FILETYPE,format=FILETYPE) 
    ######### Effects of Phasing ###################
    # pylab.subplot(1,3,3)
    # width=.8/3
    # pos=np.arange(3)
    # success=np.asarray([[99,95,90], 
    #                     [99,98, 89], 
    #                     [95,90, 87]])
    # h1=pylab.bar(pos, success[:,0], width=width)
    # h2=pylab.bar(pos+width, success[:,1], width=width, color='g')
    # h3=pylab.bar(pos+2*width, success[:,2], width=width, color='r')
    # pylab.xticks(pos+1.5*width, ['YRI-CEU', 'CEU-HAN', 'CEU-HAN'])
    # pylab.axis([-.25, 3.15, 50, 100]) 
    # pylab.legend([h1[0],h2[0],h3[0]], ['Trio', 'Beagle','fastPhase'], 4)
    # pylab.text(.8, 75, 'CARTOON GRAPHICS',bbox=dict(facecolor='white', alpha=0.9))
    # pylab.legend(['-'.join(alphas.files[i]) for i in [0,5,10,15,20,25]], 3)

    



    




    pylab.show()
