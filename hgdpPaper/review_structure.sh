#!/bin/bash

#Generate data
python ../../../simulateAdmixture.py -c 1 -sceu_chr1,yri_chr1,admixed_ceu_yri_chr1 -n4 -g5 -f beagle --geneticMap ../../../../../human_genome_data/hapmap2/genetic_map_chr1_b36.txt ../../../../../human_genome_data/HumanPhasedData/hgdp/french/hgdp_french.chr1.bgl.phased.gz ../../../../../human_genome_data/HumanPhasedData/hgdp/yoruba/hgdp_yoruba.chr1.bgl.phased.gz 
sed '0~2d' admixed_ceu_yri_chr1_origin.tped|sed '0~2d' >admixed_ceu_yri_chr1_short_origin.tped
sed '0~2d' admixed_ceu_yri_chr1.tped       |sed '0~2d' >admixed_ceu_yri_chr1_short.tped       
sed '0~2d' ceu_chr1.tped                   |sed '0~2d' >ceu_chr1_short.tped                   
sed '0~2d' yri_chr1.tped                   |sed '0~2d' >yri_chr1_short.tped                   
cp admixed_ceu_yri_chr1_origin.tfam  admixed_ceu_yri_chr1_short_origin.tfam
cp admixed_ceu_yri_chr1.tfam         admixed_ceu_yri_chr1_short.tfam
cp ceu_chr1.tfam                     ceu_chr1_short.tfam       
cp yri_chr1.tfam                     yri_chr1_short.tfam                   

###################
#RUN STRUCTURE
###################

#Combine into one file readable by STRUCTURE
plink --tfile admixed_ceu_yri_chr1_short --noweb --recode12 --out plink1
plink --tfile ceu_chr1_short --noweb --recode12 --out plink2
plink --tfile yri_chr1_short --noweb --recode12 --out plink3
cut -f 2-3 plink1.map |awk 'p{print $1" "$2-p} {p=$2}'|transpose -s" ">plink.ped
x=`head -1 plink1.map |cut -f2`
sed -i -e "s/^rs/$x rs/" -e "s/^0/-1 0/" plink.ped
cat plink1.ped plink2.ped plink3.ped >>plink.ped
#hand edit first two lines to add a the first rsID from plink1.map and -1 for firt two lines



#prepare structure
echo "#define POPALPHAS 1 
#define LINKAGE 1  //run the linkage model instead of admix
#define PRINTLIKES 0
#define PRINTQHAT 1   //Puts alpha estimates in separate file _q
#define COMPUTEPROB 0 //does not print likelehood (speeds up 10-15%)">extraparams

echo "#define MAXPOPS 2
#define BURNIN 20000
#define NUMREPS 10000
#define ADMBURNIN 5000
#define INFILE plink.ped
#define OUTFILE  structure_french_yoruba

#define NUMLOCI 10879   //L, Might need to be changed
#define NUMINDS 43      //N
#define PLOIDY 2
#define MISSING -9
#define ONEROWPERIND 1 //1 means two columns per indiv. 0=2 rows per indiv.
#define LABEL 1
#define POPDATA 0
#define POPFLAG 0
#define LOCDATA 0
#define PHENOTYPE 0
#define EXTRACOLS 5
#define SITEBYSITE 1
#define ECHODATA 1

#define LOG10RSTART -2.0
#define LOG10RMAX 2.0
#define LOG10RMIN -4.0
#define LOG10RPROPSD 0.1


#define MARKERNAMES 1
#define RECESSIVEALLELES 0
#define MAPDISTANCES 1
#define PHASED 0
#define PHASEINFO 0
#define MARKOVPHASE 0
">mainparams

time ~/tmp/Structure/console/structure -L `cat plink1.map|wc -l` 

#evaluate results
echo "import numpy as np
Nadm=2  #wc -l plink1.ped
Nceu=22 #wc -l plink2.ped
Nyri=19 #wc -l plink3.ped
correctFile='admixed_ceu_yri_chr1_short_origin.tped'

results=[l.split()[2:]  for l in open('structure_french_yoruba_ss')]
results=np.asarray([l for l in results if len(l)>0], np.float)
results=np.argmax(results, 1)
results[results==2]=1  #Change all heterozygotes to same classification
results.shape=(43,-1)
ceu=0 if np.sum(results[Nadm:Nceu+Nadm,:]==0) > np.sum(results[Nadm:Nceu+Nadm,:]==3) else 3
yri=0 if np.sum(results[Nadm+Nceu:, :]==0) > np.sum(results[Nadm+Nceu:,:]==3) else 3



correct=np.asarray([l.split()[4:]  for l in open(correctFile)], np.float)
correct=correct[:, 0::2]+correct[:, 1::2]  #Combine to find 
idxCeu=correct==0
idxYri=correct==2
correct[idxCeu]=ceu
correct[idxYri]=yri

print 'Structure: %0.3g%%' % (100*np.sum(results[0:2]==correct.T)/float(correct.size))">evaluateResults.py
python evaluateResults.py
