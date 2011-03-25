#!/bin/bash

#######################################################################
# Set up directory of data
#######################################################################
mkdir data 
cd data
ln -s ../../../../human_genome_data/hapmap2/ 
ln -s ../../../../human_genome_data/HumanPlinkData/CrystalQatar/ 
ln -s ../../../../human_genome_data/HumanPlinkData/HGDP/ HGDP_raw_data
wget http://rosenberglab.bioinformatics.med.umich.edu/data/rosenbergEtAl2005/rosenbergEtAl2005.coordinates.txt  #This file needs hand editing as the formating not consistent throughout
cd ..

#######################################################################
#Run structure on Qatari's
#######################################################################
mkdir data/structure/
cd data/structure/
echo "#define POPALPHAS 1 
#define PRINTLIKES 0
#define PRINTQHAT 1   //Puts alpha estimates in separate file _q
#define COMPUTEPROB 0 //does not print likelehood (speeds up 10-15%)">extraparams
echo "#define MAXPOPS 2
#define BURNIN 20000
#define NUMREPS 10000
#define INFILE plink.ped
#define OUTFILE qatar_admix

#define NUMLOCI 5457   //Might need to be changed
#define NUMINDS 156
#define PLOIDY 2
#define MISSING -9
#define ONEROWPERIND 1
#define LABEL 1
#define POPDATA 0
#define POPFLAG 0
#define LOCDATA 0
#define PHENOTYPE 0
#define EXTRACOLS 5

#define MARKERNAMES 0
#define RECESSIVEALLELES 0
#define MAPDISTANCES 0
#define PHASED 0
#define PHASEINFO 0
#define MARKOVPHASE 0
">mainparams
plink --bfile ../CrystalQatar/qatar_unrelated --indep-pairwise 100 5 0.5
plink --bfile ../CrystalQatar/qatar_unrelated --maf 0.05 --geno 0.05 --hwe 0.001 --extract plink.prune.in  --recode --allele1234  --thin 0.2 --missing-genotype -9
time ~/tmp/Structure/console/structure -N 156 -L `cat plink.map|wc -l` -K 3 -o qatar_admix 
mv qatar_admix_q ../
rm plink*
cd ../../

#######################################################################
#Admix map the qatari and run PCA 
#######################################################################
mkdir data/qatarSupportMix
time python examineQatar.py


#######################################################################
#Simulate all matings between 2 hgdp populations with >20 haplotypes, NGENS=5, ALPHA=0.5
#######################################################################
mkdir data/hgdp_ancestral data/hgdp2 data/hgdp3 data/hgdp_alpha data/hgdp_generations
time python simulateHGDP.py  
find data -name "*.csv"|xargs -n1 -P2 gzip

#######################################################################
#Run SupportMix and Lamp on simulated data
#######################################################################
time python examineHGDP.py

#######################################################################
#Generate plots for paper
#######################################################################
time python plotResults.py

