#!/bin/bash

#######################################################################
#Simulate all matings between 2 hgdp populations with >20 haplotypes, NGENS=5, ALPHA=0.5
#######################################################################
# mkdir data data/hgdp_ancestral data/hgdp2 data/hgdp3 data/hgdp_alpha data/hgdp_generations
# ln -s ../../../../human_genome_data/hapmap2/ data/hapmap2
# time python simulateHGDP.py  
# find data -name "*.csv"|xargs -n1 -P2 gzip

#######################################################################
#Run SupportMix and Lamp on simulated data
#######################################################################
python examineHGDP.py

#######################################################################
#Run structure on Qatari's
#######################################################################
# Fix to change directory
# plink --bfile ../CrystalQatar/qatar_unrelated --indep-pairwise 100 5 0.5
# plink --bfile ../CrystalQatar/qatar_unrelated --maf 0.05 --geno 0.05 --hwe 0.001 --extract plink.prune.in  --recode --allele1234  --thin 0.2 --missing-genotype -9
# ~/tmp/structure/console/structure -N 156 -L `cat plink.map|wc -l` -K 3 -o qatar_admix
# mv qatar_admix_q ../
# rm plink*

#######################################################################
#Admix map the qatari and run PCA 
#######################################################################
python examineQatar.py


#######################################################################
#Generate plots for paper
#######################################################################
python plotResults.py

