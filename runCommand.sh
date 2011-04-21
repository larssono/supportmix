#./build/exe.linux-x86_64-2.6/SupportMix -w 400 ceu_yri/ancestral_ceu.chr22.csv ceu_yri/ancestral_yri.chr22.csv /home/local/CB/fja32/workspace/SupportMix/src/ceu_yri/admixed_ceu_yri.chr22.csv

#SupportMix -w100 -a data/tests/admixed_origin_ceu_yri.chr22.csv.gz \
#    data/tests/ancestral_ceu.chr22.csv.gz \
#    data/tests/ancestral_yri.chr22.csv.gz \
#    data/tests/admixed_ceu_yri.chr22.csv.gz 

#using TPED file for admixed
SupportMix -w100 -a data/tests/admixed_origin_ceu_yri.chr22.tped.gz \
    data/tests/ancestral_ceu.chr22.csv.gz \
    data/tests/ancestral_yri.chr22.csv.gz \
    data/tests/admixed_ceu_yri.chr22.csv.gz

