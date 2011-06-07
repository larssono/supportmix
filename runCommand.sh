#./build/exe.linux-x86_64-2.6/SupportMix -w 400 ceu_yri/ancestral_ceu.chr22.csv ceu_yri/ancestral_yri.chr22.csv /home/local/CB/fja32/workspace/SupportMix/src/ceu_yri/admixed_ceu_yri.chr22.csv

SupportMix -w100 -c22 -m data/hapmap2 -a data/tests/admixed_origin_ceu_yri.chr22.tped.gz \
    data/tests/ancestral_ceu.tped.gz \
    data/tests/ancestral_yri.tped.gz \
    data/tests/admixed_ceu_yri.tped.gz

echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

SupportMix -w100 -t -a data/tests/admixed_origin_ceu_yri.chr22.tped.gz \
    -m data/hapmap2/ \
    data/tests/ancestral_ceu.chr22.csv.gz \
    data/tests/ancestral_yri.chr22.csv.gz \
    data/tests/admixed_ceu_yri.chr22.csv.gz

