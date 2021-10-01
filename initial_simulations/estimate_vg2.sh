#!/bin/bash

theta=${1}
venv=${2}
vgb=${3}
vgt=${4}
pheno=${5}

filename1=admix_prop${theta}_vgb${vgb}_vgt${vgt}_venv${venv}

Rscript simpleSim_addenv2.R --vgb ${vgb} --theta ${theta} --venv ${venv}

./plink2 --import-dosage ${filename1}.dosage format=1 single-chr=1 noheader \
--fam ${filename1}.tfam \
--make-bed --out ${filename1} 

gcta64 --bfile ${filename1} \
--make-grm --out ${filename1}

gcta64 --reml \
--grm ${filename1} \
--pheno ${filename1}.pheno \
--out ${filename1}.reml \
--mpheno ${pheno} \
--thread-num 8

Rscript simVSest_reformat.R ${theta} ${vgb} ${vgt} ${venv}