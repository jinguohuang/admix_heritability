#!/bin/bash

theta=${1}
venv=${2}
vgb=${3}
pheno=${4}

filename1=admix_prop${theta}_vgb${vgb}

Rscript	initial_simulations/simpleSim_addenv.R --vgb ${vgb} --theta ${theta} --venv ${venv} --seed 1

plink2 --import-dosage ${filename1}.dosage format=1 single-chr=1 noheader \
--fam ${filename1}.tfam \
--make-bed --out ${filename1}

gcta64 --bfile ${filename1} \
--make-grm --out ${filename1}

gcta64 --reml \
--grm ${filename1} \
--pheno ${filename1}.pheno \
--out ${filename1}.reml \
--mpheno $pheno \
--thread-num 8

