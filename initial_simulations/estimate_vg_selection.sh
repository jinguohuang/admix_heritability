#!/bin/bash

theta=${1}
fstc=${2}
fstn=${3}
pcausal=${4}
pheno=${5}

filename1=admix_prop${theta}_fstc${fstc}_fstn${fstn}_pcausal${pcausal}

Rscript simpleSim_selection.R --theta ${theta} --fstc ${fstc} --fstn ${fstn} --pcausal ${pcausal}

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

Rscript simVSest_reformat.R ${theta} ${fstc} ${fstn} ${pcausal}