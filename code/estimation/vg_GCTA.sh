# estimate vgamma with gcta with and without ancestry as covariate
#!/bin/bash

# make sure have plink2 and gcta installed and executable
# one run for one replicate for each generation

model=${1} #HI or CGF
theta=${2} #0.1 0.2 0.5
gen=${3} #10, 20, 50, 100
P=${4} #0, 0.3, 0.6, 0.9
cov=${5} #pos or neg
seed=${6}
t=${7}

#source files
plinkdir=../data/theta0.5_gen20/plink_lanc/${model} 
phenodir=../data/theta0.5_gen20/pheno/${model}
gancdir=../data/theta0.5_gen20/ganc/${model}
#output files
grmdir=../data/theta0.5_gen20/gcta_lanc/grm
remldir=../data/theta0.5_gen20/gcta_lanc/reml_results/${model}
remlwdir=../data/theta0.5_gen20/gcta_lanc/reml_results_wganc/${model}


filename=admix_${model}_theta${theta}_gen${gen}_P${P}_${cov}_seed${seed}_t${t}_lanc
filename1=admix_${model}_theta${theta}_gen${gen}_P${P}_${cov}_seed${seed}_t${t}

./plink2 --import-dosage ${plinkdir}/${filename}.dosage format=1 single-chr=1 noheader \
--fam ${plinkdir}/${filename}.tfam \
--make-bed --out ${plinkdir}/${filename}

# construct grm with gcta
./gcta64 --bfile ${plinkdir}/${filename} \
--make-grm --out ${grmdir}/${filename} \
--thread-num 10

# estimate vgamma without ancestry

./gcta64 --reml \
--grm ${grmdir}/${filename} \
--pheno ${phenodir}/${filename1}.pheno \
--out ${remldir}/${filename}.reml \
--thread-num 10 \
--reml-no-constrain

# estimate Vgamma with ganc as covariate

./gcta64 --reml \
--grm ${grmdir}/${filename} \
--pheno ${phenodir}/${filename1}.pheno \
--qcovar ${gancdir}/${filename1}.ganc \
--out ${remlwdir}/${filename}.ganc.reml \
--thread-num 10 \
--reml-no-constrain

