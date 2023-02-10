# estimate vg with gcta with and without ancestry as covariate
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

filename=../data/admix_${model}_theta${theta}_gen${gen}_P${P}_${cov}_seed${seed}_t${t}

./plink2 --import-dosage ${filename}.dosage format=1 single-chr=1 noheader \
--fam ${filename}.tfam \
--make-bed --out ${filename}

# construct grm with gcta
./gcta64 --bfile ${filename} \
--make-grm --out ${filename} \
--thread-num 10

# estimate vg without ancestry

./gcta64 --reml \
--grm ${filename} \
--pheno ${filename}.pheno \
--out ${filename}.reml \
--thread-num 10 \
--reml-no-constrain

#paste -d"\t" sim_${filename}_exp_obs.txt <(echo -e "GCTA" | cat - <(sed -n '2p' ${filename}.reml.hsq | awk '{print $2}')) > summary_${filename}.txt


# estimate Vg with ganc as covariate

./gcta64 --reml \
--grm ${filename} \
--pheno ${filename}.pheno \
--qcovar ${filename}.qcovar \
--out ${filename}.ganc.reml \
--thread-num 10 \
--reml-no-constrain

#paste -d"\t" summary_${filename}.txt <(echo -e "GCTA_ganc" | cat - <(sed -n '2p' ${filename}.ganc.reml.hsq | awk '{print $2}')) > summary_${filename}_ganc.txt
