#!/bin/bash
# estimate va with gcta with Haseman-Elston regression
# one run for one replicate for each generation

model=${1} #HI or CGF
#theta=${2} #0.1 0.2 0.5
#gen=${3} #10, 20, 50, 100
P=${2} #0, 0.3, 0.6, 0.9
cov=${3} #pos or neg
seed=${4}
t=${5}

#source files
huandir=/home/aazaidi/huan2788/admix_heritability/data/theta0.5_gen20
plinkHE=/home/aazaidi/klema030/AdjustedHE/output/plink/${model}
plinkdir=${huandir}/plink/${model} 
phenodir=${huandir}/pheno/${model}
gancdir=${huandir}/ganc/${model}
#output files
grmdir=${huandir}/gcta/grm
grmHE=/home/aazaidi/klema030/AdjustedHE/output/grm/${model}
grmHEscaled=/home/aazaidi/klema030/AdjustedHE/output/grm/${model}/scaled
HEdir=/home/aazaidi/klema030/AdjustedHE/output/${model}/gcta
HEdirscaled=/home/aazaidi/klema030/AdjustedHE/output/${model}/gcta/scaled


filename=admix_${model}_theta0.5_gen20_P${P}_${cov}_seed${seed}_t${t}
filename1=admix_${model}_theta0.5_gen20_P${P}_${cov}_seed${seed}_t${t}
filenamescaled=admix_${model}_theta0.5_gen20_P${P}_${cov}_seed${seed}_t${t}_scaled

# ./plink2 --import-dosage ${plinkdir}/${filename}.dosage format=1 single-chr=1 noheader \
# --fam ${plinkdir}/${filename}.tfam \
# --make-bed --out ${plinkHE}/${filename}

# construct grm with gcta
# ./gcta64 --bfile ${plinkHE}/${filename} \
# --make-grm --out ${grmHE}/${filename} \
# --thread-num 10

# estimate va with HE using manual grm
./gcta64 --HEreg \
--grm ${grmHE}/${filename} \
--pheno ${phenodir}/${filename1}.pheno \
--out ${HEdirscaled}/${filename}.HEreg.plink \
--thread-num 10 

# # estimate va with HE
# ./gcta64 --HEreg \
# --grm ${grmHE}/${filename} \
# --pheno ${phenodir}/${filename1}.pheno \
# --out ${HEdir}/${filename}.HEreg.plink \
# --thread-num 10 
