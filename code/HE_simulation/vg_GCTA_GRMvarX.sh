#!/bin/bash
# estimate va with gcta with var(x) scaled GRM
# estimate va with gcta with wganc

# make sure have gcta installed and executable
# one run for one replicate for each generation
# put this in code folder

model=${1} #HI or CGF
# theta=${2} #0.1 0.2 0.5
# gen=${3} #10, 20, 50, 100
P=${2} #0, 0.3, 0.6, 0.9
cov=${3} #pos or neg
seed=${4}
t=${5}

#source files
huandir=/home/aazaidi/huan2788/admix_heritability/data/theta0.5_gen20
klemdir=/home/aazaidi/klema030/AdjustedHE
plinkdir=${huandir}/plink/${model} 
phenodir=${huandir}/pheno/${model}
gancdir=${huandir}/ganc/${model}
#output files
grmdir=${klemdir}/output/grm/gcta/GRMvarX/${model}
HEdir=${klemdir}/output/${model}/gcta/HE_GRMvarX


filename=admix_${model}_theta0.5_gen20_P${P}_${cov}_seed${seed}_t${t}
filename1=admix_${model}_theta0.5_gen20_P${P}_${cov}_seed${seed}_t${t}

# construct grm with my script, it will read the dosage file, and output grm
module load R
Rscript GRMvarX.R ${filename} ${plinkdir} ${grmdir}

# estimate va with gcta
# feed special GRM.GZ
./gcta64 --HEreg \
--grm ${grmdir}/${filename} \
--pheno ${phenodir}/${filename1}.pheno \
--out ${HEdir}/${filename}.HEreg \
--thread-num 10 \


