#!/bin/bash
# one run for one replicate

model=${1} #HI or CGF
theta=${2} #0.1 0.2 0.5
gen=${3} #10, 20, 50, 100
P=${4} #0, 0.3, 0.6, 0.9
cov=${5} #pos or neg
seed=${6}


echo "simulate admix for seed ${seed} and P ${P}"
Rscript sim_admix.R -F 0.2 -l 1e3 -n 1e4 -p 0 -M ${model} --theta ${theta} -t ${gen} -P ${P} -C ${cov} --seed ${seed} 

