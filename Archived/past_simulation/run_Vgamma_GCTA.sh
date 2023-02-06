#!/bin/bash
# one run for one replicate
# HI
seed=${1}
P=${2}
p=${3}


echo "simulate admix lanc for seed ${seed} and P ${P}"
Rscript Vgeno_CGF_output.R --seed ${seed} -P ${P} -F 0.2 -l 100 -n 1e4 -t 20 -p ${p} --theta 0.5

echo "estimate vg due to lanc"
# estimate vg due to lanc
for t in {0..20}; do ./vgamma_GCTA.sh ${seed} ${P} $t; done

echo "orgnize result"
filename=summary_admix_seed${seed}_P${P}_psel${p}_m10
head -1 ${filename}_t0_ganc.txt > ${filename}
for t in {0..20}; do cat <(tail -1 ${filename}_t${t}_ganc.txt) >> ${filename}; done

echo "remove tmp files"
rm ${filename}_t*.txt
rm sim_*.txt
