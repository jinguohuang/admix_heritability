#!/bin/bash

#organize results
model=${1} #HI or CGF
grm=${2} #varX or ld


#dir address
klemdir=/home/aazaidi/klema030/AdjustedHE
summarydir=${klemdir}/output/summary/${model}
outdir=${klemdir}/output

echo "organizing results"

head -1 ${summarydir}/admix_${model}_theta0.5_gen20_P0_pos.summary_GRM${grm}_all.txt > ${outdir}/admix_${model}_theta0.5_gen20.summary_GRM${grm}.txt

for P in 0 0.3 0.6 0.9
do
    for cov in pos neg zero
    do
        cat <(tail -n +2 ${summarydir}/admix_${model}_theta0.5_gen20_P${P}_${cov}.summary_GRM${grm}_all.txt) >> ${outdir}/admix_${model}_theta0.5_gen20.summary_GRM${grm}.txt
    done
done

echo "done!"
