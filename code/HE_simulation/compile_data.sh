#!/bin/bash

model=${1} #HI or CGF
P=${2} #0, 0.3, 0.6, 0.9
cov=${3} #pos or neg or zero


#dir address
HEdir=/home/aazaidi/klema030/AdjustedHE/output/${model}
tmpdir=/home/aazaidi/klema030/AdjustedHE/output/tmp
outdir=/home/aazaidi/klema030/AdjustedHE/output


# compile for HE reg scripts
touch admix_${model}_theta0.5_gen20_P${P}_${cov}.txt
head -n1 -q first_row.txt > admix_${model}_theta0.5_gen20_P${P}_${cov}.txt
for i in {1..10}; do for j in {0..20}; do tail -n 1 -q ./${model}/admix_${model}_theta0.5_gen20_P${P}_${cov}_seed${i}_t${j}.txt >> admix_${model}_theta0.5_gen20_P${P}_${cov}.txt; done; done


#compile for GCTA HE reg
filename=admix_${model}_theta0.5_gen20_P${P}_${cov}

for i in {1..10}; do echo -e "t\tvg_HEreg" > ${tmpdir}/${filename}_seed${i}.HEreg.txt; done
for i in {1..10}; do for ((t = 0; t <= 20; t++)); do printf ${t}"\t" | cat - <(sed -n '4p' ${HEdir}/${filename}_seed${i}_t${t}.HEreg.HEreg | awk '{print $2}')  >> ${tmpdir}/${filename}_seed${i}.HEreg.txt; done; done

head -n 1 ${tmpdir}/admix_CGF_theta0.5_gen20_P0_zero_seed1.HEreg.txt >> ${tmpdir}/admix_${model}_theta0.5_gen20_P${P}_${cov}.HEreg.txt 
tail -n +2 -q ${tmpdir}/admix_${model}_theta0.5_gen20_P${P}_${cov}_seed*.HEreg.txt >> ${tmpdir}/admix_${model}_theta0.5_gen20_P${P}_${cov}.HEreg.txt

(head -n 1 ${tmpdir}/admix_${model}_theta0.5_gen20_P${P}_${cov}.HEreg.txt && tail -n +2 -q ${tmpdir}/admix_${model}_theta0.5_gen20_P${P}_${cov}.HEreg.txt | sort -n -k 1) >> ${outdir}/admix_${model}_theta0.5_gen20_P${P}_${cov}.HEreg.txt


