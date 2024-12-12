#!/bin/bash

model=${1} #HI or CGF
P=${2} #0, 0.3, 0.6, 0.9
cov=${3} #pos or neg or zero
grm=${4} #varX or ld



#dir address
klemdir=/home/aazaidi/klema030/AdjustedHE
HEdir=${klemdir}/output/${model}/gcta/HE_GRM${grm}
tmpdir=${klemdir}/output/tmp/HE_GRM${grm}
outdir=${klemdir}/output/summary/${model}

filename=admix_${model}_theta0.5_gen20_P${P}_${cov}

# compile for HE reg adjusted script
touch ${tmpdir}/admix_${model}_theta0.5_gen20_P${P}_${cov}.HE_GRM${grm}_ganc.txt
head -n1 -q ${klemdir}/output/first_row_${grm}.txt > ${tmpdir}/admix_${model}_theta0.5_gen20_P${P}_${cov}.HE_GRM${grm}_ganc.txt
for i in {1..10}; do for j in {0..20}; do tail -n 1 -q ${HEdir}/admix_${model}_theta0.5_gen20_P${P}_${cov}_seed${i}_t${j}.txt >> ${tmpdir}/admix_${model}_theta0.5_gen20_P${P}_${cov}.HE_GRM${grm}_ganc.txt; done; done

#compile for HE reg unadjusted script (noa)
touch ${tmpdir}/admix_${model}_theta0.5_gen20_P${P}_${cov}.HE_GRM${grm}_noa.txt
head -n1 -q ${klemdir}/output/first_row_${grm}_noa.txt > ${tmpdir}/admix_${model}_theta0.5_gen20_P${P}_${cov}.HE_GRM${grm}_noa.txt
for i in {1..10}; do for j in {0..20}; do tail -n 1 -q ${HEdir}/admix_${model}_theta0.5_gen20_P${P}_${cov}_seed${i}_t${j}noa.txt >> ${tmpdir}/admix_${model}_theta0.5_gen20_P${P}_${cov}.HE_GRM${grm}_noa.txt; done; done

# #compile for gcta HEreg grms
echo -e "t\tvg_GRM${grm}" > ${tmpdir}/${filename}.GRM${grm}.txt

for ((s = 1; s <= 10; s++)); do
    for ((t = 0; t <= 20; t++)); do 
        value1=$(sed -n '4p' "${HEdir}/${filename}_seed${s}_t${t}.HEreg.HEreg" | awk '{print $2}')
        echo -e "${t}\t${value1}" >> "${tmpdir}/${filename}.GRM${grm}.txt"
    done
done

#Join files together into one
join ${tmpdir}/admix_${model}_theta0.5_gen20_P${P}_${cov}.HE_GRM${grm}_ganc.txt ${tmpdir}/${filename}.GRM${grm}.txt > ${outdir}/${filename}.summary_GRM${grm}.txt

#join unadjusted to summary
cat ${tmpdir}/admix_${model}_theta0.5_gen20_P${P}_${cov}.HE_GRM${grm}_noa.txt | cut -f 1,6 > ${tmpdir}/admix_${model}_theta0.5_gen20_P${P}_${cov}.HE_GRM${grm}.txt

join ${outdir}/${filename}.summary_GRM${grm}.txt ${tmpdir}/admix_${model}_theta0.5_gen20_P${P}_${cov}.HE_GRM${grm}.txt > ${outdir}/${filename}.summary_GRM${grm}_all.txt

echo "done!"
