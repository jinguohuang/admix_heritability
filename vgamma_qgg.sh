#!/bin/bash
# qgg estimate vgamma
seed=${1}
P=0.9
p=0


echo "simulate admix lanc for seed ${seed} and P ${P}"
Rscript Vgeno_CGF_output.R --seed ${seed} -P ${P} -F 0.2 -l 100 -n 1e3 -t 20 -p ${p} --theta 0.5

echo "estimate vg due to lanc"
# estimate vg due to lanc
for t in {0..20}; do Rscript vgamma_qgg.R ${seed} ${t} $P; done

# add vgammahat results
for t in {0..20}; do 
	filename1=admix_seed${seed}_P${P}_psel0_m10_t${t}
	paste -d"\t" sim_${filename1}_exp_obs.txt <(echo -e "vgammahat" |\
		cat - ${filename1}_vgammahat.txt ) > summary_${filename1}.txt
done

echo "orgnize result"
filename=summary_admix_seed${seed}_P${P}_psel${p}_m10
head -1 ${filename}_t0.txt > ${filename}
for t in {0..20}; do cat <(tail -1 ${filename}_t${t}.txt) >> ${filename}; done


echo "remove tmp files"
rm ${filename}_t*.txt
rm sim_*.txt
rm admix_*