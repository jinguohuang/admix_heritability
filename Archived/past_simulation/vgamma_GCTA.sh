# estimate vg due to lanc and organize result
#!/bin/bash

filename1=admix_seed${1}_P${2}_psel0_m10_t${3}

./plink2 --import-dosage ${filename1}.dosage format=1 single-chr=1 noheader \
--fam ${filename1}.tfam \
--make-bed --out ${filename1}

gcta64 --bfile ${filename1} \
--make-grm --out ${filename1} \
--thread-num 10

gcta64 --reml \
--grm ${filename1} \
--pheno ${filename1}.pheno \
--out ${filename1}.reml \
--thread-num 10 

paste -d"\t" sim_${filename1}_exp_obs.txt <(echo -e "GCTA" | cat - <(sed -n '2p' ${filename1}.reml.hsq | awk '{print $2}')) > summary_${filename1}.txt


# calculate global ancestry from dosage file
# and prep the qcovar file for GCTA
# and estimate Vgamma with ganc as correction

paste -d'\t' <(cut -f1-2 ${filename1}.fam) \
<(cut -f4- ${filename1}.dosage |\
awk '{for (i=1;i<=NF;i++){
  sums[i]+=$i}
}END{
  for(i=1;i<=NF;i++){
    print(sums[i]/200)
  }}')\
>${filename1}.qcovar

gcta64 --reml \
--grm ${filename1} \
--pheno ${filename1}.pheno \
--qcovar ${filename1}.qcovar \
--out ${filename1}.ganc.reml \
--thread-num 10 

paste -d"\t" summary_${filename1}.txt <(echo -e "GCTA_ganc" | cat - <(sed -n '2p' ${filename1}.ganc.reml.hsq | awk '{print $2}')) > summary_${filename1}_ganc.txt
