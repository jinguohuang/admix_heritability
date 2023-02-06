# estimate vg due to lanc and organize result
#!/bin/bash

filename1=admix_seed${1}_P${2}_t${3}

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

paste -d"\t" sim_${filename1}_exp_obs.txt <(echo -e "Vbetween" | cat - <(sed -n '2p' ${filename1}.reml.hsq | awk '{print $2}')) > summary_${filename1}.txt
