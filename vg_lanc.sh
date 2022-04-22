# estimate vg due to lanc
#!/bin/bash

filename1=$1
./plink2 --import-dosage ${filename1}.dosage format=1 single-chr=1 noheader \
--fam ${filename1}.tfam \
--make-bed --out ${filename1}

gcta64 --bfile ${filename1} \
--make-grm --out ${filename1}

gcta64 --reml \
--grm ${filename1} \
--pheno ${filename1}.pheno \
--out ${filename1}.reml \
--thread-num 8

