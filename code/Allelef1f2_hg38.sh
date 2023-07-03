#!/bin/bash
# extract trait specific allele and calculate freq
# 

trait=$1 # RBC, WBC, etc
pop=$2  # CEU or YRI
plink2="../../code/plink2"

echo "remove D/I in the ${trait} reference allele"
grep -v 'D\|I' ${trait}_allele_beta.txt > ${trait}_allele_beta_noDI.txt

${plink2} --bfile ../../Ref_1KG_hg38/${pop} --extract ${trait}_allele_beta_noDI.txt --make-bed --out ${pop}_${trait} --allow-extra-chr

${plink2} --bfile ${pop}_${trait} --rm-dup exclude-mismatch --make-bed --out ${pop}_${trait}_rmdup


plink --bfile ${pop}_${trait}_rmdup --reference-allele ${trait}_allele_beta_noDI.txt --freq --out ${pop}_${trait}


echo "get exclude SNP list if they don't have matched allele"
grep "Warning" ${pop}_${trait}.log | awk '{print $8}' | sed 's/\.$//'> ${pop}_${trait}.log_exclude 
plink --bfile ${pop}_${trait}_rmdup \
--exclude ${pop}_${trait}.log_exclude  \
--reference-allele ${trait}_allele_beta_noDI.txt \
--freq \
--out ${pop}_${trait}


