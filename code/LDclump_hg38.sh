#!/bin/bash
# trait data processing for traits in hg38
## prep for LD clumping: get SNP, P ready
## LD clumping: pvalue and physical distance
## process to get allele freq and calculate vgp vgvalue


trait=$1 #HDL LDL TG TC
# process trait file for GWAS catlog data
echo "processing trait ${trait} file"
awk -F"\t" '{OFS="\t"; print $12":"$13, $21, $27, $28, $31, $32}' ${trait}.tsv |\
 awk 'BEGIN {OFS="\t"} NR>1 {split($2, a, "-"); $2=a[2]} 1' |\
  grep -v "?" |\
   awk 'NR>1{OFS="\t"; print $1,$2,$3,$4,$5,$6=($8=="increase")?1:(($8=="decrease")?-1:$8)}' |\
      awk '{OFS="\t";print $1,$2,$3,$4,$5=$5*$6}' |\
         cat <(echo -e "SNP\tA1\tA1_FREQ\tP\tBETA") - > ${trait}.txt

wc -l ${trait}.txt

# ld clump 1
plink --bfile ../../Ref_1KG_hg38/GBR_nodup \
--clump ${trait}.txt \
--clump-p1 5e-8 \
--clump-r2 0.05 \
--clump-kb 250 \
--allow-extra-chr \
--out ${trait}


# ld clump 2
plink --bfile ../../Ref_1KG_hg38/GBR_nodup \
--clump ${trait}.clumped \
--clump-kb 100 \
--allow-extra-chr \
--out ${trait}_100kb

# output the clumped snp only
echo "extracting clumped SNPs from ${trait} file..."
awk 'NR==FNR{A[$1]=$0;next}$1=A[$3]' ${trait}.txt <(awk 'NF' ${trait}_100kb.clumped) |\
awk 'OFS="\t"{print $1,$2, $5}' > ${trait}_allele_beta.txt


# calculate f1 f2
../../code/Allelef1f2_hg38.sh ${trait} CEU
../../code/Allelef1f2_hg38.sh ${trait} YRI

# calculate vgp
Rscript ../../code/vgp_calculator.R ${trait}
cat ../${trait}_vgp.txt

# calculate vgvalue
../../code/gvalue_hg38.sh ${trait}
cat ../${trait}_vgvalueCI.txt
