#!/bin/bash
# extract trait specific allele and calculate gvalue in ASW

trait=$1
plink2="../../code/plink2"

plink --bfile ../../Ref_1KG_hg38/ASW --extract ${trait}_allele_beta.txt --make-bed --out ASW_${trait} --allow-extra-chr
# REMOVE DUPLICATED ID
${plink2} --bfile ASW_${trait} --rm-dup exclude-mismatch --make-bed --out ASW_${trait}_rmdup
# calculate gvalue in ASW,  The 'sum' modifier causes sums to be reported instead.
plink --bfile ASW_${trait}_rmdup --score ${trait}_allele_beta.txt header sum --out ASW_${trait}
#calculate the variance
Rscript ../../code/vgvalue_CI_calculator.R ${trait}