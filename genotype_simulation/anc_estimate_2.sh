#!/bin/bash
##For processing and ancestry estimation of msprime generated genotype data

name=AA_100Mb_adm_hap

echo "Change REF/ALT to A/T in the exported VCF file"
# save header to a separate file
head -n6 ${name}.vcf > header_${name}.txt
#use awk to replace the ref/alt alleles with A/T, then concatenate with header 
#also, add ID column to each variant
cat header_${name}.txt <(cat ${name}.vcf \
	| awk -v OFS="\t" 'NR>6 {$3=$1"_"$2"_A_T"; $4="A"; $5="T"; print ;}') > ${name}_AT.vcf

echo "converting ${name} to PLINK format"
plink \
--vcf ${name}_AT.vcf \
--make-bed \
--out ${name}_AT

echo "MAF filtering"
plink \
--bfile ${name}_AT \
--maf 0.05 \
--make-bed \
--out ${name}_AT_maf0.05


echo "LD pruning"
plink \
--bfile ${name}_AT_maf0.05 \
--indep-pairwise 100 10 0.1 \
--out ${name}_AT_maf0.05

plink \
--bfile ${name}_AT_maf0.05 \
--extract ${name}_AT_maf0.05.prune.in \
--make-bed \
--out ${name}_AT_maf0.05_pruned

# Run admixture and plot
echo "run ancestry estimation"
admixture --cv -j6 ${name}_AT_maf0.05_pruned.bed 2


# Plot the result
echo "plotting for the ancestry estimation"
Rscript plot_admixture.R ${name}_AT_maf0.05_pruned.2.Q ${name}_AT_maf0.05_pruned.fam

echo "Done ancestry estimation!"
