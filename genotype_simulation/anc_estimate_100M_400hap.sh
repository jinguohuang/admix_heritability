#!/bin/bash
##For processing and ancestry estimation of msprime generated genotype data
##To use this script, one need to put random or 10deme to run

for chr in 1
do
	echo "Processing chr${chr}"
	name=AA_100Mb_$1_400hap_chr${chr}
	echo "Change REF/ALT to A/T in the exported VCF file"
# save header to a separate file
	head -n6 ${name}.vcf > header_${chr}.txt
#use awk to replace the ref/alt alleles with A/T, then concatenate with header 
#also, add ID column to each variant
	cat header_${chr}.txt <(cat ${name}.vcf \
		| awk -v OFS="\t" 'NR>6 {$3=$1"_"$2"_A_T"; $4="A"; $5="T"; print ;}') \
	| bgzip > ${name}_AT.vcf.gz

	#Index the vcf file
	bcftools index ${name}_AT.vcf.gz
done

# Merging vcf files and convert to PLINK format and ancestry estimate
name=AA_100Mb_$1_400hap_chr1
echo "The name is ${name}"

#echo "merging vcf files"
#bcftools concat AA_100Mb_random_10000ppl_chr*_AT.vcf.gz \
#-o ${name}_AT.vcf.gz \
#-O z

echo "converting to PLINK format"
plink \
--vcf ${name}_AT.vcf.gz \
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
echo "plotting for the ancestry estimation, AFR anc distribution for ADM, and density plot of ASW_AA"
Rscript plot_admixture_dist_density_kstest_2.R ${name}_AT_maf0.05_pruned.2.Q ${name}_AT_maf0.05_pruned.fam ASW_dip_anc1.txt

echo "Done ancestry estimation!"











