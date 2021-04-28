#!/bin/bash
#PBS -l nodes=1:ppn=20 ## Requests 1 processor on 20 node
#PBS -l walltime=24:00:00 ## Requests 24 hours of walltime
#PBS -l pmem=8gb ## Requests 8 gigabytes of memory per process
#PBS A ## Specifies the allocation. Use A open for open queue
#PBS -j oe ## Requests that regular output and terminal output go to the same file
## The following is the body of the script. By default PBS scripts execute in your home directory, not the
## directory from which they were submitted. The following line places you in the directory from which the job
## was submitted.
cd $PBS_O_WORKDIR
##For processing and ancestry estimation of msprime generated genotype data
##To use this script, one need to put admixture proportion (0.1, 0.2, ..., 0.5)

prop=$1
filename=simADM_${prop}prop
echo "The name is ${filename}"

for chr in {1..10}
do
	echo "Processing chr${chr} of ${filename}"
	name=${filename}_chr${chr}
	echo "Change REF/ALT to A/T in the exported VCF file"
# save header to a separate file
	head -n6 ${name}.vcf > header_${name}.txt
#use awk to replace the ref/alt alleles with A/T, then concatenate with header 
#also, add ID column to each variant
	cat header_${name}.txt <(cat ${name}.vcf \
		| awk -v OFS="\t" 'NR>6 {$3=$1"_"$2"_A_T"; $4="A"; $5="T"; print ;}') \
	| bgzip > ${name}_AT.vcf.gz

	#Index the vcf file
	bcftools index ${name}_AT.vcf.gz
done

# Merging vcf files and convert to PLINK format and ancestry estimate


echo "merging vcf files"
bcftools concat ${filename}_chr*_AT.vcf.gz \
-o ${filename}_AT.vcf.gz \
-O z

echo "converting to PLINK format"
plink \
--vcf ${filename}_AT.vcf.gz \
--make-bed \
--out ${filename}_AT


echo "MAF filtering"
plink \
--bfile ${filename}_AT \
--maf 0.05 \
--make-bed \
--out ${filename}_AT_maf0.05


echo "LD pruning"
plink \
--bfile ${filename}_AT_maf0.05 \
--indep-pairwise 100 10 0.1 \
--out ${filename}_AT_maf0.05

plink \
--bfile ${filename}_AT_maf0.05 \
--extract ${filename}_AT_maf0.05.prune.in \
--make-bed \
--out ${filename}_AT_maf0.05_pruned

# Run admixture and plot
echo "run ancestry estimation"
admixture --cv -j6 ${filename}_AT_maf0.05_pruned.bed 2


# Plot the result
echo "plotting for the ancestry estimation, AFR anc distribution for ADM, and density plot of ASW_AA"
Rscript plot_admixture_dist.R ${filename}_AT_maf0.05_pruned.2.Q ${filename}_AT_maf0.05_pruned.fam

echo "Done ancestry estimation!"












