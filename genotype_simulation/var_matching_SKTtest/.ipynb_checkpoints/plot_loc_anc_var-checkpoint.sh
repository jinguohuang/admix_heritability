##!/bin/bash
# script to processing local ancestry and calculate variance and mean  
# step 1: Extract admixed ppl from all, previous 400 are AFR, EUR haplotypes
# step 2: calculate AFR anc percentage for each locus
# step 3: add SNP map info to the list
# step 4: plot and calculate variance, mean
model=$1
echo "Processing for ${model} model ..."
for chr in {1..100}
#for chr in 1
do
	echo "extract admixed ppl for chr${chr}"
	locfile=loc_anc_simAA_1Mb_1Kppl_${model}_chr${chr}.reformat
	mapfile=loc_anc_simAA_1Mb_1Kppl_${model}_chr${chr}.map
	cut -d' ' -f401-2400 ${locfile} > ${locfile}_adm-tmp
	echo "calculate AFR anc percentage for each locus"
	awk -F '0' '{print (NF-1)/2400}' ${locfile}_adm-tmp > ${locfile}_adm_anc0-tmp
	echo "add SNP map info to the list"
	paste -d ' ' ${mapfile} ${locfile}_adm_anc0-tmp | awk -v awkvar=${chr} '{print awkvar"_"$1,awkvar,$1,$2}' > ${locfile}_adm_anc0_annotated-tmp
	echo "Done for chr${chr}"
done

echo "finish calculating ancestry at each locus"
echo "put every chromosome together"
cat loc_anc_simAA_1Mb_1Kppl_${model}_chr*.reformat_adm_anc0_annotated-tmp >> simAA_1Mb_1Kppl_${model}_ancBYlocus.txt
echo "remove all temporary files"
rm loc_anc_simAA_1Mb_1Kppl_${model}_*-tmp
echo "plot and calculate variance, mean"
Rscript plot_loc_anc.R simAA_1Mb_1Kppl_${model}_ancBYlocus.txt
echo "Done plotting!" 
