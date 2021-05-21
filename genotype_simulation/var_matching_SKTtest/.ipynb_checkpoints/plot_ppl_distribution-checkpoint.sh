#!/bin/bash
# for plotting local to global distribution of 1Kppl 
# from simulated local ancestry data
# starting file: local ancestry file of each chromosome (x100)
# step 1: Extract admixed ppl from all, previous 400 are AFR, EUR haplotypes
# step 2: calculate AFR anc percentage for each person (by column)
# step 3: concatenate all chromosomes and calculate average by column
# step 4: plot for distribution 

model=$1
echo "Processing for ${model} model ..."
for chr in {1..100}
#for chr in 1
do
	echo "extract admixed ppl for chr${chr}"
	locfile=loc_anc_simAA_1Mb_1Kppl_${model}_chr${chr}.reformat
	mapfile=loc_anc_simAA_1Mb_1Kppl_${model}_chr${chr}.map
	cut -f401-2400 ${locfile} > ${locfile}_adm-tmp
	echo "calculate AFR anc percentage for each person"
	#Iterate over number of fields and if the field is equal to 0 
	#increment the array. Then on the end print the array devided
	#by the number of rows.
	line=$(awk 'END{print NR}' ${locfile}_adm-tmp)
	awk -v line=${line} '{ for (i = 1; i <= NF; ++i) { if($i == 0) { ++c[i]; } } } 
	END{ for (i = 1; i <= NF; ++i) { printf "%.3f%s", c[i]/line, i!=NF ? OFS : ORS; } }' ${locfile}_adm-tmp > ${locfile}_adm_anc0-tmp
	echo "Done for chr${chr}"
done

echo "finish calculating ancestry by person for each chromosome"
echo "put every chromosome together"
cat loc_anc_simAA_1Mb_1Kppl_${model}_chr*.reformat_adm_anc0-tmp >> simAA_1Mb_1Kppl_${model}_ancBYppl.txt
echo "remove all temporary files"
rm loc_anc_simAA_1Mb_1Kppl_${model}_*-tmp
echo "plot for distribution and calculate mean"
Rscript plot_anc_dist.R simAA_1Mb_1Kppl_${model}_ancBYppl.txt
echo "Done plotting!" 



