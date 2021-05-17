##!/bin/bash
##For calculating multiple chromosomes of ASW anc proportion
#for chr in {1..22}
for chr in 2
do
	#echo "gunzip for chr${chr}"
	#gunzip AMR_chr${chr}.rfmix.1.Viterbi.txt.gz 
	echo "Extract ASW from AMR"
	cut -d' ' -f2293-2424 AMR_chr${chr}.rfmix.1.Viterbi.txt > ASW_chr${chr}.rfmix.1.Viterbi
	echo "Calculate AFR anc percentage for each haplotype"
	#Iterate over number of fields and if the field is equal to 1 
	#increment the array. Then on the end print the array devided
	#by the number of rows.
	# ASW: 1-AFR
	line=$(awk 'END{print NR}' ASW_chr${chr}.rfmix.1.Viterbi)
	awk -v line=${line} '{ for (i = 1; i <= NF; ++i) { if($i == 1) { ++c[i]; } } } 
	END{ for (i = 1; i <= NF; ++i) { printf "%.3f%s", c[i]/line, i!=NF ? OFS : ORS; } }' ASW_chr${chr}.rfmix.1.Viterbi > ASW_chr${chr}.anc1
	echo "Done for chr${chr}"
done

echo "finish calculating ancestry by person for each chromosome"
echo "put every chromosome together"
cat ASW_chr*.anc1 >> ASW_ancBYppl.txt
echo "run SKT test for ASW and plot for distribution"
Rscript plot_skt_dist_ASW.R ASW_ancBYppl.txt
echo "Done plotting!"

