##!/bin/bash
# script to reformat local ancestry 
# step 1: reformat the local ancestry to bed file for bedtools
# step 2: bedtools intersect to find overlapping variant and their local ancestry
# step 3: remove duplicates unmelt the output to 2 files: variants x individual format and the variants list, and calculate genetic distance
# step 4: remove all temporary files

# input1 is the local ancestry file from simulation
# input2 is the vcf file from simulation

bedfile=$1
# ATTENTION! THIS BEDFILE HAS NO HEADER.
# extract the chromosome number 
## loc_anc_simAA_1Mb_1Kppl_10deme_chr1.vcf.txt
## after "chr" before ".vcf"
chr=$(echo "$bedfile" | grep -oP '(?<=chr)\d+(?=\.vcf)')
## (?<=) is lookbehind, (?<=chr) matches "chr" before pattern.
## \d+ matches one or more number.
## (?=) is lookahead, (?=\.vcf) matches ".vcf" after pattern.

# step 1
echo "prepare ${bedfile} to bed file for bedtools ..."
echo "add chromosome ${chr} column, change delimiter to tab, change positions to integer ..."
echo "sort the file by individual and left positions ..."
#awk -v awkvar=${chr} -v OFS='\t' 'NR>1 {print awkvar,int($1),int($2),$3,$4}' ${bedfile} |\
# do not skip for first row since the local ancestry output has no header anymore
awk -v awkvar=${chr} -v OFS='\t' '{print awkvar,int($1),int($2),$3,$4}' ${bedfile} |\
sort -k5,5n -k2,2n > ${bedfile}.bed-tmp
echo "Done preparation for bedtools"

# step 2
vcffile=$2
echo "bedtools intersect to find overlapping positions of ${bedfile} and ${vcffile} ..."
bedtools intersect -a ${bedfile}.bed-tmp -b ${vcffile} > ${bedfile}.bed_intersect-tmp
echo "Done bedtools intersect"

# step 3 
echo "reformat the bedtools output to variants file and local ancestry file..."
# break the file by individual number, for each subfile
mkdir split_folder
awk -F"\t" '{print>"split_folder/"$5".txt"}' ${bedfile}.bed_intersect-tmp 
# this will create file with indivdiual number as filename
# sort for each file by end and extract line for local ancestry only
# and remove duplicate
for i in {0..599}; do sort -u -k3,3n split_folder/${i}.txt |\
awk -F"\t" '{print $4}' > split_folder/${i}.txt-anc; done
# and combine with individual as column number 
paste -d " " split_folder/{0..599}.txt-anc > ${bedfile}.reformat
# output variants list and calculated for genetic distance
sort -u -k3,3n split_folder/1.txt | awk -F"\t" '{print $3, $3/1000000}' > ${bedfile}.map
# remove the subfolder
#rm -r split_folder
echo "Done reformating!"


# step 5
#echo "removing all temporaty files..."
#rm ${bedfile}.*-tmp
#echo "Done!"
