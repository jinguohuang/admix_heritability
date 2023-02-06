##!/bin/bash
# script to reformat local ancestry 
# step 1: reformat the local ancestry to bed file for bedtools
# step 2: bedtools intersect to find overlapping variant and their local ancestry
# step 3: remove duplicates
# step 4: python to unmelt the output to 2 files: variants x individual format and the variants list
# step 5: remove all temporary files

# input1 is the local ancestry file from simulation
# input2 is the vcf file from simulation
# input3 is the python script to reformat the local ancestry
bedfile=$1
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
awk -v awkvar=${chr} -v OFS='\t' 'NR>1 {print awkvar,int($1),int($2),$3,$4}' ${bedfile} |\
sort -k5,5n -k2,2n > ${bedfile}.bed-tmp
echo "Done preparation for bedtools"

# step 2
vcffile=$2
echo "bedtools intersect to find overlapping positions of ${bedfile} and ${vcffile} ..."
bedtools intersect -a ${bedfile}.bed-tmp -b ${vcffile} > ${bedfile}.bed_intersect-tmp
echo "Done bedtools intersect"

# step 3
echo "remove duplicate..."
uniq -u ${bedfile}.bed_intersect-tmp > ${bedfile}.bed_intersect_uniq-tmp
dup=$(diff -y --suppress-common-lines ${bedfile}.bed_intersect_uniq-tmp ${bedfile}.bed_intersect-tmp |\
grep '^' | wc -l) 
echo "Done remove duplicate, ${dup} duplicates have been removed"

# step 4 
pyfile=$3
echo "reformat the bedtools output to variants file and local ancestry file..."
python ${pyfile} ${bedfile}.bed_intersect_uniq-tmp
echo "Done reformating!"

# step 5
echo "removing all temporaty files..."
rm ${bedfile}.*-tmp
echo "Done!"
