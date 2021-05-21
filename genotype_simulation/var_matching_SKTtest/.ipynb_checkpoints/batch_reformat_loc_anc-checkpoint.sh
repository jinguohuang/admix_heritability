#!/bin/bash
#PBS -l nodes=1:ppn=2 ## Requests 1 processor on 20 node
#PBS -l walltime=6:00:00 ## Requests 24 hours of walltime
#PBS -l pmem=4gb ## Requests 8 gigabytes of memory per process
#PBS A ## Specifies the allocation. Use A open for open queue
#PBS -j oe ## Requests that regular output and terminal output go to the same file
## The following is the body of the script. By default PBS scripts execute in your home directory, not the
## directory from which they were submitted. The following line places you in the directory from which the job
## was submitted.
cd $PBS_O_WORKDIR
##For processing and ancestry estimation of msprime generated genotype data
##To use this script, one need to add either random or 10deme in the command

model=$1
echo "Reformating for ${model} model ..."

for chr in {1..100}
do
	echo "Processing chr${chr} of ${model}"
	vcffile=simAA_1Mb_1Kppl_${model}_chr${chr}.vcf
	bedfile=loc_anc_${vcffile}.txt
	./reformat_loc_anc.sh ${bedfile} ${vcffile} reformat_bedtools.py
done


