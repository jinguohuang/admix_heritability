#!/bin/bash
#PBS -l nodes=10:ppn=1 ## Requests 1 processor on 20 node
#PBS -l walltime=5:00:00 ## Requests 24 hours of walltime
#PBS -l pmem=50gb ## Requests 8 gigabytes of memory per process
#PBS A ## Specifies the allocation. Use A open for open queue
#PBS -j oe ## Requests that regular output and terminal output go to the same file
## The following is the body of the script. By default PBS scripts execute in your home directory, not the
## directory from which they were submitted. The following line places you in the directory from which the job
## was submitted.
#cd $PBS_O_WORKDIR
# to get estimated vgb with different vgb

vgb=$1
#theta=$1

echo "start simulation for vgb${vgb}"
for theta in $(seq 0.1 0.1 0.5)
do
	echo "start for prop${theta} vgb${vgb}"
	echo "Simulate!"
	Rscript simpleSim_addenv.R --vgb ${vgb} --theta ${theta} 
	echo "Done simulation"
	filename1=admix_prop${theta}_vgb${vgb}
	filename2=admix_geno_prop${theta}_vgb${vgb}
	for filename in ${filename1} ${filename2}
	do
		echo "the filename is ${filename}"
		echo "Convert tped to bfile"
		plink \
		--tfile ${filename} \
		--make-bed \
		--out ${filename}
		echo " start calculating grm "
		gcta64 \
		--bfile ${filename} \
		--make-grm-bin \
		--thread-num 20 \
		--out ${filename}
		echo " start estimate vgb"
		gcta64 \
		--reml \
		--grm ${filename} \
		--pheno pheno_prop${theta}_vgb${vgb}.txt \
		--thread-num 20 \
		--out ${filename}
	done
	echo "done estimate for prop${theta} vgb${vgb}"
	echo "start reformat for prop${theta} vgb${vgb}"
	Rscript simVSest_reformat.R ${theta} ${vgb}
	echo "Done reformat for vgb${vgb}"
done
echo "Done simulation for vgb${vgb}"
echo "Start organise results"
head -1 simVSest_prop0.1_vgb${vgb}.txt > simVSest_vgb${vgb}.txt
for i in $(seq 0.1 0.1 0.5)
do 
	tail -1 simVSest_prop${i}_vgb${vgb}.txt >> simVSest_vgb${vgb}.txt 
done



