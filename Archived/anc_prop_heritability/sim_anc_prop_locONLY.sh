#!/bin/bash
#PBS -l nodes=1:ppn=20 ## Requests 1 processor on 20 node
#PBS -l walltime=01:30:00 ## Requests 24 hours of walltime
#PBS -l pmem=4gb ## Requests 8 gigabytes of memory per process
#PBS A ## Specifies the allocation. Use A open for open queue
#PBS -j oe ## Requests that regular output and terminal output go to the same file
## The following is the body of the script. By default PBS scripts execute in your home directory, not the
## directory from which they were submitted. The following line places you in the directory from which the job
## was submitted.
cd $PBS_O_WORKDIR
prop=$1
chrom=$2

python sim_anc_prop_locONLY.py -o simADM -c ${chrom} -p ${prop} -L 20000000 -adm 20000 -ref 2000

