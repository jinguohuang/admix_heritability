#!/bin/bash
#PBS -l nodes=1:ppn=20 ## Requests 1 processor on 20 node
#PBS -l walltime=24:00:00 ## Requests 24 hours of walltime
#PBS -l feature=rhel7  ##job will run on a RHEL 7 node
#PBS -l pmem=8gb ## Requests 8 gigabytes of memory per process
#PBS –A open ## Specifies the allocation. Use –A open for open queue
#PBS -j oe ## Requests that regular output and terminal output go to the same file
## The following is the body of the script. By default PBS scripts execute in your home directory, not the
## directory from which they were submitted. The following line places you in the directory from which the job
## was submitted.
cd $PBS_O_WORKDIR
# for run local ancestry for each sample
i=$1
python local_anc_single.py --outpre local_anc --sample $i

