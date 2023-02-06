#!/bin/bash
#PBS -l nodes=1:ppn=2 ## Requests 1 processor on 2 node
#PBS -l walltime=0:00:10 ## Requests hours of walltime
#PBS -l pmem=4gb ## Requests 4 gigabytes of memory per process
#PBS A open ## Specifies the allocation. Use A open for open queue
#PBS -j oe ## Requests that regular output and terminal output go to the same file
## The following is the body of the script. By default PBS scripts execute in your home directory, not the
## directory from which they were submitted. The following line places you in the directory from which the job
## was submitted.
#cd $PBS_O_WORKDIR

# run with chromosome number and model


python var_match_loc_anc.py -o simAA_1Mb_1Kppl -c $1 -model $2 -L 1000000 -ref 200 -adm 2000 --nadmix 1000






