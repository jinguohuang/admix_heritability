#!/bin/bash -l        
#SBATCH --time=01:30:00
#SBATCH --ntasks=10
#SBATCH --mem=20g
#SBATCH --tmp=20g
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=klema030@umn.edu 
#SBATCH --output=/home/aazaidi/klema030/AdjustedHE/output/slurm/%j.out  #save the slurm output .out file at


model=${1} #HI or CGF
#theta=0.5 #0.1 0.2 0.5
#gen=20 #10, 20, 50, 100
P=${2} #0, 0.3, 0.6, 0.9
cov=${3} #pos or neg or zero
seed=${4}
t=${5}

./vg_GCTA_GRMvarX.sh ${model} ${P} ${cov} ${seed} ${t}

