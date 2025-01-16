#!/bin/bash -l        
#SBATCH --time=00:00:20
#SBATCH --ntasks=1
#SBATCH --mem=2g
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=klema030@umn.edu 
#SBATCH --output=/home/aazaidi/klema030/AdjustedHE/output/slurm/%j.out  #save the slurm output .out file at


model=${1} #HI or CGF
#theta=${2} #0.1 0.2 0.5
#gen=${3} #10, 20, 50, 100
P=${2} #0, 0.3, 0.6, 0.9
cov=${3} #pos or neg or zero
seed=${4} #1-10
t=${5} #0-20

echo "simulate HE for model ${model} and theta 0.5 and gen 20 and P ${P} and ${cov} and seed ${seed} and t ${t}"
./gcta_HE.sh ${model} ${P} ${cov} ${seed} ${t}