#!/bin/bash -l        
#SBATCH --time=01:00:00
#SBATCH --ntasks=10
#SBATCH --mem=20g
#SBATCH --tmp=20g
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=klema030@umn.edu 
#SBATCH --output=/home/aazaidi/klema030/AdjustedHE/output/slurm/%j.out  #save the slurm output .out file at


module load R/4.1.0

model=${1} #HI or CGF
#theta=${2} #0.1 0.2 0.5
#gen=${3} #10, 20, 50, 100
P=${2} #0, 0.3, 0.6, 0.9
cov=${3} #pos or neg
seed=${4} #1-10
t=${5} #0-20
grm=${6} #GRMld, GRMvarX, GRM

grmdir=/home/aazaidi/klema030/AdjustedHE/output/grm/gcta/${grm}/${model}

echo "simulate adjusted HE for model ${model} and theta 0.5 and gen 20 and P ${P} and ${cov} and seed ${seed} and t ${t}"
Rscript /home/aazaidi/klema030/AdjustedHE/AdjHE_residuals.R ${model} ${P} ${cov} ${seed} ${t} ${grm} ${grmdir}
