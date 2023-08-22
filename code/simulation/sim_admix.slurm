#!/bin/bash -l        
#SBATCH --time=01:00:00 
#SBATCH --ntasks=1
#SBATCH --mem=20g
#SBATCH --tmp=20g
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=youremailaddress@umn.edu # report failed jobs to your email
#SBATCH --output=/admix_heritability/output/%j.out  #save the slurm output .out file

cd $SLURM_SUBMIT_DIR #Submit Directory 

module load R # msi needs this

model=${1} #HI or CGF
theta=${2} #0.1 0.2 0.5
gen=${3} #10, 20, 50, 100
P=${4} #0, 0.3, 0.6, 0.9
cov=${5} #pos or neg
seed=${6}


echo "simulate admix for model ${model} and theta ${theta} and gen ${gen} and P ${P} and ${cov} and seed ${seed}"
Rscript sim_admix.R -F 0.2 -l 1e3 -n 1e4 -p 0 -M ${model} --theta ${theta} -t ${gen} -P ${P} -C ${cov} --seed ${seed} 
