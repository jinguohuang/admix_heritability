#!/bin/bash

#SBATCH -p moria
#SBATCH -c 1
#SBATCH -t 00:00:10
#SBATCH --mem=1g
#SBATCH --mail-type=END
#SBATCH --output=./log/slurm_%A_%a.out
#SBATCH --error=./log/slurm_%A_%a.out

bash ~/projects/admix_heritability/code/simulation/sim_vg.sh \
	$SLURM_ARRAY_TASK_ID \
	$1 \
	$2
