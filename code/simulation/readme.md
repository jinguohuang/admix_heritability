## Simulating complex traits in admixed populations
* Script: **sim_admix.R, main_admix.R**
* Requirement: R version 4.2.3
* Usage:
  ```
  Rscript sim_admix.R -F ${FST} -l ${nloci} -n ${nindv} -M ${model} --theta ${theta} -t ${gen} -P ${P} -C ${cov} --seed ${seed}
  ```

  | Flag | Description | Default 
  | :---: | :---  | :---: |
  | `--fst, -F` | Desired Fst between parental groups at causal loci | 0.2 | 
  | `--nloci, -l` | Number of loci | 1e3 |
  | `--psel, -p` | Selection strength | 0 |
  | `-n` | Sample size of population | 1e3 |
  | `--seed` | Random seed | 120 |
  | `--pganc, -P` | Strength of ancestry-based assortative mating strength | 0.9 |
  | `--gen, -t` | Number of generations (time) since admixture | 1 |
  | `--theta` | Total admixture proportion from pop1, range 0 to 1 | 0.5 |
  | `--model, -M` | Admixture model: HI or CGF  | HI |
  | `--cov, -C` | Covariance sign of trait architecture: pos or neg, pos <br />for divergent selection, neg for stabilizing selection | pos |


* Example: To simulate a trait under divergent selection with Fst 0.2 between source populations, causal loci 1e3, with 1e4 effective population size, under hybrid isolation (HI) model, with 50% admixture proportion at t=0, with assortative mating strength P=0.9 and track for 20 generations for a replicate under seed 1:
  ```
  Rscript sim_admix.R -F 0.2 -l 1e3 -n 1e4 -p 0 -M HI --theta 0.5 -t 20 -P 0.9 -C pos --seed 1
  ```

* Output: A table records the essential summary statistics of the simulated admixed population for each generation, such as mean and variance of ancestry, genetic variance, etc. The output columns are explained as follows:
  | Colomn | Description  
  | :---: | :---  |
  | _t_ | Time (generation) since admixture (t)|
  | _mean.theta_ | Mean of ancestry in the admixed population ($E_{\theta}$)|
  | _var.theta_ | Variance of ancestry in the admixed population ($V_{\theta}$) |
  | _var.prs.geno_ | Genetic variance due to genotype ($V_g$) |
  | _var.prs.lanc_ | Genetic variance due to local ancestry ($V_{\gamma}$) |
  | _var.pheno_ | Variance in phenotype |
  | _vg.term1_ | (1.1) of $V_g$: genetic variance comes from source populations |
  | _vg.term2, vg.term3_  | (1.2), (1.3) of $V_g$: genetic variance due to allele frequency difference |
  | _vg.term4_ | (1.4) of $V_g$: genetic variance come from LD |
  | _vg.sum_ | Expected $V_g$: (1.1)+(1.2)+(1.3)+(1.4) |
  | _cov.ganc.prsgeno, cor.ganc.prsgeno_ | Covariance, correlation of ancestry and PRS|
  | _cov.ganc.prslanc, cor.ganc.prslanc_ | Covariance, correlation of ancestry and PRS due to local ancestry |
  | _cov.ganc.pheno, cor.ganc.pheno_ | Covariance, correlation of ancestry and phenotype |
  | _theta_ | Proportion of total admixture proportion ${\theta}$|
  | _P_ | Strength of ancestry-based assortative mating |
  | _seed_ | Random seed for replication |
  | _vbetween.genic_ | The genic genetic variance between source populations |
  | _vbetween.obs_ | The genetic variance between source populations calculated with genetic values |
  | _cov_ | Sign of the LD contribution of two simulated trait architectures |
  | _model_ | Admixture model used: hybrid isolation (HI) or continuous gene flow (CGF) |
  | _m1_ | Migration rate at each generation for CGF model |
  | _prs.diff_ | Mean PRS difference in source populations |

* Output (optional, required files for GREML estimation): output genotype/local ancestry, PRS, phenotype to PLINK format.
  | File extension | Description  
  | :---: | :---  |
  | _.dosage_ | Genotype or local ancestry matrix |
  | _.tfam_ | Invidividual info |
  | _.pheno_ | Phenotype |
  | _.prs_ | Genetic value or PRS |
  | _.ganc_ | Global ancestry of individuals |

* To replicate the experiments in this paper, submit jobs with the script **sim_admix.slurm** to your computing cluster as follows. Our job scheduler is Slurm, if you use PBS, you can convert it with this [cheatsheet](https://www.msi.umn.edu/slurm/pbs-conversion).
  ```
  # HI
  for s in {1..10}; do for i in 0.1 0.2 0.5; do for j in 0 0.3 0.6 0.9; do for k in pos neg; do for g in 100; do sbatch sim_admix.slurm HI $i $g $j $k $s; done; done; done; done; done
  # CGF
  for s in {1..10}; do for i in 0.1 0.2 0.5; do for j in 0 0.3 0.6 0.9; do for k in pos neg; do for g in 10 20 50 100; do sbatch sim_admix.slurm CGF $i $g $j $k $s; done; done; done; done; done
  ```


  
