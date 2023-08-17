## Simulating complex traits in admixed populations
* Script: **sim_admix.R, main_admix.R**
* Requirement: R version 4.2.3
* Usage: ```Rscript sim_admix.R -F ${FST} -l ${nloci} -n ${nindv} -M ${model} --theta ${theta} -t ${gen} -P ${P} -C ${cov} --seed ${seed} ```

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


* Example: To simulate a trait under divergent selection with Fst 0.2 between source populations, causal loci 1e3, with 1e4 effective population size, under hybrid isolation (HI) model, with 50% admixture proportion at t=0, with assortative mating strength P=0.9 and track for 20 generations for a replicate under seed 1: ```Rscript sim_admix.R -F 0.2 -l 1e3 -n 1e4 -p 0 -M HI --theta 0.5 -t 20 -P 0.9 -C pos --seed 1 ```
* Output: A table records the essential summary statistics of the simulated admixed population for each generation, such as mean and variance of ancestry, genetic variance, etc.
* Output (optional, required files for GREML estimation): output genotype/local ancestry, PRS, phenotype to PLINK acceptable format: .dosage, .fam, .pheno, .ganc (global ancestry as a covariate).
