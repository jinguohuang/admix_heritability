# Admix heritability code
This repository contains codes used to carry out the analyses and create the main figures in "Interpreting SNP heritability in admixed populations".

## Simulating complex traits in admixed populations
* Script: **sim_admix.R, main_admix.R**
* Requirement: R version 4.2.3
* Usage: ```Rscript sim_admix.R -F ${FST} -l ${nloci} -n ${nindv} -M ${model} --theta ${theta} -t ${gen} -P ${P} -C ${cov} --seed ${seed} ```
* Example: To simulate a trait under divergent selection with Fst 0.2 between source populations, causal loci 1e3, with 1e4 effective population size, under hybrid isolation (HI) model, with 50% admixture proportion at t=0, with assortative mating strength P=0.9 and track for 20 generations for a replicate under seed 1: ```Rscript sim_admix.R -F 0.2 -l 1e3 -n 1e4 -p 0 -M HI --theta 0.5 -t 20 -P 0.9 -C pos --seed 1 ```
* Output: A table records the essential information of the simulated admixed population for each generation, such as mean and variance of ancestry, genetic variance, etc.
* Output (optional, required files for GREML estimation): output genotype/local ancestry, PRS, phenotype to PLINK acceptable format: .dosage, .fam, .pheno, .ganc (global ancestry as a covariate).

## Estimating ${V}_g$ with GREML 
* Script: **vg_GCTA.sh**
* Requirement: PLINK 2.0, and GCTA version 1.94.1
* Usage: ```./vg_GCTA.sh ${model} ${theta} ${gen} ${P} ${cov} ${seed} ${t}```
* Example: To estimate the genetic variance of the complex trait in the admixed population at generation 10 simulated above using GREML implemented in GCTA: ```./vg_GCTA.sh HI 0.5 20 0.9 pos 1 10```
* Output: genetic relationship matrix of the population, reml file of estimation results with and without global ancestry as covariate.

## Decomposing ${V}_g$ in traits in African Americans
*

## Plot

* Fig2: ${\theta}$


* Fig3: ${V}_g$


* Fig4: $\hat{V}_g$


* Fig5: $\hat{V}_{\gamma}$


* Fig6: Decomposing variance explained

