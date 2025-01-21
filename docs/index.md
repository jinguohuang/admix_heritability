# Admix heritability code
This repository contains codes used to carry out the analyses and create the main figures in "Interpreting SNP heritability in admixed populations".

## Simulating complex traits in admixed populations
* Script: **sim_admix.R, main_admix.R**
* Requirement: R version 4.2.3
* Usage: ```Rscript sim_admix.R -F ${FST} -l ${nloci} -n ${nindv} -M ${model} --theta ${theta} -t ${gen} -P ${P} -C ${cov} --seed ${seed} ```
* Example: To simulate a trait under divergent selection with Fst 0.2 between source populations, causal loci 1e3, with 1e4 effective population size, under hybrid isolation (HI) model, with 50% admixture proportion at t=0, with assortative mating strength P=0.9 and track for 20 generations for a replicate under seed 1: ```Rscript sim_admix.R -F 0.2 -l 1e3 -n 1e4 -p 0 -M HI --theta 0.5 -t 20 -P 0.9 -C pos --seed 1 ```
* Output: A table records the essential summary statistics of the simulated admixed population for each generation, such as mean and variance of ancestry, genetic variance, etc.
* Output (optional, required files for GREML estimation): output genotype/local ancestry, PRS, phenotype to PLINK acceptable format: .dosage, .fam, .pheno, .ganc (global ancestry as a covariate).

## Estimating ${V}_g$ with GREML 
* Script: **vg_GCTA.sh**
* Requirement: PLINK 2.0, and GCTA version 1.94.1
* Usage: ```./vg_GCTA.sh ${model} ${theta} ${gen} ${P} ${cov} ${seed} ${t}```
* Example: To estimate the genetic variance of the complex trait in the admixed population at generation 10 simulated above using GREML implemented in GCTA: ```./vg_GCTA.sh HI 0.5 20 0.9 pos 1 10```
* Output: genetic relationship matrix of the population, reml file of estimation results with and without global ancestry as covariate.

## Estimating ${V}_g$ with HE regression 
* Script: **AdjHE_residual_noa.slurm, AdjHE_residual.slurm, gcta_HE.slurm, run_gcta_GRMld.slurm, run_gcta_GRMvarX.slurm**
* Requirement: R version 4.2.3, PLINK 2.0, and GCTA version 1.94.1
* Usage: see README.md in [**HE_simulation**](https://github.com/jinguohuang/admix_heritability/tree/master/code/HE_simulation)
* Example: To estimate the genetic variance of the complex trait in the admixed population at generation 10 simulated above using HE regression without ancestry as a covariate implemented in R: ```sbatch AdjHE_residual_noa.slurm HI 0 pos 1 0```
* Output: .txt file of estimation results with global ancestry as covariate.

## Decomposing ${V}_g$ in traits in African Americans
* Script: **LDclump_hg38.sh, Allelef1f2_hg38.sh, vgCI_estimate.R**
* Requirement: R version 4.2.3, PLINK 1.9, PLINK 2.0
* Usage: ```./LDclump_hg38.sh ${trait}``` and ```vgCI_estimate.R``` 
* Example: To process GWAS Catalog downloaded summary statistics for trait HDL: LD clumping, calculate allele frequency of associated SNPs in CEU and YRI, calculate the proportion of variance explained: ```LDclump_hg38.sh HDL```, ```vgCI_estimate.R``` (input is Table S1) 
* Output: the proportion of variance explained for the trait with 95% confidence intervals ```vgpCI.txt```.

## Plot

* Fig2: ${\theta}$: ```Rscript plot_theta.R ```

* Fig3: ${V}_g$: ```Rscript plot_vg.R```

* Fig4: $\hat{V}_g^{GREML}$: ```Rscript plot_GREML_vg.R```

* Fig5: $\hat{V}_g^{HE}$: ```Rscript plot_HE_vg_CI.R```

* Fig7: $\hat{V}_{\gamma}$: ```Rscript plot_GREML_vgamma.R```

* Fig8: Decomposing variance explained: ```Rscript plot_vgCI.R```
