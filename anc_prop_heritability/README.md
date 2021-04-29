### Simulate genotype with different admixture proportion
* Simulate admixed population under different admixture proportion with script: __sim_anc_prop.py__  
  Could change input like: ouput prefix (-o), chromosome number (random seed, -c), admixture proportion (-p),   chromosome lenghth (-L), admixed population sample size (-adm), reference population sample size (-ref),    etc.
  could be used like:
  ```
  python sim_anc_prop.py -o simADM -c 1 -p 0.1 -L 20000000 -adm 20000 -ref 2000 
  ```
  this will simulate 10K admixed people and 1K reference parental people with 0.1 admixture proportion and 200M  chromosome.
* Simulate admixed population with the script above, with input variables of admixture proportion and     chromosome number: __sim_anc_prop.sh__ 
  could be used like:
  ```
  ./sim_anc_prop.sh 0.1 1
  ```
  this will simulate the same population described above.
* After genotype simulation, this script __anc_estimate_plot_1.sh__ will merge multiple chromosomes vcf files and convert to PLINK format and do ancestry estimate with ADMIXTURE and plot the ancestry estimation with Rscript __plot_admixture_dist.R__  
  could be submitted to ics and run for admixture proportion 0.1~0.5:
  ```
  for i in 0.1 0.2 0.3 0.4 0.5; do qsub -A open anc_estimate_plot_1.sh -F $i; done
  ```

### Simulate phenotype with causal variants of Fst=0.2
* Uniformly sampled causal variants 1 per 100kb across 10 chromomes, among variants with parental populations 0.15<Fst<0.25, then apply (Schoech et al. 2019; Zaidi et al. Elife 2020) model to simulate effect sizes for causal variants. This script __CausalVariant_freq_fst.R__ needs input files: frequency file of variants in admixed population and Fst file of parental groups. Output of this script including CausalVariantScore file that used for calculating genetic score, and CausalVariantID that used for extracting causal variants from genotype data.  
  Could be used like:
  ```
  Rscript CausalVariant_freq_fst.R ADM_0.2prop.frq ref_0.2prop_updateid.fst
  ```
  
* After genetic scores got calculated, the .profile file can be used as input for the script __sim_pheno_V1_ancprop.R__ to simulate phenotype data with heritability = 0.8.  
  Could be used like:
  ```
  Rscript sim_pheno_V1_ancprop.R ADM_0.2prop_2kCausal.profile
  ```
  


