## Estimating ${V}_g$ with GREML 
* Script: **vg_GCTA.sh**
* Requirement: PLINK 2.0, and GCTA version 1.94.1
* Usage: ```./vg_GCTA.sh ${model} ${theta} ${gen} ${P} ${cov} ${seed} ${t}```
* Example: To estimate the genetic variance of the complex trait in the admixed population at generation 10 simulated above using GREML implemented in GCTA: ```./vg_GCTA.sh HI 0.5 20 0.9 pos 1 10```
* Output: genetic relationship matrix of the population, reml file of estimation results with and without global ancestry as covariate.

* Script: **GRMvarX.R**
* Usage: ```Rscript GRMvarX.R ${filename} ${plinkdir} ${grmdir}```
* Example: To feed GCTA with GRM with updated scaling: variance of genotype, instead of the standard scaling of GREML 2f(1-f), substitute this line with command to construct grm with gcta in **vg_GCTA.sh**
  
* Script: **GRMld.R**
* Usage: ```Rscript GRMld.R ${filename} ${plinkdir} ${grmdir}```
* Example: To feed GCTA with GRM with updated scaling: LD matrix, instead of the standard scaling of GREML 2f(1-f), substitute this line with command to construct grm with gcta in **vg_GCTA.sh**
