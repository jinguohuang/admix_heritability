## HE estimation of $V_g$

**HE regression without any covariates**
* Script: **AdjHE_residual_noa.slurm**
* Usage: `sbatch AdjHE_residual_noa.slurm ${model} ${P} ${cov} ${seed} ${t} ${grm}`
* Example: `sbatch AdjHE_residual_noa.slurm HI 0 pos 1 0 grmld`

**HE regression with ancestry as a covariate**
* Script: **AdjHE_residual.slurm**
* Usage: `sbatch AdjHE_residual.slurm ${model} ${P} ${cov} ${seed} ${t} ${grm}`
* Example: `sbatch AdjHE_residual.slurm HI 0 pos 1 0 grmld`

**HE regression with GCTA**
* Requirement: PLINK 2.0, and GCTA version 1.94.1
* Scripts: **gcta_HE.slurm**, with scaled GRMs: **run_gcta_GRMld.slurm** and **run_gcta_GRMvarX.slurm**
* Usage:
  * **gcta_HE** `sbatch gcta_HE.slurm ${model} ${P} ${cov} ${seed} ${t}`
  * **GRMld:** `sbatch run_gcta_GRMld.slurm ${model} ${P} ${cov} ${seed} ${t}`
  * **GRMvarX:** `sbatch run_gcta_GRMvarX.slurm ${model} ${P} ${cov} ${seed} ${t}`
* Example: `sbatch run_gcta_GRMld.slurm CGF 0.9 neg 10 20`

**To compile outputs from scipts and GCTA without scaled GRMs**
* Script: **compile_data.sh**
* Usage: `./compile_data.sh ${model} ${P} ${cov}`
* Example: `.compile_data.sh HI 0 zero`

**To compile outputs from scripts and GCTA with scaled GRMs**
* Script: first run **compile_HE_GRM.sh** and then run **organize_GRM_summary.sh**
* Usage:
  * compile_HE_GRM.sh: `./compile_HE_GRM.sh ${model} ${P} ${cov} ${grm}`
  * organize_GRM_summary.sh: `./organize_GRM_summary.sh ${model} ${grm}`
* Example:
  * compile_HE_GRM.sh: `./compile_HE_GRM.sh CGF 0 pos varX`
  * organize_GRM_summary.sh: `./organize_GRM_summary.sh CGF ld` 
