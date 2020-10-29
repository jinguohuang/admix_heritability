# admix_simulation

This repo is for simulating the trait correlation under ancestry stratification

### Simulate admixed population - African American, Cape Verdean

* V1: __sim_AfrAmr_20200824.py__
With ASN admixture, 100Mb chromosome 
* V2: __sim_AfrAmr_20200830.py__
Delete ASN events, 100Mb chromosome 
* V3: __anc_estimate_10Mb.sh sim_AA_v3_10Mb.py__
Seeds {1..22} as different chromosomes, 22 10Mb chromosomes
* 10demes: __AA_sim_10demes.ipynb__
Simulate 10 demes with different admixture proportion,works fine with certain random seeds in 10MB, but tends to error if 100MB
* 10demes/11demes: __AA_sim_10demes_1.ipynb anc_estimate.sh plot_admixture.R__
Simulate 10 demes with different admixture proportion and different population size. Solved error in last script. 
