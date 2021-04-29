
### Simulate admixed population - African American
Folder: simulate_AdmPop/

* V1: __sim_AfrAmr_20200824.py__  
With ASN admixture, 100Mb chromosome 
* V2: __sim_AfrAmr_20200830.py__  
Delete ASN events, 100Mb chromosome  
* V3: __anc_estimate_10Mb.sh sim_AA_v3_10Mb.py__  
Seeds {1..22} as different chromosomes, 22 10Mb chromosomes
* 10demes: __AA_sim_10demes.ipynb__  
Simulate 10 demes with different admixture proportion,works fine with certain random seeds in 10MB, but tends to error if 100MB. 
Found out the error was caused by the number of admixed individuals, deleted this file because git commit error.
* 10demes/11demes: __AA_sim_10demes_1.ipynb anc_estimate.sh plot_admixture.R__  
Simulate 10 demes with different admixture proportion and different population size. Solved error in last script. 


### Get local ancestry of simulated individuals with msprime
Folder: local_anc_trail1/

* Initial test the code: __Local_anc_msprime.ipynb__    
  Function credit: https://gist.github.com/rwaples/10336129f75239465279dbfe163bf3c1  
  Test with 11 demes model, 100Mb chromosome, test with AFR,EUR,ADM individuals respectively at 20, 5000 generations
* Test the code with simpler model: __Local_anc_msprime_2.ipynb anc_estimate_2.sh local_to_global.py__  
  Test with 1 admixed deme model, 100Mb chromosome, remove growth rate and migration rate before 16 generations, 200 haplotypes for AFR, EUR, ADM.   
  Export their genotype and estimate global ancestry with ADMIXTURE.  
  Calculate their global ancestry (True) from the local ancestry and compare with the predicted.
* Test the code with 10 demes model: __Local_anc_msprime_3.ipynb anc_estimate_3.sh local_to_global_anc_10deme.py local_anc_10deme_plot.py__  
  Test with 10 admixed (decile, 0.05-0.95) deme model, 100Mb chromosome, remove growth rate and migration rate before 16 generations, 200 haplotypes for AFR, EUR,     2000 ADM.  
  Export their genotype and estimate global ancestry with ADMIXTURE.  
  Since it took too long to calculate local anc for all individuals, pick sample (AFR,EUR,AFR), Calculate their global ancestry (True) from the local ancestry for     each and plot their local anc at generation 20 and 2030 (10generations after AFR-EUR split).
* Output local ancestry of 400 hap model below: __local_anc_single.py local_anc_single.sh__  
  Ouput local ancestry of each person in 400 hap model in dictionary format (.py file).
  The shell file will submit job to ics with the command like: 
  ```
  for i in {501..600}; do qsub -A open -l feature=rhel7 local_anc_single.sh -F $i; done
  ```

### Test with different migration rate for better model fit
Folder: mig_rate_model/

* 10 demes model with different migration rate vs random mating model: __mig_rate_random.ipynb AA_10deme_mig_sim.py__  
  Simulate 10 demes model(with pop structure) 10Mb, time of admixture: 16 gen ago, AFR,EUR: 100 hap ADM: 400 hap in total. Average anc: 0.76525. No growth rate,       migration rate before 16 gen. 10 demes follow decile admixture proportion(0.05~0.95). Allow neighboring mig_rate between demes 0.025, 0.05, 0.075, 0.1, 0.15,0.2.  
  Simulate admixed pipulation (without pop structure) 10 MB, with same parameters with 10 demes model except for no migration rate.  
  Export their genotype and estimate global ancestry with ADMIXUTRE. Plot the density plot to compare simulated ancestry distribution with ASW haplotype      distribution. And KS test to see if distribution same. __plot_admixture_dist_density_kstest.R anc_estimate_plot.sh__  

* Debug the better fit simulation model above:__sim_AA_10deme_migrate.py sim_AA_random_mate.py__  
  -Adjust the migration matrix edge (double the end migration rate only)  
  -Change to diploid output  
  -Change to 10Mb x 10.  
  -Random mate model with equal average admixture proportion with 10 deme model.  
  Run the above script and merge 10 chromosomes, plot and KS test:__anc_estimate_10Mb_migrate.sh anc_estimate_10Mb_random.sh__  
  Plot script: __plot_admixture_dist_density_kstest_1.R__  
  ASW AFRICAN ancestry distribution data: __ASW_dip_anc1.txt ASW_ind_anc1.txt__   
  

* Follow the better fit simulation model above and simulate 10000 haplotypes: __sim_AA_10deme_10000hap.py sim_AA_random_mate_10000hap.py__  
  Run the above script and merge 10 chromosomes, plot and KS test: __anc_estimate_10Mb_10000hap_migrate.sh anc_estimate_10Mb_random_10000hap.sh__  
  Plot script: __plot_admixture_dist_density_kstest_2.R__  
  
### New method to Get local ancestry of simulated individuals with msprime 
Folder: local_anc_new/

* Reference code: https://gist.github.com/gtsambos/3165c24aebb50d2ca4b743b25673c413
* Test with new local ancestry inference method, test with my 400hap, 10,000hap random mating model, and compare the local to global vs ADMIXTURE extimated global ancestry results __Loc_anc.ipynb__  
* Test with the new method in script __loc_anc_10M_20000hap.py__ to get the local ancestry of 10Mb, 10K admixed people and 50 people per reference population, can finish within an hour.
* Process the simulated genotype and do ancestry estimation __anc_estimate_10Mb_random_20000hap.sh__ 

### New method to get Local ancestry variance
Folder: var_local_anc_new/

* Convert the local ancestry to format: each row as variants, each comlumn as individual (100000000 * 400) and plot for AFR ancestry percentage at each locus, and plot local to global percentage against admixture global for 10deme, random mating 400 hap model: __var_local_anc.ipynb__
* Simulate 400hap 100MB genotype under random mating and 10 deme model and get local ancestry for each person: __loc_anc_400hap_100M.py__, __loc_anc_400hap_100M_10deme.py__ 
* Process the simulated genotype and do ancestry estimation with either random or 10deme as input for processing the corresponding one: __anc_estimate_100M_400hap.sh__  
  Could be used like:
  ```
  ./anc_estimate_100M_400hap.sh random
  ```
* Convert the local ancestry to format: each row as variants, each comlumn as individual (100000000 * 400), output as a file saving the local ancestry: __reformat_loc_anc_10deme.py__ , __reformat_loc_anc_random.py__
* Plot the the distribution of AFR ancestry calculated from the output of last step in a Manhattan plot, and calculate the mean and variance of AFR ancestry: __plot_loc_anc_10deme.R__ , __plot_loc_anc_random.R__ 






