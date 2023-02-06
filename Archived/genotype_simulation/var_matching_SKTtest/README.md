### Variance matching and SKT test for simulated data and ASW
Process of variance calculation of local ancestry and SKT test:  
* Step1: simulate genotype and local ancestry  
  code: __var_match_loc_anc.py  var_match_loc_anc.sh__  
  usage: simulate 100 chromosome for 10deme and random mating model:
  ```
   for i in {1..100}; do ./var_match_loc_anc.sh $i 10deme; done 
   for i in {1..100}; do ./var_match_loc_anc.sh $i random; done
  ```
* Step2: reformat local ancestry to row as variants and column as individual:  
  code: __reformat_bedtools.py  reformat_loc_anc.sh   batch_reformat_loc_anc.sh__   
  usage: submit to ics to get reformatted local ancestry of 10deme and random mating model:
  ```
    for i in random 10deme; do qsub -A open batch_reformat_loc_anc.sh -F $i; done
  ```
* Step3: calculate variance and plot for ancestry distribution by locus  
  code: __plot_loc_anc.R  plot_loc_anc_var.sh__  
  usage: 
  ```
    ./plot_loc_anc_var.sh random
    ./plot_loc_anc_var.sh 10deme
  ```
* Step4: plot for ancestry distribution by person and calcualte variance  
  code: __plot_anc_dist.R  plot_ppl_distribution.sh__  
  usage: 
  ```
    ./plot_ppl_distribution.sh 10deme
    ./plot_ppl_distribution.sh random
  ```
* Step5: SKT test for simulated data and ASW  
  code: __plot_skt_test.R  skt_test_ASW.sh  plot_skt_dist_ASW.R__    
  usage:  
  ```
    Rscript plot_skt_test.R
    ./skt_test_ASW.sh
  ```
* Proportion matching update for 10deme model with ASW; reformat for local ancestry output of bedtools intersect.  
  code: __var_matching.ipynb__


#### Update variance matching and SKT test code
* Change the number of initial population size of admixed population to 1K in total, so it will be 100 per deme for 10 deme model.  
  code: __var_match_loc_anc_v2.sh__  
  Usage is the same as in step 1 above.

* Update the SKT test code, split by odd and even group and then sample with replacement or without replacement.  
  code: __skt_test_plot_v2.R__  It will plot for both simulated data and ASW.  
  usage:  
  ```
    Rscript skt_test_plot_v2.R
  ```
  
  
