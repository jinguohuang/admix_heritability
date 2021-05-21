### Admixture LD curve
Plot the genetic distance (cM) vs local ancestry correlation of simulated data.

* Step 1: Simulate genotype: 200Mb, 400hap, initial pop size=1K, random/10deme.  
  code: __sim_AA_geno_lanc.py__  
  usage:  
  ```
  python sim_AA_geno_lanc.py -o AA_200Mb_400hap -c 1 -model random -L 200000000
  python sim_AA_geno_lanc.py -o AA_200Mb_400hap -c 1 -model 10deme -L 200000000
  ```
* Step 2: Reformat for local ancestry and calculate genetic distance from physical distance (1cM=1Mb).  
  This is updated reformat to match the output of the simulated local ancestry above.
  code: __reformat_loc_anc.sh  batch_reformat_loc_anc.sh__  
  usage:  
  ```
  qsub -A open batch_reformat_loc_anc.sh -F 10deme
  qsub -A open batch_reformat_loc_anc.sh -F random
  ```
* Step 3: Sample 50cM from 200cM chromosome and plot for genetic distance vs local ancestry correlation.  
  Sample 100~150cM, treat every cM as a group.  
  Sample 10 SNPs from each group. Calculate genetic distance and correlation of group 1 with every group.  
  code: __sample_500pos.R  plot_ALDcurve.R__   
  usage:  
  ```
  Rscript sample_500pos.R loc_anc_AA_200Mb_400hap_10deme_chr1.vcf.txt.map loc_anc_AA_200Mb_400hap_10deme_chr1.vcf.txt.reformat_adm
  Rscript sample_500pos.R loc_anc_AA_200Mb_400hap_random_chr1.vcf.txt.map loc_anc_AA_200Mb_400hap_random_chr1.vcf.txt.reformat_adm 
  Rscript plot_ALDcurve.R map_sample500_10deme.txt lanc_sample500_10deme.txt 
  Rscript plot_ALDcurve.R map_sample500_random.txt lanc_sample500_random.txt
  ```
  


