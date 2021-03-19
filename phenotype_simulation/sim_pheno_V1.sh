#!/bin/bash
# Apply the Schoech et al. 2019; Zaidi et al. Elife 2020 model 
# to effect size to get phenotype data
# so effect size is frequency dependent

#1,2 Simple model with 10,000 ppl random mating model
#Simulating genotype
./anc_estimate_10Mb_random_10000ppl.sh 
# killed at ancestry estimation

# submit to ics to do ancestry estimate
qsub -A open -l feature=rhel7 anc_est.pbs

# plot the ancestry estimate result
Rscript plot_admixture_dist_density_kstest_2.R AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05_pruned.2.Q AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05_pruned.fam ASW_dip_anc1.txt

#3  Extract causal variants
#4 Simulate effect sizes for these variants 
# so they are dependent to frequency
# calculate frequency
plink --bfile AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05 --freq --out AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05
# sample causal variants
Rscript CausalVariantScore_V1.R 
# extract their genotype
plink --bfile AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05 --extract CausalVariant1000_v1.txt --make-bed --out AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05_1kSNP
#1000 variants and 10100 people pass filters and QC.

#5 Calculate genetic score
plink --bfile AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05_1kSNP --score CausalVariantScore_v1.txt --out AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05_1kSNP
#--score: Results written to AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05_1kSNP.profile .

# get genetic score of admixed ppl only
tail -10000 AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05_1kSNP.profile > sim_10Kppl_geneticscore_V1.txt

# remove 100 non-admixed ppl
head -100 AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05_extract.fam | awk '{print $1,$2}' > remove_100ppl.txt  
plink --bfile AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05_1kSNP --remove remove_100ppl.txt --make-bed --out AA_random_10Kppl_1kSNP
# 1000 variants and 10000 people pass filters and QC.

#6,7,8 Simulate phenotype data
Rscript sim_pheno_V1.R
# get pheno data sim_pheno_10Kppl_random_2.txt

#9 Estimate heritability
# construct grm
gcta64 --bfile AA_random_10Kppl_1kSNP --autosome --make-grm-bin --out AA_random_10Kppl_1kSNP
# estimate heritability
gcta64 --reml --grm-bin AA_random_10Kppl_1kSNP --pheno sim_pheno_10Kppl_random_V1.txt --out GCTAherit_sim_pheno_10Kppl_random_1kSNP
# V(G)/Vp	0.831198	0.006775



