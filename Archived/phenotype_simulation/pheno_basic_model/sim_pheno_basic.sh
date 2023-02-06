#!/bin/bash
# Apply the basic model to simulated genotype data 
# to get phenotype data

#1,2 Simple model with 10,000 ppl random mating model
#Simulating genotype
./anc_estimate_10Mb_random_10000ppl.sh 
# killed at ancestry estimation

# submit to ics to do ancestry estimate
qsub -A open -l feature=rhel7 anc_est.pbs

# plot the ancestry estimate result
Rscript plot_admixture_dist_density_kstest_2.R AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05_pruned.2.Q AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05_pruned.fam ASW_dip_anc1.txt

#3 Extract causal variants
# calculate frequency
plink --bfile AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05 --freq --out AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05
# sample causal variants
Rscript CausalVariantExtract.R AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05.frq
# extract their genotype
plink --bfile AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05 --extract CausalVariantExtract.txt --make-bed --out AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05_extract
#1000 variants and 10100 people pass filters and QC.

#4 Simulate effect sizes for these variants 
# with a normal distribution (μ=0, var=0.8).
Rscript CausalVariantScore.R AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05.frq

#5 Calculate genetic score
plink --bfile AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05_extract --score CausalVariantScore.txt 
#--score: Results written to plink.profile .

# get genetic score of admixed ppl only
tail -10000 plink.profile > sim_10Kppl_geneticscore.txt

# remove 100 non-admixed ppl
head -100 AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05_extract.fam | awk "{print $1,$2}" > remove_100ppl.txt   
plink --bfile AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05_extract --remove remove_100ppl.txt --make-bed --out AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05_extract_1kSNP
# 1000 variants and 10000 people pass filters and QC.

#6,7,8 Simulate phenotype data
Rscript sim_pheno.R
# get pheno data sim_pheno_10Kppl_random_2.txt

#9 Estimate heritability
# construct grm
gcta64 --bfile AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05_extract_1kSNP --autosome --make-grm-bin --out AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05_extract_1kSNP
# estimate heritability
gcta64 --reml --grm-bin AA_10Mb_random_10000ppl_chr1_10_AT_maf0.05_extract_1kSNP --pheno sim_pheno_10Kppl_random_2.txt --out GCTAherit_sim_pheno_10Kppl_random_2_1kSNP

