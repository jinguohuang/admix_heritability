# Phenotype simulation
### Basic model
code: **simulation_phenotype.R**  
y=βX+ϵ  
freq\~unif(0.1,0.9)  
X\~B(1000,p=freq)  
β\~N(0,var=0.8)  
ϵ\~N(0,var=0.2)

### Apply the basic model to simulated genotype data to get phenotype data
code: **sim_pheno_basic.sh**
1) Simulating genetic data of 10,000ppl under random mating model: 10M for each chromosome, simulate 10 chromosomes. 
2) QC the genotype with PLINK: --mind 0.1 --geno 0.1 --maf 0.05
3) Sample 1,000 variants at random such that they are uniformly spaced across the genome with ~ 100Kb between them. We do this for each chromosome separately.  
4) Simulate effect sizes for these variants with a normal distribution (μ=0, var=0.8).  
5) Use these effects and the genotypes (not pruned) to calculate the simulated genetic value for each individual in PLINK.
6) Simulate environmental value for these variants with a normal distribution (μ=0, var=0.2).  
7) Normalize the genetic value and add environmental value together to get phenotype data.
8) Normalize phenotype data. 
9) Estimate heritability with gcta.


### Apply (Schoech et al. 2019; Zaidi et al. Elife 2020) model to effect sizes 
code: **sim_pheno_V1.sh**
1) Simulating genetic data of 10,000ppl under random mating model: 10M for each chromosome, simulate 10 chromosomes. 
2) QC the genotype with PLINK: --mind 0.1 --geno 0.1 --maf 0.05
3) Sample 1,000 variants at random such that they are uniformly spaced across the genome with ~ 100Kb between them. We do this for each chromosome separately.  
4) Apply the model to effect sizes the effects we assign to causal variants will depend on their frequency with rare variants having larger effects and common variants having smaller effects.
5) Use these effects and the genotypes (not pruned) to calculate the simulated genetic value for each individual in PLINK.
6) Simulate environmental value for these variants with a normal distribution (μ=0, var=0.2).  
7) Normalize the genetic value and add environmental value together to get phenotype data.
8) Normalize phenotype data. 
9) Estimate heritability with gcta.
