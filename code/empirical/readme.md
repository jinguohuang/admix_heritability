## Decomposing ${V}_g$ in traits in African Americans
* Script: **LDclump_hg38.sh, Allelef1f2_hg38.sh, vgCI_estimate.R**
* Requirement: R version 4.2.3, PLINK 1.9, PLINK 2.0
* Usage: ```./LDclump_hg38.sh ${trait}``` and ```vgCI_estimate.R``` 
* Example: To process GWAS Catalog downloaded summary statistics for trait HDL: LD clumping, calculate allele frequency of associated SNPs in CEU and YRI, calculate the proportion of variance explained: ```LDclump_hg38.sh HDL```, ```vgCI_estimate.R``` (input is Table S1) 
* Output: the proportion of variance explained for the trait with 95% confidence intervals ```vgpCI.txt```.
