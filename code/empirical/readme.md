## Decomposing ${V}_g$ in traits in African Americans
* Script: **LDclump_hg38.sh, Allelef1f2_hg38.sh, vgp_calculator.R, gvalue_hg38.sh**
* Requirement: R version 4.2.3, PLINK 1.9, PLINK 2.0
* Usage: ```./LDclump_hg38.sh ${trait}```
* Example: To process GWAS Catalog downloaded summary statistics for trait HDL: LD clumping, calculate allele frequency of associated SNPs in CEU and YRI, calculate the proportion of variance explained, calculate the variance of polygenic scores: ```LDclump_hg38.sh HDL```
* Output: the proportion of variance explained for the trait ```(${trait}_vgp.txt)```, the variance and 95%CI of polygenic scores for the trait ```(${trait}_vgvalueCI.txt)```.
