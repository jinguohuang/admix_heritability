#!/bin/bash

module load plink/2.00-alpha-091019
module load anaconda/miniconda3_4.8.3-jupyter


p=${1}
rep=${2}

plink2 --pfile ~/lab/aazaidi/projects/admix_heritability/data/ldsc/plink/admix_chrall \
--extract ~/lab/aazaidi/projects/admix_heritability/data/ldsc/plink/admix_chrall.cm.snplist \
--allow-extra-chr \
--seed ${rep} \
--thin ${p} --export A --out ~/lab/aazaidi/projects/admix_heritability/data/ldsc/plink/admix_chrall.cm_thin${p}_${rep}

# plink2 --pfile ~/lab/aazaidi/projects/admix_heritability/data/ldsc/plink/admix_chrall.cm_thin${p}_${rep} \
# --allow-extra-chr \
# --freq --out ~/lab/aazaidi/projects/admix_heritability/data/ldsc/plink/admix_chrall.cm_thin${p}_${rep}_freq


Rscript generate_effects.R ${p} ${rep}

# plink2 --pfile ~/lab/aazaidi/projects/admix_heritability/data/ldsc/plink/admix_chrall \
# --allow-extra-chr \
# --score ~/lab/aazaidi/projects/admix_heritability/data/ldsc/effects/admix_chrall.cm_thin${p}_${rep}.effects \
# --out ~/lab/aazaidi/projects/admix_heritability/data/ldsc/pgs/admix_chrall.cm_thin${p}_${rep}.pgs


#run gwas
plink2 --pfile ~/lab/aazaidi/projects/admix_heritability/data/ldsc/plink/admix_chrall \
--glm hide-covar \
--pheno ~/lab/aazaidi/projects/admix_heritability/data/ldsc/phenotypes/admix_chrall_p${p}_${rep}.pgs \
--pheno-name y \
--covar ~/lab/aazaidi/projects/admix_heritability/data/ldsc/plink/admix_chrall.cm.50k.pca.eigenvec \
--covar-name PC1 \
--allow-extra-chr \
--out ~/lab/aazaidi/projects/admix_heritability/data/ldsc/gwas/admix_chrall_p${p}_${rep}

#compute ldscores







