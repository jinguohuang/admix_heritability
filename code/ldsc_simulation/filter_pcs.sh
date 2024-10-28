#!/bin/bash

module load plink/2.00-alpha-091019


plink2 \
  --maf 0.01 \
  --out ~/lab/aazaidi/projects/admix_heritability/data/ldsc/plink/admix_chrall.cm \
  --pfile ~/lab/aazaidi/projects/admix_heritability/data/ldsc/plink/admix_chrall \
  --allow-extra-chr \
  --write-snplist

plink2 \
  --extract ~/lab/aazaidi/projects/admix_heritability/data/ldsc/plink/admix_chrall.cm.snplist \
  --out ~/lab/aazaidi/projects/admix_heritability/data/ldsc/plink/admix_chrall.cm.50k \
  --pfile ~/lab/aazaidi/projects/admix_heritability/data/ldsc/plink/admix_chrall \
  --thin-count 50000 \
  --allow-extra-chr \
  --write-snplist

plink2 \
  --pfile ~/lab/aazaidi/projects/admix_heritability/data/ldsc/plink/admix_chrall \
  --extract ~/lab/aazaidi/projects/admix_heritability/data/ldsc/plink/admix_chrall.cm.50k.snplist \
  --allow-extra-chr \
  --pca 20 \
  --out ~/lab/aazaidi/projects/admix_heritability/data/ldsc/plink/admix_chrall.cm.50k.pca