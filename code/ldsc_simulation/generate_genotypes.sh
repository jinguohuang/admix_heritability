#!/bin/bash

module load conda
source activate msprime-env
module load plink/2.00-alpha-091019

chrom=${1}

python generate_vcf.py -c ${chrom} \
-o ~/lab/aazaidi/projects/admix_heritability/data/ldsc/vcf/admix

gzip -f ~/lab/aazaidi/projects/admix_heritability/data/ldsc/vcf/admix_chr${chrom}.vcf


