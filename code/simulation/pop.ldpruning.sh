#!/bin/bash

rep=${1}
pop=${2}

## code block already done

# plink2 --memory 8000 --threads 2 \
# --pfile ~/lab/riahi006/Dependencies/Data/1kG.plink/all_hg38 \
# --keep ~/projects/admix_heritability/data/1kg.${pop}.samples \
# --make-pgen \
# --out ~/projects/admix_heritability/data/${pop}/1kg.${pop} \
# --maf 0.01 --chr 1-22 --allow-extra-chr

#  remove duplicates
# plink2 --pfile ~/projects/admix_heritability/data/${pop}/1kg.${pop} \
# --threads 2 \
# --memory 8000 \
# --rm-dup exclude-all \
# --make-pgen \
# --out ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup

# thin for LD
plink2 --pfile ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup \
--indep-pairwise 100 10 0.5 \
--out ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.ldpr
--threads 2 \

plink2 --pfile ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup \
--extract ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.ldpr.prune.in \
--make-pgen 
--out ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.ldpr

# variant thinning
# plink2 --pfile ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup \
# --threads 2 \
# --memory 8000 \
# --thin-count 200000 \
# --make-pgen \
# --out ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.thinned \

# compute allele frequency
# plink2 --pfile ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.thinned \
# --freq --out ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.thinned.f \