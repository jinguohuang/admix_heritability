#!/bin/bash

module load R
rep=${1}
pop=${2}

## code block already done

#plink2 --memory 8000 --threads 2 \
# --pfile ~/lab/riahi006/Dependencies/Data/1kG.plink/all_hg38 \
# --keep ~/projects/admix_heritability/data/1kg.${pop}.samples \
# --make-pgen \
# --out ~/projects/admix_heritability/data/${pop}/1kg.${pop} \
# --maf 0.01 --chr 1-22 --allow-extra-chr

 #remove duplicates
# plink2 --pfile ~/projects/admix_heritability/data/${pop}/1kg.${pop} \
# --threads 2 \
# --memory 8000 \
# --rm-dup exclude-all \
# --make-pgen \
# --out ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup

# # thin for LD
# # plink2 --pfile ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup \
# # --threads 2 \
# # --memory 8000 \
# # --indep-pairwise 100 10 0.5 \
# # --out ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.ldpr
# # --threads 2 \

# # plink2 --pfile ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup \
# # --threads 2 \
# # --threads 2 \
# # --memory 8000 \
# # --extract ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.ldpr.prune.in \
# # --threads 2 \
# # --make-pgen 
# # --out ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.ldpr


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

#select a set of list of 1000 randomly chosen snps and output genotypes

tail -n+2 ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.thinned.f.afreq | \
shuf -n1000 > ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.thinned.${rep}.snplist

plink2 --pfile ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.thinned \
--extract ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.thinned.${rep}.snplist \
--export A \
--out ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.thinned.${rep}

echo "outputting variance components"

Rscript ~/projects/admix_heritability/code/simulation/sim_effects.R ${rep} ${pop}



