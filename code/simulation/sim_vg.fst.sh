#!/bin/bash

module load R
rep=${1}
pop=${2}

# plink2 --pfile ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup \
# --extract ~/projects/admix_heritability/data/ceu.yri/1kg.ceu.yri.rmdup.fst0.05.fst \
# --make-pgen \
# --out ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.fst0.05


# # variant thinning
# plink2 --pfile ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.fst0.05 \
# --threads 2 \
# --memory 8000 \
# --thin-count 200000 \
# --make-pgen \
# --out ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.fst0.05.thinned \

# # compute allele frequency
# plink2 --pfile ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.fst0.05.thinned \
# --freq --out ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.fst0.05.thinned.f \


tail -n+2 ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.fst0.05.thinned.f.afreq | \
shuf -n1000 > ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.fst0.05.thinned.${rep}.snplist

plink2 --pfile ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.fst0.05.thinned \
--extract ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.fst0.05.thinned.${rep}.snplist \
--export A \
--out ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.fst0.05.thinned.${rep}

echo "outputting variance components"

Rscript ~/projects/admix_heritability/code/simulation/sim_effects.R \
~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.fst0.05.thinned.${rep}.raw \
${rep} \
 ${pop} \
 ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.fst0.05.thinned.${rep}.vg

