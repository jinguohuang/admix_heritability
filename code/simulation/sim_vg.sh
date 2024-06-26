#!/bin/bash

module load R
rep=${1}
pop=${2}
n=${3}


#select a set of list of 1000 randomly chosen snps and output genotypes

tail -n+2 ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.thinned.f.afreq | \
shuf -n1000 > ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.thinned.${n}.${rep}.snplist

plink2 --pfile ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.thinned \
--extract ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.thinned.${n}.${rep}.snplist \
--export A \
--out ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.thinned.${n}.${rep}

echo "outputting variance components"

Rscript ~/projects/admix_heritability/code/simulation/sim_effects.R \
~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.thinned.${n}.${rep}.raw \
${rep} \
${pop} \
~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.thinned.${n}.${rep}.vg







