#!/bin/bash

module load R
rep=${1}
pop=${2}
n=${3}


#select a set of list of 1000 randomly chosen snps and output genotypes

mkdir -p ~/projects/admix_heritability/data/${pop}/snplist
mkdir -p ~/projects/admix_heritability/data/${pop}/vg
mkdir -p ~/projects/admix_heritability/data/${pop}/raw

tail -n+2 ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.ldpr.f.afreq | \
shuf -n${n} > ~/projects/admix_heritability/data/${pop}/snplist/1kg.${pop}.rmdup.ldpr.${n}.${rep}.snplist

plink2 --pfile ~/projects/admix_heritability/data/${pop}/1kg.${pop}.rmdup.ldpr \
--extract ~/projects/admix_heritability/data/${pop}/snplist/1kg.${pop}.rmdup.ldpr.${n}.${rep}.snplist \
--export A \
--out ~/projects/admix_heritability/data/${pop}/raw/1kg.${pop}.rmdup.ldpr.${n}.${rep}

echo "outputting variance components"

Rscript ~/projects/admix_heritability/code/simulation/sim_effects.R \
~/projects/admix_heritability/data/${pop}/raw/1kg.${pop}.rmdup.ldpr.${n}.${rep}.raw \
${rep} \
${pop} \
${n} \
~/projects/admix_heritability/data/${pop}/vg/1kg.${pop}.rmdup.ldpr.${n}.${rep}.vg







