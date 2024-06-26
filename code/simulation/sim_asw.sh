#!/bin/bash

rep=${1}

## code block already done

# plink2 --memory 8000 --threads 2 \\
# --pfile ~/lab/riahi006/Dependencies/Data/1kG.plink/all_hg38 \\
# --keep 1kg.asw.samples \\
# --make-pgen \\
# --out test \\
# --maf 0.01 --chr 1-22 --allow-extra-chr

# #remove duplicates
# plink2 --pfile test \
# --threads 2 \
# --memory 8000 \
# --rm-dup exclude-all \
# --make-pgen \
# --out 1kg.asw.rmdup

# thin for LD
# plink2 --pfile 1kg.asw.rmdup \
# --threads 2 \
# --memory 8000 \
# --indep-pairwise 100 10 0.5 \
# --out 1kg.asw.rmdup.ldpr

# plink2 --pfile 1kg.asw.rmdup \
# --threads 2 \
# --memory 8000 \
# --extract 1kg.asw.rmdup.ldpr.prune.in \
# --make-pgen \
# --out 1kg.asw.rmdup.ldpr


#variant thinning
# plink2 --pfile 1kg.asw.rmdup \
# --threads 2 \
# --memory 8000 \
# --thin-count 200000 \
# --make-pgen \
# --out 1kg.asw.rmdup.thinned

#compute allele frequency
# plink2 --pfile 1kg.asw.rmdup.thinned \
# --freq --out 1kg.asw.rmdup.thinned.f

#select a set of list of 1000 randomly chosen snps and output genotypes

tail -n+2 ~/projects/admix_heritability/data/asw/1kg.asw.rmdup.thinned.f.afreq | \
shuf -n1000 > ~/projects/admix_heritability/data/asw/1kg.asw.rmdup.thinned.${rep}.snplist

plink2 --pfile ~/projects/admix_heritability/data/asw/1kg.asw.rmdup.thinned \
--extract ~/projects/admix_heritability/data/asw/1kg.asw.rmdup.thinned.${rep}.snplist \
--export A \
--out ~/projects/admix_heritability/data/asw/1kg.asw.rmdup.thinned.${rep}

echo "outputting variance components"

Rscript ~/projects/admix_heritability/code/simulation/sim_effects_asw.R ${rep}



