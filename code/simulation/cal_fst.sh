#/!bin/bash

grep -E '(YRI|CEU)' ~/projects/admix_heritability/data/1kg.samples > ~/projects/admix_heritability/data/1kg.ceu.yri.samples


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
