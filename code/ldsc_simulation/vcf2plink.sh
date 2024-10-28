

#!/bin/bash

module load plink/2.00-alpha-091019
module load bcftools

echo "merging vcf files"

bcftools concat ~/lab/aazaidi/projects/admix_heritability/data/ldsc/vcf/admix_chr[1-100].vcf.gz \
-o ~/lab/aazaidi/projects/admix_heritability/data/ldsc/vcf/admix_chrall.vcf.gz \
-O z

bcftools annotate --rename-chrs ~/lab/aazaidi/projects/admix_heritability/data/ldsc/vcf/chrmap.txt \
-Oz -o ~/lab/aazaidi/projects/admix_heritability/data/ldsc/vcf/admix_chrall_recontig.vcf.gz \
~/lab/aazaidi/projects/admix_heritability/data/ldsc/vcf/admix_chrall.vcf.gz 

bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
-Oz -o ~/lab/aazaidi/projects/admix_heritability/data/ldsc/vcf/admix_chrall_recontig_reid.vcf.gz \
~/lab/aazaidi/projects/admix_heritability/data/ldsc/vcf/admix_chrall_recontig.vcf.gz

bcftools view --max-alleles 2 -Oz -o ~/lab/aazaidi/projects/admix_heritability/data/ldsc/vcf/admix_chrall_recontig_reid_snps.vcf.gz \
~/lab/aazaidi/projects/admix_heritability/data/ldsc/vcf/admix_chrall_recontig_reid.vcf.gz


plink2 --vcf ~/lab/aazaidi/projects/admix_heritability/data/ldsc/vcf/admix_chrall_recontig_reid_snps.vcf.gz \
--double-id \
--make-pgen \
-allow-extra-chr \
--out ~/lab/aazaidi/projects/admix_heritability/data/ldsc/plink/admix_chrall


