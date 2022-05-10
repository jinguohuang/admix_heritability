#!/bin/bash

# simulate admix lanc
Rscript AssortMate.R --seed ${1} -P ${2}

# estimate vg due to lanc
./vg_lanc.sh admix_fst0.2_psel0_seed${1}_P${2}

# reformat the output
Rscript reformat_expobs.R admix_fst0.2_psel0_seed${1}_P${2}

# remove tmp files
rm admix_fst0.2_psel0_seed${1}_P${2}.* 

