#!/bin/bash

# simulate admix lanc
Rscript AssortMate.R --seed ${1} -P ${2} -t ${3}

filename=admix_fst0.2_psel0_seed${1}_P${2}_t${3}

# estimate vg due to lanc 
./vg_lanc.sh ${filename}

# estimate vg due to lanc in meta pop 
./vg_lanc.sh meta_${filename}

# reformat the output
Rscript reformat_expobs.R ${filename}

# remove tmp files
rm ${filename}.* 
rm meta_${filename}.*
rm sim_${filename}* 
