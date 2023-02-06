#!/bin/bash

# simulate admix lanc
Rscript NeutralEvo.R --seed ${1}

# estimate vg due to lanc
./vg_lanc.sh admix0.5_fst0.2_psel0_seed${1}

# reformat the output
Rscript reformat_expobs.R admix0.5_fst0.2_psel0_seed${1}

