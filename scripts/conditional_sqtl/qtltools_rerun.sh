#!/bin/bash

mkdir -p ../output/01.qtltools_re/cond_norm_pbs/
mkdir -p ../output/01.qtltools_re/cond_perm_pbs/
sbatch /work/users/s/e/seyoun/CQTL_sQTL/scripts/conditional_sqtl/qtltools_rerun.sbatch $1 $2 $3 $4 $5 $6 $7
