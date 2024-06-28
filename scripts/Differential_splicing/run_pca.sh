#!/bin/bash
  
ml python/3.9.6

count_gz=$1
pca_numb=$2



sbatch -t 10:00:00 --mem=4G --wrap="bash -c 'python /work/users/s/e/seyoun/CQTL_sQTL/tools/leafcutter/scripts/prepare_phenotype_table.py ${count_gz} -p ${pca_numb}'"
