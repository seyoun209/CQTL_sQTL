#!/bin/bash

BASE_DIR="/work/users/s/e/seyoun/CQTL_sQTL"
OUTPUT_DIR="${BASE_DIR}/output/Enrichment/eclip_tissue_sep"
DATA_PREP_DIR="${BASE_DIR}/external_data/encode_rbp/rbp_prep"

# List of tissues
TISSUES=("K562" "HepG2" "SM-9MVZL")

# Loop through tissues and RBPs
for tissue in "${TISSUES[@]}"; do
    for bed_file in ${DATA_PREP_DIR}/${tissue}_*.bed; do
        rbp=$(basename ${bed_file} .bed | cut -d'_' -f2-)
        for condition in PBS FNF; do
            echo "Submitting job for ${condition} - ${tissue} - ${rbp}"
            sbatch --job-name=fenrich_${condition}_${tissue}_${rbp} \
                   --error=${BASE_DIR}/output/logs/fenrich_${condition}_${tissue}_${rbp}_e.%j \
                   --output=${BASE_DIR}/output/logs/fenrich_${condition}_${tissue}_${rbp}_o.%j \
                   fenrich_run.sbatch ${condition} ${tissue} ${rbp} ${bed_file}
        done
    done
done

echo "All jobs submitted"
