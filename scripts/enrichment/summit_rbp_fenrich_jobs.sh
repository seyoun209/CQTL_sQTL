#!/bin/bash
BASE_DIR="/work/users/s/e/seyoun/CQTL_sQTL"
OUTPUT_DIR="${BASE_DIR}/output/Enrichment/eclip_rbp"
DATA_PREP_DIR="${BASE_DIR}/external_data/encode_rbp/rbp_prep_rbpOnly"

# Loop through RBP bed files
for bed_file in ${DATA_PREP_DIR}/*.bed; do
    rbp=$(basename ${bed_file} .bed)
    for condition in PBS FNF; do
        echo "Submitting job for ${condition} - ${rbp}"
        sbatch --job-name=fenrich_${condition}_${rbp} \
               --error=${BASE_DIR}/output/logs/fenrich_${condition}_${rbp}_e.%j \
               --output=${BASE_DIR}/output/logs/fenrich_${condition}_${rbp}_o.%j \
               fenrich_rbp_all_run.sbatch ${condition} ${rbp} ${bed_file}
    done
done

echo "All jobs submitted"
