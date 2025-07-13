#!/bin/bash

BASE_DIR="/work/users/s/e/seyoun/CQTL_sQTL"
OUTPUT_DIR="${BASE_DIR}/output/Enrichment/chromhmm"
DATA_PREP_DIR="${OUTPUT_DIR}/data_prep"

# List of states
STATES=("1_TssA" "2_TssAFlnk" "3_TxFlnk" "4_Tx" "5_TxWk" "6_EnhG" "7_Enh" "8_ZNF_Rpts" "9_Het" "10_TssBiv" "11_BivFlnk" "12_EnhBiv" "13_ReprPC" "14_ReprPCWk" "15_Quies")

# Loop through states
for state in "${STATES[@]}"; do
    for condition in PBS FNF; do
        echo "Submitting job for ${condition} - ${state}"
        sbatch --job-name=fenrich_chromhmm_${condition}_${state} \
               --error=${BASE_DIR}/output/logs/fenrich_chromhmm_${condition}_${state}_e.%j \
               --output=${BASE_DIR}/output/logs/fenrich_chromhmm_${condition}_${state}_o.%j \
               chromhmm_fenrich_run.sbatch ${condition} ${state} "${DATA_PREP_DIR}/${state}.bed"
    done
done

echo "All jobs submitted"
