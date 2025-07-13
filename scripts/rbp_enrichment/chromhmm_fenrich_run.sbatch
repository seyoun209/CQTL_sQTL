#!/bin/bash
#SBATCH -p general
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --mem=10g
#SBATCH -t 1-00:00:00

# Load required modules
ml qtltools/1.3.1

# Set variables based on input arguments
CONDITION=$1
STATE=$2
BED_FILE=$3

BASE_DIR="/work/users/s/e/seyoun/CQTL_sQTL"
CHROMHMM_OUTPUT_DIR="${BASE_DIR}/output/Enrichment/chromhmm"
RBP_DATA_PREP_DIR="${BASE_DIR}/output/Enrichment/rbp/data_prep"

# Create output directory if it doesn't exist
mkdir -p "${CHROMHMM_OUTPUT_DIR}/${CONDITION}"

# Run QTLtools fenrich
QTLtools fenrich --qtl ${RBP_DATA_PREP_DIR}/significant_${CONDITION,,}_rank0.bed \
        --tss ${RBP_DATA_PREP_DIR}/significant_${CONDITION,,}_quantified_rank0.bed \
        --bed ${BED_FILE} \
        --out ${CHROMHMM_OUTPUT_DIR}/${CONDITION}/${STATE}_${CONDITION,,}_sQTL_enrichment.txt
