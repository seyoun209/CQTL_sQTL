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
RBP=$2
BED_FILE=$3

BASE_DIR="/work/users/s/e/seyoun/CQTL_sQTL"
OUTPUT_DIR="${BASE_DIR}/output/Enrichment/eclip_rbp"
RBP_DATA_PREP_DIR="${BASE_DIR}/output/Enrichment/rbp/data_prep"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}/${CONDITION}"

# Echo commands for debugging
echo "Running QTLtools fenrich with the following parameters:"
echo "CONDITION: ${CONDITION}"
echo "RBP: ${RBP}"
echo "BED_FILE: ${BED_FILE}"
echo "OUTPUT_DIR: ${OUTPUT_DIR}"
echo "RBP_DATA_PREP_DIR: ${RBP_DATA_PREP_DIR}"

# Run QTLtools fenrich
QTLtools fenrich --qtl "${RBP_DATA_PREP_DIR}/significant_${CONDITION,,}_rank0.bed" \
        --tss "${RBP_DATA_PREP_DIR}/significant_${CONDITION,,}_quantified_rank0.bed" \
        --bed "${BED_FILE}" \
        --out "${OUTPUT_DIR}/${CONDITION}/${RBP}_${CONDITION,,}_sQTL_enrichment.txt"

# Check if the output file was created
if [ -f "${OUTPUT_DIR}/${CONDITION}/${RBP}_${CONDITION,,}_sQTL_enrichment.txt" ]; then
    echo "Output file successfully created."
else
    echo "Error: Output file was not created."
    echo "Command used:"
    echo "QTLtools fenrich --qtl ${RBP_DATA_PREP_DIR}/significant_${CONDITION,,}_rank0.bed --tss ${RBP_DATA_PREP_DIR}/significant_${CONDITION,,}_quantified_rank0.bed --bed ${BED_FILE} --out ${OUTPUT_DIR}/${CONDITION}/${RBP}_${CONDITION,,}_sQTL_enrichment.txt"
fi
