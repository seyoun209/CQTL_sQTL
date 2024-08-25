#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <base_ld_dir> <output_dir> <ld_threshold>"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$2"

# Create a temporary sbatch script
temp_script=$(mktemp)
cat << EOF > "$temp_script"
#!/bin/bash
#SBATCH -J qtltools
#SBATCH -p general
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --mem=15g
#SBATCH -t 3-00:00:00
#SBATCH -e ../../output/logs/combined_ld_qtls_e.%j
#SBATCH -o ../../output/logs/combined_ld_qtls_o.%j

cd /work/users/s/e/seyoun/CQTL_sQTL
ml r/4.3.1

# Run the R script with the provided arguments
Rscript /work/users/s/e/seyoun/CQTL_sQTL/scripts/coloc/ld05_combined_sQTL.R "$1" "$2" "$3"
EOF

# Submit the job using the temporary script
sbatch "$temp_script"

# Remove the temporary script
rm "$temp_script"

echo "Job submitted for processing LD data with base directory $1, output directory $2, and LD threshold $3"
