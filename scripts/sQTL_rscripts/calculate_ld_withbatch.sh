#!/bin/bash
max_jobs=500
max_running_jobs=100
ml plink/1.90b3
ml qtltools/1.3.1
ml samtools/1.20

# Read snps into a file
awk '{print $1}' /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/sigSNPs_all.list > temp_snps.txt

while IFS= read -r i; do
    while [ $(squeue -u $USER -h | wc -l) -ge $max_jobs ]; do
        echo "Reached $max_jobs jobs. Waiting for jobs to complete..."
        sleep 60  # Wait for 1 minute before checking again
    done

    chrnm=$(cut -d':' -f1 <<< ${i})
    ref_EUR="/proj/phanstiel_lab/References/genomes/1000G/GRCh38/EUR/biallelic_snps/EUR_1000G.GRCh38.20190312_${chrnm}_biallelic"
    ref_ALL="/proj/phanstiel_lab/References/genomes/1000G/GRCh38/ALL/biallelic_snps/1000G.GRCh38.20181129_${chrnm}_biallelic"
    echo ${ref_EUR}
    echo ${ref_ALL}
    sbatch ld_calc.sbatch $chrnm $i ${ref_EUR} ${ref_ALL}
    
    current_jobs=$(squeue -u $USER -h | wc -l)
    echo "Current number of jobs: $current_jobs"
done < temp_snps.txt

# Wait for any remaining jobs
while [ $(squeue -u $USER -h | wc -l) -gt $max_running_jobs ]; do
    echo "Waiting for final jobs to complete..."
    sleep 60
done

echo "All jobs completed."

# Clean up
rm temp_snps.txt
