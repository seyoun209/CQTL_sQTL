#!/bin/bash

max_jobs=500
max_running_jobs=100


ml plink/1.90b3
ml qtltools/1.3.1
ml samtools/1.20

#mkdir -p /work/users/s/e/seyoun/CQTL_sQTL/output/geno/vcf_bed/

#vcf_in="/work/users/s/e/seyoun/CQTL_sQTL/output/geno/pbs_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz"
#output="/work/users/s/e/seyoun/CQTL_sQTL/output/geno/vcf_bed/"
#for i in {1..22}; do
#  echo "Processing chromosome $1..."
#  echo ${vcf_in}
#  bcftools view ${vcf_in} --regions chr${i} -o ${output}/chr${i}.no_indel.vcf.gz -Oz
#  tabix -p vcf ${output}/chr${i}.no_indel.vcf.gz
#  plink --vcf ${output}/chr${i}.no_indel.vcf.gz --make-bed --double-id --out ${output}/chr${i}
#
#done

#mkdir -p /work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld
#mkdir -p /work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld05
#chr_nm="$(awk '{split($1,a,":"); print a[1]}' /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/sigSNPs_all.list)"

snps="$(awk '{print $1}' /work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/sigSNPs_all_includeCond)"


while IFS= read -r i || [[ -n "$i" ]]; do
    chrnm=$(cut -d':' -f1 <<< ${i})
    ref_EUR="/proj/phanstiel_lab/References/genomes/1000G/GRCh38/EUR/biallelic_snps/EUR_1000G.GRCh38.20190312_${chrnm}_biallelic"
    ref_ALL="/proj/phanstiel_lab/References/genomes/1000G/GRCh38/ALL/biallelic_snps/1000G.GRCh38.20181129_${chrnm}_biallelic"
    echo ${ref_EUR}
    echo ${ref_ALL}
    sbatch ld_calc.sbatch $chrnm $i ${ref_EUR} ${ref_ALL}
    
    if [ $(squeue -u $USER -h | wc -l) -ge $max_jobs ]; then
        echo "Reached $max_jobs jobs. Waiting for jobs to complete..."
        
        while [ $(squeue -u $USER -h | wc -l) -gt $max_running_jobs ]; do
            sleep 10  # Check every minute
        done
        
        echo "Batch completed. Continuing with next batch."
    fi
done < <(echo "$snps")

# Wait for any remaining jobs
if [ $(squeue -u $USER -h | wc -l) -gt $max_running_jobs ]; then
    echo "Waiting for final jobs to complete..."
    while [ $(squeue -u $USER -h | wc -l) -gt $max_running_jobs ]; do
        sleep 10
    done
fi



#for i in $snps
#do
#	chrnm=$(cut -d':' -f1 <<< ${i})
#	#chr=${chrnm##chr}
#	echo $chrnm
#	echo $i
#	ref_EUR="/proj/phanstiel_lab/References/genomes/1000G/GRCh38/EUR/biallelic_snps/EUR_1000G.GRCh38.20190312_${chrnm}_biallelic"
#	ref_ALL="/proj/phanstiel_lab/References/genomes/1000G/GRCh38/ALL/biallelic_snps/1000G.GRCh38.20181129_${chrnm}_biallelic"
#
#	echo ${ref_EUR}
#	echo ${ref_ALL}
#	sbatch ld_calc.sbatch $chrnm $i ${ref_EUR} ${ref_ALL}
#done
