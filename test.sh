for i in /work/users/s/e/seyoun/CQTL_sQTL/output/trim/*fq.gz
do
    sampleName=$(basename "$i" | cut -d '_' -f 2)
    
    echo "$sampleName"
    
    if [[ $sampleName == *"CTL"* ]]; then
        vcf_path="output/geno/pbs_geno/pbs.rename_wCHR.vcf"
    elif [[ $sampleName == *"FNF"* ]]; then
        vcf_path="output/geno/fnf_geno/fnf.rename_wCHR.vcf"
    elif [[ $sampleName == *"OA"* ]]; then
        vcf_path="output/geno/oa_geno/oa.rename_wCHR.vcf"
    else
        echo "Error: Unrecognized sampleName $sampleName"
        exit 1
    fi
    
    echo "$vcf_path"
done

