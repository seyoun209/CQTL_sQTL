for i in /work/users/s/e/seyoun/CQTL_sQTL/output/align/*sorted.bam
do
outdir="/work/users/s/e/seyoun/CQTL_sQTL/output/snrnp70_subset_bam"
mkdir -p ${outdir}
nm=${i##*/}
sample_nm=${nm%%.sorted.bam}
echo $sample_nm
bam_dir=${outdir}/${sample_nm}_snrnp70.bam
sam_dir=${outdir}/${sample_nm}_snrnp70.sam
bai_dir=${outdir}/${sample_nm}_snrnp70.bai
echo $bam_dir
sbatch snrnp70_subset.sbatch $i ${bam_dir} ${sam_dir} ${bai_dir}
done

