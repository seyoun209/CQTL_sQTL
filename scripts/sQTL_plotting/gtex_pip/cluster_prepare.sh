ml python/3.9.6
ml r/4.3.1

python3 cluster_perpare_fastqtl.py \
	/work/users/s/e/seyoun/07.Doug_thesis/09.cqtl_hic/oa_tissue/RNAPIPE-Splice/output/junc_filtered/juncfiles_gz.txt \
	/work/users/s/e/seyoun/07.Doug_thesis/09.cqtl_hic/oa_tissue/RNAPIPE-Splice/output/03.qtltools_noindel/gencode.v33.chr_patch_hapl_scaff.genes.exons.txt.gz \
	/work/users/s/e/seyoun/07.Doug_thesis/09.cqtl_hic/oa_tissue/RNAPIPE-Splice/output/03.qtltools_noindel/gencode.v33.genes.gtf \
	ctl_fnf \
	/work/users/s/e/seyoun/07.Doug_thesis/09.cqtl_hic/oa_tissue/RNAPIPE-Splice/output/junc_filtered/sample_participant_lookup.txt \
	--min_clu_reads 100 \
	--min_clu_ratio 0.001 \
	--max_intron_len 500000 \
	--num_pcs 20 \
	--leafcutter_dir /work/users/s/e/seyoun/07.Doug_thesis/09.cqtl_hic/RNAPIPE-Splice/tool/leafcutter \
	-o /work/users/s/e/seyoun/07.Doug_thesis/09.cqtl_hic/oa_tissue/RNAPIPE-Splice/output/03.qtltools_noindel/03.pre_cluster
