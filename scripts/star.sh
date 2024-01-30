
ml star/2.7.11a
mkdir -p output/align

STAR --genomeDir $6 \
	--runThreadN $7 \
	--sjdbFileChrStartEnd $8 \
	--outFileNamePrefix $9 \
	--quantMode TranscriptomeSAM \
	--outSAMstrandField intronMotif \
	--readFilesCommand zcat \
	--outSAMtype BAM SortedByCoordinate \
	--readFilesIn $1 $2 \
	--outFilterType BySJout \
	--outFilterMultimapNmax 20 \
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 \
	--outFilterMismatchNmax 999 \
	--outFilterMismatchNoverReadLmax 0.04 \
	--alignIntronMin 20 \
	--alignIntronMax 1000000 \
	--alignMatesGapMax 1000000 \
	--waspOutputMode SAMtag \
	--varVCFfile <(zcat $3) 1> $4 2> $5

