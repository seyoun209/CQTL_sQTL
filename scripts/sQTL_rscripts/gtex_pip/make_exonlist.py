import pandas as pd
import qtl.annotation

annot = qtl.annotation.Annotation('/work/users/s/e/seyoun/Ref/genome/gencode.v45.annotation.genes.gtf')
exon_df = pd.DataFrame([[g.chr, e.start_pos, e.end_pos, g.strand, g.id, g.name]
                        for g in annot.genes for e in g.transcripts[0].exons],
                       columns=['chr', 'start', 'end', 'strand', 'gene_id', 'gene_name'])
exon_df.to_csv('/work/users/s/e/seyoun/Ref/genome/gencode.v45.annotation.genes.exons.txt.gz', sep='\t', index=False)
