# CQTL_sQTL

This pipeline is used for processing and analyzing differential splicing, sQTLs, re-sQTLs, and colocalizations.

The associated paper is available on [BioRxiv: doi: 10.1101/2024.11.11.622754](https://www.biorxiv.org/content/10.1101/2024.11.11.622754v1)
#### Abstract
Osteoarthritis affects millions worldwide, yet effective treatments remain elusive due to poorly understood molecular mechanisms. While genome-wide association studies (GWAS) have identified over 100 OA-associated loci, identifying the genes impacted at each locus remains challenging. Several studies have mapped expression quantitative trait loci (eQTL) in chondrocytes and colocalized them with OA GWAS variants to identify putative OA risk genes; however, the degree to which genetic variants influence OA risk via alternative splicing has not been explored. We investigated the role of alternative splicing in OA pathogenesis using RNA-seq data from 101 human chondrocyte samples treated with PBS (control) or fibronectin fragment (FN-f), an OA trigger. We identified 590 differentially spliced genes between conditions, with FN-f inducing splicing events similar to those in primary OA tissue. We used CRISPR/Cas9 to mimic an SNRNP70 splicing event observed in OA and FN-f-treated chondrocytes and found that it induced an OA-like expression pattern. Integration with genotyping data revealed 7,188 splicing quantitative trait loci (sQTL) affecting 3,056 genes. While many sQTLs were shared, we identified 738 and 343 condition-specific sQTLs for control and FN-f, respectively. We identified 15 RNA binding proteins whose binding sites were enriched at sQTL splice junctions and found that expression of those RNA binding proteins correlated with exon inclusion. Colocalization with OA GWAS identified 6 putative risk genes, including a novel candidate, PBRM1. Our study highlights the significant impact of alternative splicing in OA and provides potential therapeutic targets for future research.

## Workflow

Clone workflow into the working directory

```bash
git clone git@github.com:seyoun209/CQTL_sQTL.git .
```
### Differential alternative splicing

1. Prepare a tap-separated 'samplesheet.txt' file with the names of Rea1 and Read2 gzipped fastq files and the paths to these files under the Seqeuncing_Directory column. An example is shown below:

	| Project   | Cell_Type | Genotype	| Bio_Rep	| Tech_Rep	| Seq_Rep	| Read1 | Read2 | Sequencing_Directory |
	|---------|-----------|----------|---------|----------|---------|-------------------|-------------------|---------------------------| 
	| PROJ  | CELL  | WT	| 1 | 1 | 1 | sample1_R1.fq.gz  | sample1_R2.fq.gz	| /path/to/fastq/directory/ |
	| PROJ  | CELL  | WT	| 2 | 1 | 1 | sample2_R1.fq.gz  | sample2_R2.fq.gz	| /path/to/fastq/directory/ |
	| PROJ  | CELL  | MUT	| 1 | 1 | 1 | sample3_R1.fq.gz  | sample3_R2.fq.gz	| /path/to/fastq/directory/ |
	| PROJ  | CELL  | MUT	| 2 | 1 | 1 | sample4_R1.fq.gz  | sample4_R2.fq.gz	| /path/to/fastq/directory/ |

2. Edit the appropriate 'rna_prcoess.yaml'

3. Submit to SLURM with sbatch scripts:

```bash
sbatch run_RNAprocessing
sbatch run_diffsplicing
```

### sQTLs

4. Map sQTLs with:

```bash
sbatch run_sqtl
```

### response sQTLs

```bash
sbatch scripts/sQTL_rscripts/response.sbatch
```

### Colocalization
```R
Rscript scripts/coloc/coloc_calculate_all_condition_h4.R
```

Colocalization was performed with [Boer et al. 2021](https://www.cell.com/cell/fulltext/S0092-8674(21)00941-7?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421009417%3Fshowall%3Dtrue) GWAS, which includes data for 11 defined phenotypes encompassing major sites of OA.  

 
