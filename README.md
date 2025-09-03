# CQTL_sQTL
[![DOI](https://zenodo.org/badge/748754087.svg)](https://doi.org/10.5281/zenodo.15875183)

This pipeline is used for processing and analyzing differential splicing, sQTLs, re-sQTLs, and colocalizations.

The associated paper is available on [Nat Comm: doi:https://doi.org/10.1038/s41467-025-63299-0](https://doi.org/10.1038/s41467-025-63299-0)
#### Abstract
Osteoarthritis affects millions worldwide, yet effective treatments remain elusive due to poorly understood molecular mechanisms. While genome-wide association studies (GWAS) have identified hundreds of osteoarthritis-associated loci, identifying the genes impacted at each locus remains challenging. We investigate alternative splicing using RNA-sequencing data from 101 human chondrocyte samples treated with phosphate-buffered saline or fibronectin fragment, an osteoarthritis trigger. We identified 590 differentially spliced genes between conditions, with FN-f inducing splicing events similar to those in primary osteoarthritis tissue. CRISPR/Cas9 mimicking of an SNRNP70 splicing event observed in osteoarthritis induced an osteoarthritis-like expression pattern. Integration with genotyping data revealed 7188 splicing quantitative trait loci (sQTL) affecting 3056 genes, including 738 and 343 condition-specific sQTLs for resting and fibronectin fragment, respectively. Colocalization with osteoarthritis GWAS identified 6 putative risk genes. Our study highlights the significant impact of alternative splicing in osteoarthritis and provides potential therapeutic targets.

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

 
