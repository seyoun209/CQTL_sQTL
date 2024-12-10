# CQTL_sQTL

This pipeline is used for processing and analyzing differential splicing, sQTLs, re-sQTLs, and colocalizations.
The associated paper is available on BioRxiv: doi: 10.1101/2024.11.11.622754

##Workflow

Clone workflow into working directory

```bash
git clone git@github.com:seyoun209/CQTL_sQTL.git .
```
### Differential alternative splicing

1. Prepare tap-seperated 'samplesheet.txt' file with the names of Rea1 and Read2 gzipped fastq files, and the paths to theses files under the Seqeuncing_Directory column. An example is shown below:

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
sbatch run_RNAprocessing
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

Colocalization were performed with [Boer et al. 2021](https://www.cell.com/cell/fulltext/S0092-8674(21)00941-7?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421009417%3Fshowall%3Dtrue) GWAS, which includes data for 11 defined phenotypes encompassing major sites of OA.  

 
