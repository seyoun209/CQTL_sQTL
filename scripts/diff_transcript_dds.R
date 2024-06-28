setwd("/work/users/s/e/seyoun/CQTL_sQTL/")

## Load required libraries
library(data.table)
library(tximeta)
library(DESeq2)
library(dplyr)
library(plyranges)
library(liftOver)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)



aligned_samplesheet_path <- "aligned_samplesheet.txt"
donor_samples_path <- "./donor_samples.txt"
rna_extraction_path <- "./rna_extraction.txt"
conditions_to_include <- unlist(strsplit("CTL FNF", " ")) # Split the first argument into condition values
#use_wasp_id <- tail(args, 1) == "wasp" # Check if the last argument is "wasp"


# Load configuration from YAML file
config <- yaml::read_yaml("config/rna_prcoess.yaml")
aligned_samples <- fread(aligned_samplesheet_path)
donor_samples <- fread(donor_samples_path)
rna_extraction <- fread(rna_extraction_path)

# Identify common columns between aligned_samples and donor_samples, excluding 'Donor'
common_cols_donor <- setdiff(intersect(colnames(aligned_samples), colnames(donor_samples)), "Donor")

# Remove these common columns from donor_samples before the join
donor_samples_clean <- dplyr::select(donor_samples, -all_of(common_cols_donor))

# Perform the first left join
combined_data <- aligned_samples %>%
  left_join(donor_samples_clean, by = "Donor")

common_cols_rna <- setdiff(intersect(colnames(combined_data), colnames(rna_extraction)), "Read2")
rna_extraction_clean <- dplyr::select(rna_extraction, -all_of(common_cols_rna))
combined_data <- combined_data %>%
  left_join(rna_extraction_clean, by = "Read2")

# Process ID column with and without "_wasp"
combined_data <- combined_data %>%
  mutate(ID = paste(Donor, Condition, Tech_Rep, Sex, sep = "_"),
         ID_wasp = paste(Donor, Condition, Tech_Rep, Sex, "wasp", sep = "_"))

# Conditionally adjust FragmentBatch for 'OA' condition
combined_data <- combined_data %>%
  mutate(FragmentBatch = ifelse(Condition == "OA", 0, FragmentBatch))

# Omit samples and filter rows based on configuration and conditions
combined_data <- combined_data[!combined_data$Donor %in% config$samples_to_omit, ] %>%
  dplyr::filter(Condition %in% conditions_to_include)
selected_columns <- c("ID","Donor","Condition","Sex", "Age","Race","OAGradeAvg","CauseOfDeath","FragmentBatch","RIN","RNAextractionKitBatch","RNAshippedDate")
coldata <- combined_data %>% dplyr::select(all_of(selected_columns))

## Add quant paths and names
coldata$files <- file.path("output", "quant", coldata$ID, "quant.sf")
colnames(coldata) <- gsub("ID", "names", colnames(coldata))
file.exists(coldata$files)

## Import data with tximeta & summarize to gene
se <- tximeta(coldata)
gse <- summarizeToGene(se)
#-------------------------------------------------------------------------------
#  Transcript LEVEL
#-------------------------------------------------------------------------------
## Convert to factors (avoids a warning)
colData(se)[] <- lapply(colData(se), factor)

## Build DESeq object
dds <- DESeqDataSet(se, design = ~ Condition + Donor)

## Filter out lowly expressed genes
## (at least 10 counts in at least 2 samples)
keep <- rowSums(counts(dds) >= 10) >= ceiling(nrow(colData(se))*0.10)
dds <- dds[keep,]

## Fit model
dds <- DESeq(dds)

## Save dds
save(dds, file = "output/quant/differential_transcript_expression_dds.rda")

load("output/quant/differential_gene_expression_dds.rda")

vsd <- vst(dds, blind=FALSE)
ntd <- normTransform(dds)
vsd_update <- rownames_to_column(as.data.frame(assay(vsd)),"ENST")
ntd_update <- rownames_to_column(as.data.frame(assay(ntd)),"ENST")

fwrite(vsd_update,file="output/quant/normalized_vst_transcript.txt",sep="\t",quote=F,row.names=T,col.names=T)
fwrite(ntd_update,file="output/quant/normalized_log2_transcript.txt",sep="\t",quote=F,row.names=T,col.names=T)


#Differential gene expression
de_genes_shrink <- lfcShrink(dds_gene,
                             coef = "Condition_FNF_vs_CTL", format = "GRanges") |>
  plyranges::names_to_column("gene_id")

## outputs
# res <- results(dds_gene, name ="Condition_WD_vs_KD" )
# summary(res)
res <- lfcShrink(dds, coef="Condition_FNF_vs_CTL")
summary(res)
DESeq2::plotMA(res,ylim=c(-2,2),main = "Differential RNAseq CTL vs FNF Analysis",
               ylab = "LFC",
               xlab = "mean of norm. counts")

## plot PCA
#plotPCA(normTransform(dds_gene), intgroup = "Condition",ntop=23036 ) + ggplot2::theme(aspect.ratio = 1)
#plotPCA(normTransform(dds_gene), intgroup = "Condition",ntop=500) + ggplot2::theme(aspect.ratio = 1)

# Join results with gene info
de_genes_shrink <-
  inner_join(x = as.data.frame(de_genes_shrink),
             y = as.data.frame(rowData(gse)) %>%
               dplyr::select(c("gene_id", "tx_ids")),
             by = "gene_id") %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  keepStandardChromosomes(pruning.mode = "coarse") %>%
  as.data.frame()

write_csv(de_genes_shrink, file = "output/quant/de_transcript_results.csv")
#de_genes_shrink <- read.table("/work/users/s/e/seyoun/CQTL_sQTL/output/quant/de_genes_results.csv",header=T,sep=",") |> as.data.frame()
# Get significant genes
#sig_deGenes_pval01_l2fc15 <- de_genes_shrink |>
#  dplyr::filter(padj < 0.01 & abs(log2FoldChange) > 1.5) |>
#  write_csv("/work/users/s/e/seyoun/CQTL_sQTL/output/quant/sig_deGenes_pval01_l2fc15.csv")
#sig_deGenes_pval01_l2fc1.2 <- de_genes_shrink |>
#  dplyr::filter(padj < 0.01 & abs(log2FoldChange) > 1.2) |>
#  write_csv("/work/users/s/e/seyoun/CQTL_sQTL/output/quant/sig_deGenes_pval01_l2fc1.2.csv")




