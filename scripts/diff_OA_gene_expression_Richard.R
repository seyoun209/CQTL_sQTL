# Getting read for the Richard for his reference
setwd("/work/users/s/e/seyoun/CQTL_sQTL")
load("output/quant/deGene_OA_expression_dds.rda")
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(openxlsx)

#-------------------------------------------------------------------------------
#Differential gene expression
de_genes_shrink <- lfcShrink(dds_gene,
                             coef = "Condition_OA_vs_CTL", format = "GRanges") |>
  plyranges::names_to_column("gene_id")

de_genes_df <- as.data.frame(de_genes_shrink)
normalized_counts <- counts(dds_gene, normalized=TRUE)

# Get the condition information
condition <- colData(dds_gene)$Condition

# Calculate the baseMean for each condition
baseMean_CTL <- rowMeans(normalized_counts[, condition == "CTL"])
baseMean_OA <- rowMeans(normalized_counts[, condition == "OA"])


# Add gene symbols to DE results
de_genes_df <- de_genes_df |>
  dplyr::mutate(ensgID = gsub("\\..*$", "" , gene_id))

annots <- AnnotationDbi::select(org.Hs.eg.db,
                                de_genes_df$ensgID,
                                columns = c("SYMBOL"),
                                keytype = "ENSEMBL")
annots <- dplyr::rename(annots, gene_id = ENSEMBL)
annots_unique <- annots %>% distinct(gene_id, .keep_all = TRUE)
de_results_with_symbol <- left_join(de_genes_df, annots_unique, by = c("ensgID" = "gene_id"))

library(dplyr)

# Convert baseMean_CTL and baseMean_OA to data frames
baseMean_CTL_df <- data.frame(gene_id = names(baseMean_CTL), baseMean_PBS = as.numeric(baseMean_CTL), row.names = NULL)
baseMean_OA_df <- data.frame(gene_id = names(baseMean_OA), baseMean_OA = as.numeric(baseMean_OA), row.names = NULL)

# Merge with de_results_with_symbol
de_results_with_symbol <- de_results_with_symbol %>%
  left_join(baseMean_CTL_df, by = "gene_id") %>%
  left_join(baseMean_OA_df, by = "gene_id")

# Select specific columns
selected_df <- de_results_with_symbol %>%
  dplyr::select(ensgID, SYMBOL, log2FoldChange, pvalue, padj, baseMean_PBS, baseMean_OA)

# Sort by padj and log2FoldChange
sorted_df <- selected_df %>%
  arrange(padj)
write.csv(sorted_df, file = "OAvsPBS_gene_set.csv", col.names = TRUE, row.names = FALSE, quote = FALSE)


# View the result
head(sorted_df)
de_results_with_symbol |> filter(SYMBOL %in% c('KMT2D','KMT2C'))


