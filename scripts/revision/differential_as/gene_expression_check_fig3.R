# Check the gene expression for the condition specific. 
setwd("/work/users/s/e/seyoun/CQTL_sQTL/output")

library(rtracklayer)
library(dplyr)
library(DESeq2)
library(tibble)

#------------------------------------------------------------
#load the data
response_pbs_results_nolmer <- fread("01.qtltools_re/response_qtl/reseponsesQTL_PBS_significant.csv")
response_fnf_results_nolmer <- fread("01.qtltools_re/response_qtl/reseponsesQTL_FNF_significant.csv")

#Gene expression data

load("quant/differential_gene_expression_dds.rda") # dds_gene

# Count data
counts_gene <- counts(dds_gene, normalized = TRUE)

counts_slc26a4 <- as.data.frame(counts_gene) |>
  rownames_to_column(var = "ENSG") |> 
    filter(ENSG == "ENSG00000091137.14") |>
  dplyr::select(-ENSG)

counts_mapk8 <- as.data.frame(counts_gene) |>
  rownames_to_column(var = "ENSG") |> 
    filter(ENSG == "ENSG00000107643.17") |>
  dplyr::select(-ENSG)

de_genes_shrink <- lfcShrink(dds_gene,
                             coef = "Condition_FNF_vs_CTL",
                             format = "GRanges") |>
  plyranges::names_to_column("gene_id")

# Filter for SLC26A4 (ENSG00000091137.14)
de_slc26a4 <- de_genes_shrink %>%
  filter(gene_id == "ENSG00000091137.14")

## MAPK8

de_mapk8 <- de_genes_shrink %>%
  filter(gene_id == "ENSG00000107643.17") # MAPK8


vsd_gene_update <- fread("quant/normalized_vst_gene.txt")

#ENSG00000091137.14 (SLC26A4)


vsd_slc26a4 <- vsd_gene_update %>%
  filter(ENSG == "ENSG00000091137.14") %>%
  dplyr::select(-ENSG)

expression_slc26a4 <- vsd_slc26a4 %>%
  pivot_longer(cols = everything(), names_to = "sample", values_to = "expression") %>%
  mutate(group = if_else(str_detect(sample, "FNF"), "FNF", "CTL"))

# Create boxplot comparing CTL vs FNF expression
ggplot(expression_slc26a4, aes(x = group, y = expression, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("CTL" = "skyblue", "FNF" = "orange")) +
  labs(x = "Group", y = "Expression (VST)", 
       title = "Expression of SLC26A4 (ENSG00000091137.14)") +
  theme_minimal()

sqtl_slc26a4 <- response_pbs_results_nolmer |> filter(gene_name == "SLC26A4")

create_intron_usage_plot("SLC26A4", sqtl_slc26a4$phe_id, all_psi_inv,all_geno_transpose_df,meta_catl_simple,response_pbs_results)

# Function 

create_gene_expression_boxplot <- function(gene, vsd_geneExp, meta_catl_simple, highConf_resQtL, all_geno_transpose_df) {
  
  # Filter highConf results for the gene of interest.
  # If multiple rows exist, choose the one with the lowest interaction p-value.
  test_boxplotInfo <- highConf_resQtL %>%
    dplyr::filter(gene_name == gene)
  
  if(nrow(test_boxplotInfo) != 1){
    test_boxplotInfo <- test_boxplotInfo %>% slice_min(interaction_pval)
  }
  
  # Extract allele and variant info
  minorAllele <- unique(test_boxplotInfo$minor_allele)
  variantID <- unique(test_boxplotInfo$var_id)
  ensg <- unique(test_boxplotInfo$ensg)

  
  # Process allele information: assume variantID is in the form "chr:pos:ref:alt"
  alleles <- unlist(strsplit(variantID, ":"))
  ref_allele <- alleles[3]
  alt_allele <- alleles[4]
  
  # Define protective allele (the allele not matching the minor allele)
  protective_allele <- ifelse(minorAllele == ref_allele, alt_allele, ref_allele)
  
  # Process gene expression for the gene of interest
  gene_expr <- vsd_geneExp %>%
    dplyr::filter(ENSG == ensg) %>%
    tidyr::pivot_longer(cols = -c(ENSG), names_to = "sampleID", values_to = "expression")
  
  # Process genotype data for the variant based on sample IDs
  geno_data <- all_geno_transpose_df[rownames(all_geno_transpose_df) %in% variantID, ] %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sampleID")
  colnames(geno_data)[2] <- "genotype"
  
  # Join gene expression with genotype and metadata.
  meta_data_fi <- gene_expr %>%
    left_join(geno_data, by = "sampleID") %>%
    left_join(meta_catl_simple, by = c("sampleID" = "ID"))
  
  # Map genotype (originally 0, 1, 2) into allele labels.
  # Here we assume that genotype "0" means homozygous for the protective allele,
  # "2" means homozygous for the minor allele, and "1" is heterozygous.
  meta_combined_all <- meta_data_fi %>%
    mutate(
      genotype = factor(genotype, levels = c("0", "1", "2")),
      genotype = dplyr::case_when(
        genotype == "2" ~ paste(minorAllele, minorAllele, sep = "/"),
        genotype == "1" ~ paste(minorAllele, protective_allele, sep = "/"),
        genotype == "0" ~ paste(protective_allele, protective_allele, sep = "/")
      )
    ) %>%
    # Set levels manually so that protective (0) comes first, then heterozygous (1) then homozygous minor (2)
    mutate(
      genotype = factor(genotype, levels = c(
        paste(protective_allele, protective_allele, sep = "/"),
        paste(minorAllele, protective_allele, sep = "/"),
        paste(minorAllele, minorAllele, sep = "/")
      ))
    )
  
  # Optionally, you can compute y-axis range from the expression data
  maxrange <- range(meta_combined_all$expression, na.rm = TRUE)
  yaxis_range_min <- min(0, round(maxrange[1], 1))
  yaxis_range_max <- round(maxrange[2], 1)
  if (yaxis_range_min == 0 & yaxis_range_max == 0) {
    yaxis_range_min <- maxrange[1]
    yaxis_range_max <- maxrange[2]
  }
  
  # Create the gene expression boxplot
  p <- ggplot(meta_combined_all, aes(x = genotype, y = expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.25, alpha = 0.7) +
    geom_jitter(position = position_jitter(width = 0.2), size = 0.25, color = "grey30") +
    labs(x = "Genotype", y = "Normalized Expression", title = gene_name) +
    scale_fill_manual(values = c('#0067B9', '#FCCE52')) +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(strip.placement = "outside",
          axis.line.y = element_line(linewidth = 0.25),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(color = "black", linewidth = 0.25),
          axis.ticks.length.y = unit(-0.1, "cm"),
          axis.title.x = element_blank(),
          axis.title.y = element_markdown(size = 6, family = "Helvetica", margin = margin(r = -0.5)),
          text = element_text(family = "Helvetica"),
          axis.text.y = element_text(color = "black", size = 6),
          axis.text.x = element_text(color = "black", size = 6, margin = margin(b = 0)),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.background = element_rect(fill = "transparent", color = "transparent"),
          plot.background = element_rect(fill = "transparent", color = "transparent"),
          panel.grid = element_blank(),
          legend.position = "none",
          panel.spacing.x = unit(0.5, "cm"))
  
  return(p)
}

create_gene_expression_boxplot("SLC26A4", vsd_gene_update, meta_catl_simple, 
response_pbs_results_nolmer, all_geno_transpose_df)
