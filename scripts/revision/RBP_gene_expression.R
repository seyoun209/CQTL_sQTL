setwd("/work/users/s/e/seyoun/CQTL_sQTL/")
# Load libraries
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(scales)
library(ggtext)
library(plotgardener)
library(httpgd)

load("external_data/encode_rbp/rbp_prep_rbpOnly/pbs_rbp_elip.Rdata") #pbs_rbp_eclip 
load("external_data/encode_rbp/rbp_prep_rbpOnly/fnf_rbp_elip.Rdata") #fnf_rbp_eclip

add_gene_symbols <- function(vsd_genes_exp, org_db = org.Hs.eg.db) {

  # Add gene symbols to DE results
  annots <- AnnotationDbi::select(org_db,
                                  vsd_genes_exp$gene_id,
                                  columns = c("SYMBOL"),
                                  keytype = "ENSEMBL")
  annots <- dplyr::rename(annots, gene_id = ENSEMBL)
  annots_unique <- annots %>% distinct(gene_id, .keep_all = TRUE)
  de_results_with_symbol <- left_join(vsd_genes_exp, annots_unique, by = "gene_id")



  return(de_results_with_symbol)
}


pbs_rbp_eclip_subset <- pbs_rbp_eclip |> group_by(rbp) %>%
  dplyr::filter(observed > 5) |>
  ungroup()  |>
  arrange(desc(odd_med))

fnf_rbp_eclip_subset <- fnf_rbp_eclip |> group_by(rbp) %>%
  dplyr::filter(observed > 5) |>
  ungroup()  |>
  arrange(desc(odd_med))

combined_data_rbp_eclip <- rbind( pbs_rbp_eclip,fnf_rbp_eclip )
all_rbps  <- unique(combined_data_rbp_eclip$rbp)

vsd_gene <- fread("output/quant/normalized_vst_gene.txt") |>
  dplyr::mutate(gene_id = gsub("\\..*$", "" , ENSG))

vsd_geneExp <- add_gene_symbols(vsd_gene)

# Get RBP expression data
  rbp_exp_data <- vsd_geneExp %>%
    filter(SYMBOL %in% all_rbps) 

rbp_exp_data_numeric <- rbp_exp_data %>% 
  select_if(is.numeric)

rbp_exp_data <- rbp_exp_data %>%
  mutate(expression_sum =  rowMeans(rbp_exp_data_numeric, na.rm = TRUE))


ggplot(rbp_exp_data, aes(x = expression_sum)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  labs(x = "Row Mean Expression",
       y = "Frequency",
       title = "Distribution of Row Mean Expression for RBPs") +
    theme_minimal()

combined_rbp_data <- combine_rbp_data(pbs_rbp_corr_all, fnf_rbp_corr_all)

combined_rbp_data_dist <- combined_rbp_data %>%
  mutate(
    chromosome = str_extract(variantID, "chr\\d+|chrX|chrY"),
    position = as.numeric(str_extract(variantID, "(?<=:)\\d+")),
    intron_chr = str_extract(Intron_junction_id, "chr\\d+|chrX|chrY"),
    intron_start = as.numeric(str_extract(Intron_junction_id, "(?<=:)\\d+")),
    intron_end = as.numeric(str_split_fixed(Intron_junction_id, ":", 4)[,3]),
    clu_id = str_extract(Intron_junction_id, "clu_\\d+"),
    strand = str_extract(Intron_junction_id, "[+-]$"),
    distance_from_start = ifelse(chromosome == intron_chr, 
                                 abs(position - intron_start), NA),
    distance_from_end = ifelse(chromosome == intron_chr, 
                               abs(position - intron_end), NA),
    dist_intron_junction_var = pmin(distance_from_start, distance_from_end, na.rm = TRUE)
  ) %>%
  dplyr::select(Gene, `Ensembl ID`, Intron_junction_id, ClusterID, dist_intron_junction_var, 
         rsID, variantID, Significant_overlapping_RBPs, Source)


combined_rbp_data_onlySig <- combined_rbp_data_dist %>%
  filter(Significant_overlapping_RBPs != "") 


sig_enriched_ensg <- combined_rbp_data_onlySig %>%
  distinct(`Ensembl ID`) %>%
  filter(!is.na(`Ensembl ID`)) %>%
  pull(`Ensembl ID`)


