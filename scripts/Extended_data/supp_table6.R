#Supp table 6 for the Pathway enriched sQTL binding to RNA binding sites ( RBP overall binding table (Overall output with the var_id, phe_id , source(pbs or FNF)),
#fenrich output (PBS, FNF), Pathway of the sGenes associated with RBP)

library(dplyr)
library(tidyverse)
library(data.table)
library(openxlsx)


load("output/Enrichment/fnf_results_rbp_overlapping.rda")
load("output/Enrichment/pbs_results_rbp_overlapping.rda")

# pbs_data_df <- pbs_results %>%
#   mutate(overlapping_rbps = map_chr(overlapping_rbps, ~ {
#     # Remove the c() wrapper and quotes
#     clean_string <- str_remove_all(.x, "^c\\(|\\)$|\"|'")
#     # Split the string, trim whitespace, and rejoin with commas
#     str_split(clean_string, ",\\s*") %>% 
#       unlist() %>% 
#       str_trim() %>% 
#       paste(collapse = ", ")
#   })) |> data.frame()
# 
# 
# fnf_data_df <- fnf_results %>%
#   mutate(overlapping_rbps = map_chr(overlapping_rbps, ~ {
#     # Remove the c() wrapper and quotes
#     clean_string <- str_remove_all(.x, "^c\\(|\\)$|\"|'")
#     # Split the string, trim whitespace, and rejoin with commas
#     str_split(clean_string, ",\\s*") %>% 
#       unlist() %>% 
#       str_trim() %>% 
#       paste(collapse = ", ")
#   })) |> data.frame()
# 
# 
# 
# # In this excel, I need a list of RBP , Var_id, phe_id ,  var_chr  var_from, ensg, SYMBOL.y
# 
# fnf_rbp_overlaps <- fnf_data_df  %>%
#   dplyr::mutate(
#     `Ensembl ID` = ifelse(is.na(ensg) | ensg == "", "Not-annotated", ensg),
#     Gene = ifelse(is.na(SYMBOL.y) | SYMBOL.y == "", `Ensembl ID`, SYMBOL.y)
#   ) %>%
#   select(
#     Gene = SYMBOL.y,
#     `Ensembl ID` = ensg,
#     intron_junction_id = phe_id,
#     clusterID,
#     rsID,
#     variantID = var_id,
#     variant_chr= var_chr,
#     variant_pos = var_from,
#     RBP_sQTL_Overlaps= overlapping_rbps,
#     minor_allele,
#     MAF,
#     qval
#   ) %>%
#   arrange(qval)
# 
# 
# pbs_rbp_overlaps <- pbs_data_df  %>%
#   dplyr::mutate(
#     `Ensembl ID` = ifelse(is.na(ensg) | ensg == "", "Not-annotated", ensg),
#     Gene = ifelse(is.na(SYMBOL.y) | SYMBOL.y == "", `Ensembl ID`, SYMBOL.y)
#   ) %>%
#   select(
#     Gene = SYMBOL.y,
#     `Ensembl ID` = ensg,
#     intron_junction_id = phe_id,
#     clusterID,
#     rsID,
#     variantID = var_id,
#     variant_chr= var_chr,
#     variant_pos = var_from,
#     RBP_sQTL_Overlaps= overlapping_rbps,
#     minor_allele,
#     MAF,
#     qval
#   ) %>%
#   arrange(qval)


write_formatted_sheet <- function(wb, sheet_name, data) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, data, startRow = 1, startCol = 1)
  
  headerStyle <- createStyle(textDecoration = "bold", halign = "center", valign = "center", fgFill = "grey90")
  dataStyle <- createStyle(wrapText = TRUE, halign = "left", valign = "center")
  
  addStyle(wb, sheet_name, headerStyle, rows = 1, cols = 1:ncol(data), gridExpand = TRUE)
  addStyle(wb, sheet_name, dataStyle, rows = 2:(nrow(data)+1), cols = 1:ncol(data), gridExpand = TRUE)
  
  setColWidths(wb, sheet_name, cols = 1:ncol(data), widths = "auto")
  freezePane(wb, sheet_name, firstRow = TRUE)
}


# Function to prepare data (rename columns and select subset)
prepare_data <- function(data) {
  data %>%
    dplyr::rename(
      Observed = observed,
      Total_QTLs = total_qtls,
      Expected_Mean = expected_mean,
      Expected_SD = expected_sd,
      RBP = rbp,
      Empirical_P = epval,
      OR_Lower = odd_dnv,
      OR_Median = odd_med,
      OR_Upper = odd_upv
    ) %>%
    dplyr::select(Observed, Total_QTLs, Expected_Mean, Expected_SD, RBP,  Empirical_P, OR_Lower, OR_Median, OR_Upper)
}




load("external_data/encode_rbp/rbp_prep_rbpOnly/pbs_rbp_elip.Rdata") #pbs_rbp_eclip 
load("external_data/encode_rbp/rbp_prep_rbpOnly/fnf_rbp_elip.Rdata") #fnf_rbp_eclip

pbs_rbp_eclip_subset <- pbs_rbp_eclip |> group_by(rbp) %>%
  dplyr::filter(observed > 5) |>
  ungroup()  |>
  arrange(desc(odd_med))

fnf_rbp_eclip_subset <- fnf_rbp_eclip |> group_by(rbp) %>%
  dplyr::filter(observed > 5) |>
  ungroup()  |>
  arrange(desc(odd_med))





#Make the overall RBPs
load("output/Enrichment/pbs_results_allCorrelation_info.rda") # pbs_results_allCorrelation_info
load("output/Enrichment/fnf_results_allCorrelation_info.rda") #fnf_results_allCorrelation_info


prepare_rbp_data <- function(data, dataType) {
  data %>%
    mutate(Source = rep(dataType)) %>%
    dplyr::select(SYMBOL, ensg, phe_id, clusterID, rsID, var_id, all_overlapping_rbps, significant_overlapping_rbps, Source) %>% 
    dplyr::rename(
      Gene = SYMBOL,
      `Ensembl ID` = ensg,
      Intron_junction_id = phe_id,
      ClusterID = clusterID,
      variantID = var_id,
      sQTL_overlapping_RBPs = all_overlapping_rbps,
      Significant_overlapping_RBPs = significant_overlapping_rbps
    )
}

pbs_rbp_corr_all <- prepare_rbp_data(pbs_results_allCorrelation_info, "PBS")
fnf_rbp_corr_all <- prepare_rbp_data(fnf_results_allCorrelation_info, "FN-f")


combine_rbp_data <- function(pbs_data, fnf_data) {
  # Combine the datasets
  combined_rbp_data <- bind_rows(pbs_data, fnf_data)
  
  # Identify duplicates based on Intron_junction_id, ClusterID, rsID, and variantID
  duplicates <- combined_rbp_data %>%
    group_by(Intron_junction_id, ClusterID, rsID, variantID) %>%
    filter(n() > 1) %>%
    ungroup()
  
  # For duplicates, combine sQTL_overlapping_RBPs and Significant_overlapping_RBPs
  duplicates_combined <- duplicates %>%
    group_by(Intron_junction_id, ClusterID, rsID, variantID) %>%
    summarise(
      Gene = unique(Gene)[1],
      `Ensembl ID` = unique(`Ensembl ID`)[1],
      sQTL_overlapping_RBPs = paste(unique(sQTL_overlapping_RBPs), collapse = "; "),
      Significant_overlapping_RBPs = paste(unique(Significant_overlapping_RBPs), collapse = "; "),
      Source = "Both",
      .groups = "drop"
    )
  
  # Remove duplicates from the original combined data
  combined_data_unique <- combined_rbp_data %>%
    anti_join(duplicates, by = c("Intron_junction_id", "ClusterID", "rsID", "variantID","sQTL_overlapping_RBPs"))
  
  # Combine the unique data with the combined duplicates
  final_data <- bind_rows(combined_data_unique, duplicates_combined) %>%
    arrange(Intron_junction_id, ClusterID, rsID, variantID)
  
  return(final_data)
}


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
         rsID, variantID, Significant_overlapping_RBPs)

combined_rbp_data_onlySig <- combined_rbp_data_dist %>%
  filter(Significant_overlapping_RBPs != "") 


sig_enriched_ensg <- combined_rbp_data_onlySig %>%
  distinct(`Ensembl ID`) %>%
  filter(!is.na(`Ensembl ID`)) %>%
  pull(`Ensembl ID`)

write.table(sig_enriched_ensg, file = "output/Enrichment/sig_enriched_ensg.txt",sep='\t',quote=F,row.names=F,col.names=F)

#system("scripts/enrichment/run_homer.sh /work/users/s/e/seyoun/CQTL_sQTL/output/Enrichment/sig_enriched_ensg.txt /work/users/s/e/seyoun/CQTL_sQTL/output/Enrichment/homer/homer_sig_enriched")



##pathway Reacome + kegg
# Read in from Homer
reactome_data <- read_delim("output/Enrichment/homer/homer_sig_enriched/reactome.txt") |>
  mutate(pval = exp(1)^logP) |>
  dplyr::filter(pval < 0.01) |>
  distinct(Term, .keep_all = TRUE) |>
  mutate(`-log10pval` = -log10(pval)) |>
  mutate(category = "Reactome Pathway")

kegg_data <- read_delim("output/Enrichment/homer/homer_sig_enriched/kegg.txt") |>
  mutate(pval = exp(1)^logP) |>
  dplyr::filter(pval < 0.01) |>
  distinct(Term, .keep_all = TRUE) |>
  mutate(`-log10pval` = -log10(pval)) |>
  mutate(category = "KEGG Pathway")

combined_pathway <- bind_rows(reactome_data, kegg_data)

## Format and write to table
pathway_table <- combined_pathway |>
  dplyr::select(-logP, -pval, -`Entrez Gene IDs`) |>
  relocate(`-log10pval`, .after = Enrichment) |>
  arrange(desc(`-log10pval`))
write_csv(pathway_table, file = "output/Enrichment/table/pathway_table_sig.csv")






# Create a new workbook
wb <- createWorkbook()

write_formatted_sheet(wb, "PBS_RBP_enrichment", prepare_data(pbs_rbp_eclip_subset))
write_formatted_sheet(wb, "FN-f_RBP_enrichment", prepare_data(fnf_rbp_eclip_subset))
write_formatted_sheet(wb, "sQTL_RBP_sig_correlation", combined_rbp_data_onlySig)
write_formatted_sheet(wb, "Pathway (sQTL_RBP_sig)",pathway_table)

# Save the workbook

saveWorkbook(wb, "output/results_plots/Supp_Data/Supplementary_Data6.xlsx", overwrite = TRUE)

