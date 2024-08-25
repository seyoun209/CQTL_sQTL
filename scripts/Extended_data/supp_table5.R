# This is for the PBS and FN-f lead sQTLs.
setwd("/work/users/s/e/seyoun/CQTL_sQTL/")
library(data.table)
library(openxlsx)
library(dplyr)
library(tidyr)
response_pbs_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_pbs_resInv.rds")
response_fnf_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_fnf_resInv.rds")

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

# Create PBS sQTL table
pbs_sqtl_table <- response_pbs_results  %>%
  dplyr::mutate(
    `Ensembl ID` = ifelse(is.na(ensg) | ensg == "", "Not-annotated", ensg),
    Gene = ifelse(is.na(SYMBOL) | SYMBOL == "", `Ensembl ID`, SYMBOL)
  ) %>%
  select(
    Gene = SYMBOL,
    `Ensembl ID` = ensg,
    intron_junction_id = phe_id,
    intron_junction_chr = phe_chr,
    intron_junction_from = phe_from,
    intron_junction_to = phe_to,
    intron_junction_strand = phe_strd,
    clusterID,
    num_cis_variants = n_var_in_cis,
    rsID,
    variantID = var_id,
    variant_pos = var_from,
    beta = PBS_beta,
    beta_se = PBS_beta_se,
    nom_pval = PBS_p,
    minor_allele,
    MAF,
    nominal_threshold = therehod,
    qval,
    signal = rank
  ) %>%
  arrange(qval)



# Create FNF sQTL table
fnf_sqtl_table <- response_fnf_results %>%
  dplyr::mutate(
    `Ensembl ID` = ifelse(is.na(ensg) | ensg == "", "Not-annotated", ensg),
    Gene = ifelse(is.na(SYMBOL) | SYMBOL == "", `Ensembl ID`, SYMBOL)
  ) %>%
  select(
    Gene = SYMBOL,
    `Ensembl ID` = ensg,
    intron_junction_id = phe_id,
    intron_junction_chr = phe_chr,
    intron_junction_from = phe_from,
    intron_junction_to = phe_to,
    intron_junction_strand = phe_strd,
    clusterID,
    num_cis_variants = n_var_in_cis,
    rsID,
    variantID = var_id,
    variant_pos = var_from,
    beta = FNF_beta,
    beta_se = FNF_beta_se,
    nom_pval = FNF_p,
    minor_allele,
    MAF,
    nominal_threshold = therehod,
    qval,
    signal = rank
  ) %>%
  arrange(qval)


process_data <- function(data) {
  data %>%
    dplyr::filter(interaction_pval < 0.05) %>%
    dplyr::mutate(
      high_confidence = ifelse(abs(delta_beta) > 0.2 & minor_alle_count >= 5, "high-confidence", ""),
      `Ensembl ID` = ifelse(is.na(ensg) | ensg == "", "unavailable", ensg),
      Gene = ifelse(is.na(SYMBOL) | SYMBOL == "", `Ensembl ID`, SYMBOL)
    ) %>%
    dplyr::select(
      Gene,
      `Ensembl ID`,
      intron_junction_id = phe_id,
      intron_junction_chr = phe_chr,
      intron_junction_from = phe_from,
      intron_junction_to = phe_to,
      intron_junction_strand = phe_strd,
      clusterID,
      rsID,
      variantID = var_id,
      variant_pos = var_from,
      PBS_beta,
      PBS_beta_se,
      PBS_nom_pval = PBS_p,
      FNF_beta,
      FNF_beta_se,
      FNF_nom_pval = FNF_p,
      minor_allele,
      MAF,
      interaction_pval,
      high_confidence
    ) %>%
    arrange(desc(high_confidence),interaction_pval)
}


# Process PBS and FNF data
pbs_resQTL <- process_data(response_pbs_results)
fnf_resQTL <- process_data(response_fnf_results)


# Create workbook and add sheets
wb <- createWorkbook()
write_formatted_sheet(wb, "PBS_sQTL", pbs_sqtl_table)
write_formatted_sheet(wb, "FN-f_sQTL", fnf_sqtl_table)
write_formatted_sheet(wb, "PBS_re-sQTL", pbs_resQTL)
write_formatted_sheet(wb, "FN-f_re-sQTL", fnf_resQTL)

saveWorkbook(wb, "output/results_plots/Supp_Data/Supplementary_Data5.xlsx", overwrite = TRUE)
