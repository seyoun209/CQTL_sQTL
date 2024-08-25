# Supp table2 - Differentially expressed intron junction in response to FN-f treatment 

library(dplyr)
library(tidyr)
library(formattable)
library(flextable)
library(officer)
library(readxl)
library(writexl)



# This is the Dataset for the PBS vs FNF----------------------------------------
load("output/clu_fnf/introns_fnf_joinAll")
load("./output/clu_oa/introns_oa_joinAll")

process_introns_data <- function(introns_data) {
  introns_processed <- introns_data %>% 
    dplyr::filter(abs(deltapsi_batch) > 0.15) %>%
    dplyr::select(clusterID, gene, ensemblID, chr, start, end, verdict, 
                  deltapsi_batch, phe_id, p.adjust, loglr) %>%
    group_by(clusterID) %>%
    mutate(abs_deltapsi = abs(deltapsi_batch),
           max_deltaPSI = ifelse(abs_deltapsi == max(abs_deltapsi), "Highest", ""),
           max_deltaPSI = ifelse(ensemblID == ".", "", max_deltaPSI)) %>%
    ungroup() %>%
    dplyr::select(-abs_deltapsi) %>%
    dplyr::rename(deltaPSI = deltapsi_batch) %>%
    mutate(chr = factor(chr, levels = paste0("chr",c(1:22, "X", "Y")))) %>%
    arrange(chr, desc(abs(deltaPSI))) %>%
    mutate(chr = as.character(chr)) %>%
    # Rename columns to desired names
    rename(
      "Cluster ID" = clusterID,
      "Gene" = gene,
      "Ensembl ID" = ensemblID,
      "Chromosome" = chr,
      "Start" = start,
      "End" = end,
      "Verdict" = verdict,
      "Delta PSI" = deltaPSI,
      "Phenotype ID" = phe_id,
      "Adjusted P-value" = p.adjust,
      "Log Likelihood Ratio" = loglr,
      "Highest |Delta PSI|" = max_deltaPSI
    )
  
  return(introns_processed)
}

# Process the data
fnf_data <- process_introns_data(introns_fnf_pval_include)
oa_data <- process_introns_data(introns_oa_pval_include)




write_xlsx(list("PBS vs FN-f" = fnf_data,
                "PBS vs OA" = oa_data),
           path = "output/results_plots/Supp_Data/Supplementary_Data2.xlsx")
