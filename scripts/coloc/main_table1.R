# Main table coloc
#Plotting for per all the coloc FNF genes
coloc_results_fnf <- fread("output/coloc/coloc_result_table_fnf_prob_calc_.txt")
coloc_results_pbs <- fread("output/coloc/coloc_result_table_pbs_prob_calc_.txt")


# all LD > 0.5
coloc_results_pbs_organized <-  coloc_results_pbs %>%
  mutate(
    clusterID = gsub(".*:(clu_[0-9]+_[+-]).*", "\\1", phe_id),
    `Intron junction` = gsub("(.*):clu_[0-9]+_[+-]", "\\1", phe_id),
    `PBS-(PP H4)/(PP H3 + PP H4)` = formatC(PBS_coloc_H4 / (PBS_coloc_H3 + PBS_coloc_H4)),digits=3,format="f") %>%
  rename(`Boer et al (Cell 2021)` = subtype, 
         `PBS PP H3` = PBS_coloc_H3,
         `PBS PP H4` =  PBS_coloc_H4,
         `FN-f PP H3` = FNF_coloc_H3,
         `FN-f PP H4` =  FNF_coloc_H4) %>%
  select(gene, clusterID,`Intron junction`, clusterID, sQTL_rsID ,sQTL_rank, sQTL_var_id, sQTL_minor_allele, sQTL_MAF,`PBS PP H3`, `PBS PP H4`,`FN-f PP H3`,`FN-f PP H4`,
         `PBS-(PP H4)/(PP H3 + PP H4)`, `Boer et al (Cell 2021)`,  GWAS_rsID ,GWAS_pos, GWAS_ref, GWAS_alt)

coloc_results_fnf_organized <- coloc_results_fnf %>%
  mutate(
    clusterID = gsub(".*:(clu_[0-9]+_[+-]).*", "\\1", phe_id),
    `Intron junction` = gsub("(.*):clu_[0-9]+_[+-]", "\\1", phe_id),
    `FN-f-(PP H4)/(PP H3 + PP H4)` = formatC(FNF_coloc_H4 / (FNF_coloc_H3 + FNF_coloc_H4)),digits=3,format="f") %>%
  rename(`Boer et al (Cell 2021)` = subtype, 
         `PBS PP H3` = PBS_coloc_H3,
         `PBS PP H4` =  PBS_coloc_H4,
         `FN-f PP H3` = FNF_coloc_H3,
         `FN-f PP H4` =  FNF_coloc_H4) %>%
  select(gene, clusterID,`Intron junction`, clusterID, sQTL_rsID ,sQTL_rank, sQTL_var_id, sQTL_minor_allele, sQTL_MAF,`PBS PP H3`, `PBS PP H4`,`FN-f PP H3`,`FN-f PP H4`,
         `FN-f-(PP H4)/(PP H3 + PP H4)`, `Boer et al (Cell 2021)`,  GWAS_rsID ,GWAS_pos, GWAS_ref, GWAS_alt)
# Significant Main table--------------------------------------------------------
create_key <- function(df) {
  paste(df$gene, df$clusterID, df$`Intron junction`, df$sQTL_rsID, sep = "_")
}
coloc_sig_pbs_organized <-  coloc_results_pbs %>% dplyr::filter(PBS_coloc_H4 >= 0.7) |>
  mutate(
    clusterID = gsub(".*:(clu_[0-9]+_[+-]).*", "\\1", phe_id),
    `Intron junction` = gsub("(.*):clu_[0-9]+_[+-]", "\\1", phe_id),
    `(PP H4)/(PP H3 + PP H4)` = formatC(PBS_coloc_H4 / (PBS_coloc_H3 + PBS_coloc_H4)),digits=3,format="f") %>%
  rename(`Boer et al (Cell 2021)` = subtype, 
         `PP H3` = PBS_coloc_H3,
         `PP H4` =  PBS_coloc_H4) %>%
  select(gene, clusterID,`Intron junction`, clusterID, sQTL_rsID , sQTL_var_id, sQTL_minor_allele, sQTL_MAF,`PP H3`,`PP H4`, `(PP H4)/(PP H3 + PP H4)`,
         `Boer et al (Cell 2021)`,  GWAS_rsID ) 


coloc_sig_fnf_organized <- coloc_results_fnf %>% dplyr::filter(FNF_coloc_H4 >= 0.7) |>
  mutate(
    clusterID = gsub(".*:(clu_[0-9]+_[+-]).*", "\\1", phe_id),
    `Intron junction` = gsub("(.*):clu_[0-9]+_[+-]", "\\1", phe_id),
    `(PP H4)/(PP H3 + PP H4)` = formatC(FNF_coloc_H4 / (FNF_coloc_H3 + FNF_coloc_H4)),digits=3,format="f") %>%
  rename(`Boer et al (Cell 2021)` = subtype, 
         `PP H3` = FNF_coloc_H3,
         `PP H4` =  FNF_coloc_H4) %>%
  select(gene, clusterID,`Intron junction`, clusterID, sQTL_rsID , sQTL_var_id, sQTL_minor_allele, sQTL_MAF,`PP H3`,`PP H4`, `(PP H4)/(PP H3 + PP H4)`,
         `Boer et al (Cell 2021)`,  GWAS_rsID)


#coloc_sig_pbs_organized$key <-   create_key(coloc_sig_pbs_organized)
#coloc_sig_fnf_organized$key <- create_key(coloc_sig_fnf_organized)

#pbs_keys <- coloc_sig_pbs_organized$key
#fnf_keys <- coloc_sig_fnf_organized$key

#coloc_sig_pbs_organized$Specificity <- ifelse(coloc_sig_pbs_organized$key %in% fnf_keys, "Shared", "PBS-specific")
#coloc_sig_fnf_organized$Specificity <- ifelse(coloc_sig_fnf_organized$key %in% pbs_keys, "Shared", "FNF-specific")

wb <- createWorkbook()

write_formatted_sheet(wb, "PBS_coloc", coloc_sig_pbs_organized)
write_formatted_sheet(wb, "FN-f_coloc", coloc_sig_fnf_organized)

saveWorkbook(wb, "output/results_plots/coloc/main_table1.xlsx", overwrite = TRUE)
#-------------------------------------------------------------------------------
#Supp table 7
wb <- createWorkbook()
write_formatted_sheet(wb, "PBS_coloc", coloc_results_pbs_organized)
write_formatted_sheet(wb, "FN-f_coloc", coloc_results_fnf_organized)

saveWorkbook(wb, "output/results_plots/Supp_Data/Supplementary_Data7.xlsx.xlsx", overwrite = TRUE)

