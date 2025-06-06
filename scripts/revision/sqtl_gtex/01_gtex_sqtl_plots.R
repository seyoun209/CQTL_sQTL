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

gtex_dir <- "external_data/GTEx_v10_sQTL/GTEx_Analysis_v10_sQTL_updated"
#------------------------------------------------------------------------------------
# Load the datas
## load the sQTL results
response_pbs_results <- readRDS("/work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/conditional_pbs/response_pbs_resInv.rds")
response_fnf_results <- readRDS("/work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/conditional_pbs/response_fnf_resInv.rds")
load("/work/users/s/e/seyoun/CQTL_sQTL/output/revision/data/comparison_gtex.rdata")

# Makeing a table
pbs_gtex_summary_filtered <- pbs_comparison_gtex$all_matches %>%
  filter(phe_id == gtex_phe_id) %>% 
  # Group by chondrocyte sQTL (var_id, phe_id) and by each tissue separately:
  group_by(var_id, phe_id, tissue) %>%
  summarize(
    ld_info = paste(unique(paste0(ld_variant, " (R2=", ld_r2, ")")), collapse = "; "),
    .groups = "drop"
  ) %>%
  # Now, for each sQTL group, aggregate the tissues and their ld_info together.
  group_by(var_id, phe_id) %>%
  summarize(
    tissues = paste(sort(unique(tissue)), collapse = ", "),
    tissue_count = n_distinct(tissue),
    ld_info = paste(paste0(tissue, ": ", ld_info), collapse = " | "),
    .groups = "drop"
  ) %>%
  # Rename phe_id to genomicLoc for the join
  rename(genomicLoc = phe_id)
# Merge the summarized GTEx info back to your chondrocyte sQTL results.
pbs_final_table <- response_pbs_results %>%
  left_join(pbs_gtex_summary_filtered, by = c("var_id", "genomicLoc"))

gtex_exist_pbs <- pbs_final_table |> filter(!is.na(tissue_count)) 


fnf_gtex_summary_filtered <- fnf_comparison_gtex$all_matches %>%
  filter(phe_id == gtex_phe_id) %>% 
  # Group by chondrocyte sQTL (var_id, phe_id) and by each tissue separately:
  group_by(var_id, phe_id, tissue) %>%
  summarize(
    ld_info = paste(unique(paste0(ld_variant, " (R2=", ld_r2, ")")), collapse = "; "),
    .groups = "drop"
  ) %>%
  # Now, for each sQTL group, aggregate the tissues and their ld_info together.
  group_by(var_id, phe_id) %>%
  summarize(
    tissues = paste(sort(unique(tissue)), collapse = ", "),
    tissue_count = n_distinct(tissue),
    ld_info = paste(paste0(tissue, ": ", ld_info), collapse = " | "),
    .groups = "drop"
  ) %>%
  # Rename phe_id to genomicLoc for the join
  rename(genomicLoc = phe_id)
# Merge the summarized GTEx info back to your chondrocyte sQTL results.
fnf_final_table <- response_fnf_results %>%
  left_join(fnf_gtex_summary_filtered, by = c("var_id", "genomicLoc"))

gtex_exist_fnf <- fnf_final_table |> filter(!is.na(tissue_count)) 

#------------------------------------------------------------------------------------
#Table 1 Colocalization between sQTL and OA GWAS signals
#pbs_coloc_sQTL_gtex <- pbs_final_table |> 
#filter(var_id == c("chr6:18399163:C:T","chr16:69921902:T:C","chr1:184036994:A:G","chr22:37792277:C:A","chr21:39347315:C:T"))

#fnf_coloc_sQTL_gtex <- gtex_exist_fnf |> 
#filter(var_id == c("chr6:18399163:C:T","chr22:37808252:C:T","chr3:52594305:A:T"))


coloc_results_pbs <- fread("output/coloc/coloc_result_table_pbs_prob_calc_.txt")
coloc_results_fnf <- fread("output/coloc/coloc_result_table_fnf_prob_calc_.txt")


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

# For the PBS colocalization results
coloc_results_pbs_organized <- coloc_results_pbs_organized %>%
  left_join(
    pbs_final_table %>% select(var_id, genomicLoc, tissue_count, tissues),
    by = c("sQTL_var_id" = "var_id", "Intron junction" = "genomicLoc")
  ) %>%
  mutate(`GTEx tissue present` = if_else(!is.na(tissue_count) & tissue_count > 0, "Yes", "No"))

# For the FNF colocalization results
coloc_results_fnf_organized <- coloc_results_fnf_organized %>%
  left_join(
    pbs_final_table %>% select(var_id, genomicLoc, tissue_count, tissues),
    by = c("sQTL_var_id" = "var_id", "Intron junction" = "genomicLoc")
  ) %>%
  mutate(`GTEx tissue present` = if_else(!is.na(tissue_count) & tissue_count > 0, "Yes", "No"))

coloc_results_pbs_organized |> filter(`PBS PP H4` > 0.7) |> select(gene, `Intron junction`,  sQTL_rsID,tissue_count)
coloc_results_fnf_organized |> filter(`FN-f PP H4` > 0.7) |> select(gene, `Intron junction`,  sQTL_rsID,tissue_count)
#------------------------------------------------------------------------------------
#Update Supp table 7
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
#-------------------------------------------------------------------------------------
wb <- createWorkbook()
write_formatted_sheet(wb, "PBS_coloc", coloc_results_pbs_organized)
write_formatted_sheet(wb, "FN-f_coloc", coloc_results_fnf_organized)

saveWorkbook(wb, "output/results_plots/Supp_Data/Supplementary_Data7_v2.xlsx.xlsx", overwrite = TRUE)


#------------------------------------------------------------------------------------
# Make the barplot?

# For PBS:
pbs_tissue_counts_by_rank <- pbs_final_table %>%
  mutate(tissue_count = if_else(is.na(tissue_count), 0L, tissue_count)) %>% 
  group_by(rank, tissue_count) %>%
  summarise(freq = n(), .groups = "drop")

# Convert rank to character with degree symbol if desired (optional)
pbs_tissue_counts_by_rank <- pbs_final_table %>%
  mutate(tissue_count = if_else(is.na(tissue_count), 0L, tissue_count)) %>% 
  group_by(rank, tissue_count) %>%
  summarise(freq = n(), .groups = "drop") %>%
  # Create a label for rank with the degree symbol
  mutate(rank_lbl = case_when(
    rank == 0 ~ "0°",
    rank == 1 ~ "1°",
    rank == 2 ~ "2°",
    rank == 3 ~ "3°",
    TRUE ~ as.character(rank)
  ))


ggplot(pbs_tissue_counts_by_rank, aes(x = factor(tissue_count), y = freq)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  #scale_fill_manual(values = c("0°" = "#A5DEF2", "1°" = "#5BAEB7", "2°" = "#1E80C1", "3°" = "#414C6B")) +
  scale_y_continuous(
    name = "Counts of sQTL-splice intron junction",
    trans = pseudo_log_trans(base = 2),
    breaks = c(0, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096),
    labels = comma,
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_x_discrete(
    name = "GTEx tissue count",
    drop = FALSE
  ) +
  coord_cartesian(clip = "off") +
  labs(title = "Distribution of GTEx Tissue Count by Rank (PBS)") +
  theme(
    strip.placement = "outside",
    axis.line = element_line(linewidth = 0.25),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black", linewidth = 0.25),
    axis.ticks.length.y = unit(-0.1, "cm"),
    axis.title.x = element_markdown(size = 8, family = "Helvetica", margin = margin(t = 5)),
    axis.title.y = element_markdown(size = 8, family = "Helvetica", margin = margin(r = 5)),
    text = element_text(family = "Helvetica"),
    axis.text.y = element_text(color = "black", size = 6),
    axis.text.x = element_text(color = "black", size = 6, margin = margin(t = 5)),
    strip.background = element_blank(),
    strip.text.x.top = element_text(size = 8, margin = margin(b = 5)),
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.9, 0.9),
    panel.spacing.y = unit(0.5, "cm")
  )

#-------------------------------------------------------------------------------------
# This is PBS and FNF together of GTEx tissue counts (ignored the conditional)

# Summarize PBS data: replace NA with 0, add Group column, then count frequency per tissue_count
pbs_summary <- pbs_final_table %>%
  mutate(tissue_count = if_else(is.na(tissue_count), 0L, tissue_count),
         Group = "PBS") %>%
  group_by(tissue_count, Group) %>%
  summarise(freq = n(), .groups = "drop")

# Summarize FN-f data similarly
fnf_summary <- fnf_final_table %>%
  mutate(tissue_count = if_else(is.na(tissue_count), 0L, tissue_count),
         Group = "FN-f") %>%
  group_by(tissue_count, Group) %>%
  summarise(freq = n(), .groups = "drop")

# Combine the summaries
all_summary <- bind_rows(pbs_summary, fnf_summary)
all_summary <- all_summary |> mutate(Group = factor(Group, levels = c("PBS", "FN-f")))
# Create the barplot
gtex_sqtl_count_Barplot <- ggplot(all_summary, aes(x = factor(tissue_count), y = freq, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  scale_fill_manual(values = c("PBS" = "#BFDDFF", "FN-f" = "#FFDDA2")) +
  scale_y_continuous(
    name = "Counts of sQTLs",
    trans = pseudo_log_trans(base = 2),
    breaks = c(0, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096),
    labels = comma,
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_x_discrete(
    name = "GTEx tissue count",
    drop = FALSE
  ) +
  coord_cartesian(clip = "off") +
  labs(title = "Distribution of GTEx Tissue Count by Dataset") +
  theme_minimal() +
  theme(
    axis.line = element_line(linewidth = 0.25),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black", linewidth = 0.25),
    axis.ticks.length.y = unit(-0.1, "cm"),
    axis.title.x = element_markdown(size = 8, family = "Helvetica", margin = margin(t = 5)),
    axis.title.y = element_markdown(size = 8, family = "Helvetica", margin = margin(r = 5)),
    text = element_text(family = "Helvetica"),
    axis.text.y = element_text(color = "black", size = 6),
    axis.text.x = element_text(color = "black", size = 6, margin = margin(t = 5)),
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.9, 0.9),
    panel.spacing.y = unit(0.5, "cm")
  )

save(gtex_sqtl_count_Barplot,
file = "output/revision/plots/gtex_sqtl_count_Barplot.rda")


#---------------------------------------------------------
# Make the barplot with the bin (0, 1-5, 5-15, 15-50)

all_summary_binned <- all_summary %>%
  mutate(tissue_bin = case_when(
    tissue_count == 0 ~ "No overlaps",
    tissue_count >= 1 & tissue_count <= 5 ~ "Overlaps of 1-5 tissues",
    tissue_count >= 6 & tissue_count <= 15 ~ "Overlaps of 5-15 tissues",
    tissue_count >= 16 ~ "Overlaps of 15-50 tissues"  # adjust as needed
  )) %>%
  mutate(tissue_bin = factor(tissue_bin, levels = c("No overlaps", "Overlaps of 1-5 tissues", "Overlaps of 5-15 tissues", "Overlaps of 15-50 tissues")))

# Aggregate by tissue_bin and Group:
bin_summary <- all_summary_binned %>%
  group_by(Group, tissue_bin) %>%
  summarise(freq = sum(freq), .groups = "drop")


Gtex_binsize_barplot <- ggplot(bin_summary, aes(x = factor(tissue_bin), y = freq, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  scale_fill_manual(values = c("PBS" = "#BFDDFF", "FN-f" = "#FFDDA2")) +
  scale_y_continuous(
    name = "Counts of sQTL",
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.1))
  ) +
  geom_text(aes(label = freq),
            position = position_dodge(width = 0.9),
            vjust = -0.5,
            size = 2) +
  #scale_y_continuous(
  #  name = "Counts of sQTLs",
  #  trans = pseudo_log_trans(base = 2),
  #  breaks = c(0, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096),
  #  labels = comma,
  #  expand = expansion(mult = c(0, 0.1))
  #) +
  scale_x_discrete(
    name = "GTEx tissue count",
    drop = FALSE
  ) +
  coord_cartesian(clip = "off") +
  labs(title = "Distribution of GTEx Tissue Count by Dataset") +
  theme_minimal() +
  theme(
    axis.line = element_line(linewidth = 0.25),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black", linewidth = 0.25),
    axis.ticks.length.y = unit(-0.1, "cm"),
    axis.title.x = element_markdown(size = 8, family = "Helvetica", margin = margin(t = 5)),
    axis.title.y = element_markdown(size = 8, family = "Helvetica", margin = margin(r = 5)),
    text = element_text(family = "Helvetica"),
    axis.text.y = element_text(color = "black", size = 6),
    axis.text.x = element_text(color = "black", size = 6, margin = margin(t = 5)),
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.9, 0.9),
    legend.key.size = unit(0.2, "cm"),
    legend.text = element_text(size = 6),
    panel.spacing.y = unit(0.5, "cm"),
    title = element_blank()
  )

save(Gtex_binsize_barplot,
     file = "output/revision/plots/gtex_sqtl_binsize_Barplot.rda")
# Barplot in plotgardener

# plot gardener to keep both plot at the same time------------------------------
pdf(file = "output/revision/plots/gtex_sqtl_count_bin_Barplot.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 3.5)

pageCreate(width = 5, height =3.5 , default.units = "inches", showGuides = FALSE)
load("output/revision/plots/gtex_sqtl_binsize_Barplot.rda")
plotGG(Gtex_binsize_barplot, x = 0.4, y = 0.5, width = 4.5, height = 3)
dev.off()
