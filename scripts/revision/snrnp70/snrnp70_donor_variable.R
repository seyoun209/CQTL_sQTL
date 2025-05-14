#To understand the donor varialbe in SNRNP70 both gene and splicing level
setwd("/work/users/s/e/seyoun/CQTL_sQTL/output")
library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(httpgd)
#------------------------------------------------------------------

load("clu_fnf/introns_fnf_joinAll") #introns_fnf_pval_include
introns_fnf_sig <- introns_fnf_pval_include %>% dplyr::filter(p.adjust < 0.05)
introns_fnf_all_sig <- introns_fnf_sig %>% dplyr::filter(abs(deltapsi_batch) > 0.15)

snrnp70_introns_fnf <- introns_fnf_all_sig |>
  filter(genes == "SNRNP70") |>
  arrange(desc(abs(deltapsi)))


#"chr19:49101471:49104634:clu_19784_+

meta_cqtl <- fread("clu_fnf/meta_cqtl")
meta_cqtl_simple <- meta_cqtl %>%
  dplyr::filter(Condition %in% c("CTL","FNF")) %>%
  dplyr::filter(!str_detect(Donor, 'AM7352|AM7244'))

limma_psi_batchremove_junction <- fread("clu_fnf/psi_fnf_limma_batch_corrected")
psi_long <- limma_psi_batchremove_junction %>%
  filter(Junction == "chr19:49101471:49104634:clu_19784_+") %>%
  pivot_longer(cols = -Junction, names_to = "ID", values_to = "psi")

psi_long_meta <- psi_long %>%
  left_join(meta_cqtl_simple, by = "ID")
 

donor_diff <- psi_long_meta %>%
  pivot_wider(id_cols = Donor, names_from = Condition, values_from = psi) %>%
  # Calculate the absolute difference between FNF and CTL
  mutate(diff = abs(FNF - CTL)) %>%
  arrange(desc(diff))
donor_order <- donor_diff$Donor

psi_long_meta <- psi_long_meta %>%
  mutate(Donor = factor(Donor, levels = donor_order))



dodge_width <- 0.3
ggplot(psi_long_meta, aes(x = Donor, y = psi, group = Donor)) +
  # Draw a dashed line connecting CTL and FNF for each donor
  geom_line(position = position_dodge(width = dodge_width),
            color = "grey60", size = 1, linetype = "dashed") +
  geom_point(aes(color = Condition), 
             position = position_dodge(width = dodge_width),
             size = 3, alpha = 0.8) +
  labs(title = "PSI per Donor by Condition (Ordered by PSI Difference)",
       x = "Donor (ordered by PSI difference)",
       y = "PSI Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#To make sure with the gene expression level
#------------------------------------------------------------------
gene_expression_vst <- fread("quant/normalized_vst_gene.txt",sep="\t")

snrnp70_gene_vst <- gene_expression_vst |> filter(ENSG == "ENSG00000104852.15") |> 
pivot_longer(cols = -ENSG, names_to = "ID", values_to = "vst")

vst_long_meta <- snrnp70_gene_vst %>%
  left_join(meta_cqtl_simple, by = "ID")
 

donor_diff <- vst_long_meta %>%
  pivot_wider(id_cols = Donor, names_from = Condition, values_from = vst) %>%
  # Calculate the absolute difference between FNF and CTL
  mutate(diff = abs(FNF - CTL)) %>%
  arrange(desc(diff))
donor_order <- donor_diff$Donor

vst_long_meta <- vst_long_meta %>%
  mutate(Donor = factor(Donor, levels = donor_order))

dodge_width <- 0.3

ggplot(vst_long_meta, aes(x = Donor, y = vst, group = Donor)) +
  # Draw a dashed line connecting CTL and FNF for each donor
  geom_line(position = position_dodge(width = dodge_width),
            color = "grey60", size = 1, linetype = "dashed") +
  geom_point(aes(color = Condition), 
             position = position_dodge(width = dodge_width),
             size = 3, alpha = 0.8) +
  labs(title = "VST per Donor by Condition (Ordered by gene expression Difference)",
       x = "Donor (ordered by gene expression difference)",
       y = "gene expression(vsd) Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





# Calculate PSI differences per donor
donor_diff_psi <- psi_long_meta %>%
  pivot_wider(id_cols = Donor, names_from = Condition, values_from = psi) %>%
  mutate(psi_diff = abs(FNF - CTL)) %>%
  select(Donor, psi_diff)

# Calculate VST differences per donor
donor_diff_vst <- vst_long_meta %>%
  pivot_wider(id_cols = Donor, names_from = Condition, values_from = vst) %>%
  mutate(vst_diff = abs(FNF - CTL)) %>%
  select(Donor, vst_diff)

# Merge the differences by Donor
donor_diff_merge <- left_join(donor_diff_psi, donor_diff_vst, by = "Donor")

# Plot PSI difference vs VST difference to see whether donors with a big PSI change also have a large gene expression change.
ggplot(donor_diff_merge, aes(x = psi_diff, y = vst_diff)) +
  geom_point(size = 3) +
  geom_text(aes(label = Donor), vjust = -0.5, hjust = 0.5, size = 3) +
  labs(x = "PSI Difference |FNF - CTL|", 
       y = "VST Difference |FNF - CTL|", 
       title = "Comparison of PSI and Gene Expression Differences per Donor") +
  theme_minimal()
