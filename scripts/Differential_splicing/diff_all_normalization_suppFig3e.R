# Supplementary Figure 3E Plot - expression
## Author: Seyoun Byun
##Date:06.20.2024
#-------------------------------------------------------------------------------
setwd("/work/users/s/e/seyoun/CQTL_sQTL")
library(data.table)
library(ggtext)
library(ggrepel)
library(colorspace)


meta_cqtl <- fread("output/clu_fnf/meta_cqtl")
meta_cqtl_subset <- meta_cqtl[!meta_cqtl$Donor %in% config$samples_to_omit, ]
meta_cqtl_subset$files <- file.path("output", "quant", meta_cqtl_subset$ID, "quant.sf")
colnames(meta_cqtl_subset) <- gsub("ID", "names", colnames(meta_cqtl_subset))
file.exists(meta_cqtl_subset$files)

## Import data with tximeta & summarize to gene
se_all <- tximeta(meta_cqtl_subset)
gse_all <- summarizeToGene(se_all)


## Convert to factors (avoids a warning)
colData(gse_all)[] <- lapply(colData(gse_all), factor)

## Build DESeq object
dds_gene_all <- DESeqDataSet(gse_all, design = ~ Condition)

## Filter out lowly expressed genes
## (at least 10 counts in at least 22 samples)
keep <- rowSums(counts(dds_gene_all) >= 10) >= ceiling(nrow(colData(gse_all))*0.10)
dds_gene_all <- dds_gene_all[keep,]

## Fit model
dds_gene_all_fit <- DESeq(dds_gene_all)
## Save dds
save(dds_gene_all_fit, file = "output/quant/deGene_allCond_expression_dds.rda")
load("output/quant/deGene_allCond_expression_dds.rda")

de_genes_shrink_fnf <- lfcShrink(dds_gene_all_fit,
                             coef = "Condition_FNF_vs_CTL", format = "GRanges") |>
  plyranges::names_to_column("gene_id")

de_genes_shrink_oa <- lfcShrink(dds_gene_all_fit,
                             coef = "Condition_OA_vs_CTL", format = "GRanges") |>
  plyranges::names_to_column("gene_id")
res_fnf <- lfcShrink(dds_gene_all_fit, coef="Condition_FNF_vs_CTL")
res_oa <- lfcShrink(dds_gene_all_fit, coef="Condition_OA_vs_CTL")


de_genes_shrink_fnf <-
  inner_join(x = as.data.frame(de_genes_shrink_fnf),
             y = as.data.frame(rowData(gse_all)) %>%
               dplyr::select(c("gene_id", "tx_ids")),
             by = "gene_id") %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  keepStandardChromosomes(pruning.mode = "coarse") %>%
  as.data.frame()

de_genes_shrink_oa <-
  inner_join(x = as.data.frame(de_genes_shrink_oa),
             y = as.data.frame(rowData(gse_all)) %>%
               dplyr::select(c("gene_id", "tx_ids")),
             by = "gene_id") %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  keepStandardChromosomes(pruning.mode = "coarse") %>%
  as.data.frame()

write_csv(de_genes_shrink_fnf, file = "output/quant/de_genes_fromAll_fnf_results.csv")
write_csv(de_genes_shrink_oa, file = "output/quant/de_genes_fromAll_oa_results.csv")

## DE stats
snrnp70_de_fnf_l2fc <- read_csv("output/quant/de_genes_fromAll_fnf_results.csv") |>
  filter(gene_id == "ENSG00000104852.15") |>
  pull(log2FoldChange)
snrnp70_de_oa_l2fc <- read_csv("output/quant/de_genes_fromAll_oa_results.csv") |>
  filter(gene_id == "ENSG00000104852.15") |>
  pull(log2FoldChange)

snrnp70_de_fnf_padj <- read_csv("output/quant/de_genes_fromAll_fnf_results.csv") |>
  filter(gene_id == "ENSG00000104852.15") |>
  pull(padj)
snrnp70_de_oa_padj <- read_csv("output/quant/de_genes_fromAll_oa_results.csv") |>
  filter(gene_id == "ENSG00000104852.15") |>
  pull(padj)


SNRNP70_count <- get_gene_condition_Counts("ENSG00000104852.15", dds_gene_all_fit) |> 
  mutate(Condition = case_when(
    Condition == "CTL" ~ "PBS",
    Condition == "FNF" ~ "FN-f",
    TRUE ~ "OA"  # Default for other cases (e.g., "OA" remains unchanged)
  )) |> 
  mutate(Condition = factor(Condition, levels = c("PBS", "FN-f", "OA")))



snrnp70_expBoxplot <- ggplot(SNRNP70_count, aes(x = Condition,
                       y = log2(count),
                       fill = Condition)) +
  geom_boxplot(outlier.shape = NA,
               linewidth = 0.25, alpha = 0.7) +
  geom_jitter(width = 0.2, color = "grey40", size = 0.25)+
  annotate(geom = "segment", x = 1, xend = 2,
           y = 12.75, yend = 12.75, linewidth = 0.25 ,colour ='#FFB81C') +
  annotate(geom = "richtext", label = paste("padj = ",signif(snrnp70_de_fnf_padj, digits = 3),","," ","log~2~FC = ", signif(snrnp70_de_fnf_l2fc, digits = 3)),
           x = 1.5, y = 12.8, family = "Helvetica", size = 2, fill = NA, label.color = NA) +
  annotate(geom = "segment", x = 1, xend = 3,
           y = 12.9, yend = 12.9, linewidth = 0.25,colour='#E07653') +
  annotate(geom = "richtext", label =  paste("padj = ",signif(snrnp70_de_oa_padj, digits = 3),","," ","log~2~FC = ", signif(snrnp70_de_oa_l2fc, digits = 3)),
           x = 2.5, y = 12.95, family = "Helvetica", size = 2, fill = NA, label.color = NA)+
  scale_fill_manual(values = c('#1e87a5','#FFB81C','#E07653')) +
  scale_y_continuous( name = "**SNRNP70** log~2~(normalized counts)",
                      limits = c(11, 13), breaks = seq(11, 13, 0.5)) +
  coord_cartesian(clip = "off") +
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 0.25),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(size = 6, family= "Helvetica",
                                        margin = margin(r = -0.5)),
        text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.text.x = element_text(color = "black", size = 6, margin = margin(b = 0)),
        strip.background = element_blank(),
        strip.text.x.top = element_text(size = 8, margin = margin(b = 5)),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        legend.position = "none",
        panel.spacing.y = unit(0.5, "cm"))

westernblot_SNRNP70 <- data.frame(
  Donor = c("7882", "7882", "7884", "7884", "7885", "7885"),
  Condition = c("PBS", "FN-f", "PBS", "FN-f", "PBS", "FN-f"),
  raw_norm = c(0.96135099, 0.302220352, 0.324304648, 0.2491684, 0.237143421, 0.133916615)
)
westernblot_SNRNP70$Condition <- factor(westernblot_SNRNP70$Condition,levels =c("PBS", "FN-f"))

summary_data <- aggregate(raw_norm ~ Condition, data = westernblot_SNRNP70, 
                          FUN = function(x) c(mean = mean(x), 
                                              se = sd(x) / sqrt(length(x))))
summary_data <- data.frame(
  Condition = summary_data$Condition,
  mean_raw_norm = summary_data$raw_norm[,1],
  se_raw_norm = summary_data$raw_norm[,2]
)

print(summary_data)

SNRNP70_barplot <- ggplot(summary_data, aes(x = Condition, y = mean_raw_norm, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_raw_norm - se_raw_norm, 
                    ymax = mean_raw_norm + se_raw_norm), 
                width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = c('#1e87a5','#FFB81C')) +
  scale_y_continuous( name = "Relative expression of **SNRNP70**",
                      limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  coord_cartesian(clip = "off") +
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 0.25),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(size = 6, family= "Helvetica",
                                        margin = margin(r = -0.5)),
        text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.text.x = element_text(color = "black", size = 6, margin = margin(b = 0)),
        strip.background = element_blank(),
        strip.text.x.top = element_text(size = 8, margin = margin(b = 5)),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        legend.position = "none",
        panel.spacing.y = unit(0.5, "cm"))
SNRNP_box_barplot <- grid.arrange(snrnp70_expBoxplot, SNRNP70_barplot, nrow = 1)
ggsave(filename = "output/results_plots/Supplementary_figures/supp_fig3e_SNRNP70_plot_ex1.pdf",
       plot = SNRNP_box_barplot, width = 5, height = 5, units = "in")
save(SNRNP_box_barplot, file = "output/results_plots/Supplementary_figures/supp_fig3e_SNRNP70_plot_ex1.rda")



ggplot(westernblot_SNRNP70, aes(x = Condition,
                          y = raw_norm,
                          fill = Condition)) +
  geom_jitter(width = 0.2, color = "grey40", size = 0.25)+
  geom_boxplot(outlier.shape = NA,linewidth = 0.25, alpha = 0.7)+
  scale_fill_manual(values = c('#1e87a5','#FFB81C')) +
  scale_y_continuous( name = "SNRNP70/GAPDH",
                      limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  coord_cartesian(clip = "off") +
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 0.25),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(size = 6, family= "Helvetica",
                                        margin = margin(r = -0.5)),
        text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.text.x = element_text(color = "black", size = 6, margin = margin(b = 0)),
        strip.background = element_blank(),
        strip.text.x.top = element_text(size = 8, margin = margin(b = 5)),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        legend.position = "none",
        panel.spacing.y = unit(0.5, "cm"))


ggsave(filename = "output/results_plots/Supplementary_figures/supp_fig3e_SNRNP70boxplot.pdf",
       plot = snrnp70_expBoxplot, width = 3, height = 3, units = "in")
save(snrnp70_expBoxplot, file = "output/results_plots/Supplementary_figures/supp_fig3e_SNRNP70boxplot.rda")

#Individual log2FC
load("output/quant/differential_gene_expression_dds.rda")
normCounts_fnf <- counts(dds_gene, normalized = TRUE)
#load("output/quant/deGene_OA_expression_dds.rda")

fnf_snrnp70_l2fc <- get_sample_l2fc("ENSG00000104852.15",normCounts_fnf) %>%
  dplyr::select("Donor","log2FC") %>%
  mutate(group = 'Gene')






SNRNP70_FC_western <- westernblot_SNRNP70 %>%
  group_by(Donor) %>%
  mutate(log2FC = log2(raw_norm / raw_norm[Condition == "PBS"]))  %>%
  dplyr::filter(Condition  == "FNF") %>%
  dplyr::select("Donor","log2FC") %>%
  mutate(group = 'Protein') 

SNRNP70_FC_norm_fnf <- westernblot_SNRNP70 %>%
  group_by(Donor) %>%
  mutate(log2FC = log2(raw_norm / raw_norm[Condition == "PBS"])) %>%
  dplyr::filter(Condition  == "FNF")
SNRNP70_FC_norm_pbs <- westernblot_SNRNP70 %>%
  dplyr::filter(Condition  == "PBS")
l2fc_SNRNP70_plotdata <- rbind(fnf_snrnp70_l2fc,SNRNP70_FC_western)

logfc <- ggplot(l2fc_SNRNP70_plotdata ,
       mapping = aes(x = group , y = log2FC , fill = group , color = group )) +
  geom_hline(yintercept = 0, lty = 2, color = "grey25", linewidth = 0.25) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.5, alpha = 0.5) +
  geom_jitter(width = 0.2, color = "grey25", size = 0.5) +
  scale_fill_manual(values = c("grey80","grey80"))+
  scale_color_manual(values = c(darken("grey80", 0.3),
                                darken("grey80", 0.3))) +
  scale_y_continuous(name = "**SNRNP70** log~2~ fold change",
                     limits = c(-2, 2), breaks = seq(-2, 2, 0.5)) +
  coord_cartesian(clip = "off") +
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_markdown(size = 6, family = "Helvetica",
                                        margin = margin(r = -15)),
        text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.text.x = element_text(color = "black", size = 8, margin = margin(b = -1)),
        strip.background = element_blank(),
        strip.text.x.bottom = element_markdown(size = 8, margin = margin(t = 1)),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, family = "Helvetica",
                                  size = 10, margin = margin(b=-3)))

ggsave(filename = "output/results_plots/Supplementary_figures/supp_fig3e_SNRNP70baxplot.pdf",
       plot = logfc, width = 3, height = 3, units = "in")
save(logfc, file = "output/results_plots/Supplementary_figures/supp_fig3e_SNRNP70boxplot.rda")


