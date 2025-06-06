setwd("/work/users/s/e/seyoun/CQTL_sQTL/")
# Load libraries
library(dplyr)
library(stringr)
library(data.table)
library(rtpcr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(scales)
library(ggtext)
library(plotgardener)
library(httpgd)
library(rtracklayer)
library(DESeq2)
library(colorspace)



qpcr_fc_df <- data.frame(
Donor = factor(rep(c("Donor1", "Donor2", "Donor3", "Donor4", "Donor5"), each=2,
levels= c("Donor1", "Donor2", "Donor3", "Donor4", "Donor5"))),
Condition= factor(rep(c("exon7-exon8", "exon7-alt_exon8"), times=5),
levels= c("exon7-exon8", "exon7-alt_exon8")),
FoldChange= c(1.057103233,	0.703820973, 1.287496678,1.184825571,
1.07997656,0.841479482,0.940849134,0.381613984,1.1673521,0.824272456)
) 



qpcr_dotplot <- ggplot(qpcr_fc_df, aes(x = Condition, y = FoldChange)) +
  # Violin plot (with a light fill)
  geom_violin(fill = "grey90", alpha = 0.5, color = NA) +
  # Jittered data points colored by Donor
  geom_jitter(aes(color = Donor), width = 0.15, size = 1.5) +
  # Horizontal lines at the mean for each condition separately
  stat_summary(fun = mean, geom = "errorbar", color = "black", size = 0.25, width = 1, 
               fun.min = function(x) mean(x), fun.max = function(x) mean(x)) +
  # Mean values as text labels for each condition separately
  stat_summary(fun = mean, geom = "text", aes(label = paste("Mean:", round(after_stat(y), 3))), 
               vjust = -0.5, color = "black", size = 2) +
  # Perform paired t-test between conditions and annotate the p-value
  stat_compare_means(method = "t.test", paired = TRUE, label = "p.format",
                     label.x = 1.5, label.y = 1.4, size = 2) +
  # Set y-axis limits and label
  scale_y_continuous(limits = c(0, 1.5)) +
  labs(x = "Condition",
       y = "Fold Change FN-f/PBS (Normalized to junction exon3 - exon 4)") +
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
    legend.position = c(0.9, 0.1),
    legend.background = element_rect(fill = "transparent", color = "transparent"),
    legend.text = element_text(size = 4, family = "Helvetica"),
    panel.spacing.y = unit(0.5, "cm")
  )

save(qpcr_dotplot, 
     file = "output/revision/plots/qpcr_dotplot.rda")



# Fix the wb_df data frame structure for single condition
wb_df <- data.frame(
  Donor = factor(c("Donor1", "Donor2", "Donor3", "Donor4", "Donor5"),
                 levels = c("Donor1", "Donor2", "Donor3", "Donor4", "Donor5")),
  Condition = "Donor",  # Single treatment condition
  FoldChange = c(0.69780367, 0.65358955, 0.44092905, 0.5619581, 0.8677664)
)

# Perform one-sample t-test manually
t_test_result <- t.test(wb_df$FoldChange, mu = 1)
p_value <- format.pval(t_test_result$p.value, digits = 3)

# Create the plot
 westernblot_dotplot <- ggplot(wb_df, aes(x = Condition, y = FoldChange)) +
  # Violin plot (with a light fill)
  geom_violin(fill = "grey90", alpha = 0.5, color = NA) +
  # Jittered data points colored by Donor
  geom_jitter(color="black", width = 0, size = 1.5) +
  # Horizontal line at the mean
  stat_summary(fun = mean, geom = "errorbar", color = "black", size = 0.25, width = 0.5, 
               fun.min = function(x) mean(x), fun.max = function(x) mean(x)) +
  # Mean value as text label
  stat_summary(fun = mean, geom = "text", aes(label = paste("Mean:", round(after_stat(y), 3))), 
               vjust = -0.5, color = "black", size = 2) +
  # Add p-value annotation manually
  annotate("text", x = 1, y = 0.9, label = paste("p =", p_value), size = 2) +
  # Set y-axis limits and label
  scale_y_continuous(limits = c(0, 1.0)) +
  labs(x = "Treatment",
       y = "Fold Change SNRNP70 intensity (FN-f/PBS) normalized to GAPDH") +
  theme_minimal() +
  theme(
    axis.line = element_line(linewidth = 0.25),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black", linewidth = 0.25),
    axis.ticks.length.y = unit(-0.1, "cm"),
    axis.title.x = element_blank(),
    axis.title.y = element_markdown(size = 8, family = "Helvetica", margin = margin(r = 5)),
    text = element_text(family = "Helvetica"),
    axis.text.y = element_text(color = "black", size = 6),
    axis.text.x = element_text(color = "black", size = 6, margin = margin(t = 5)),
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "none", 
    #legend.background = element_rect(fill = "transparent", color = "transparent"),
    #legend.text = element_text(size = 4, family = "Helvetica"),
    panel.spacing.y = unit(0.5, "cm")
  )

save(westernblot_dotplot, 
     file = "output/revision/plots/westernblot_boxplot.rda")
#----------------------------------------------------------------
#Gene expression

load("output/quant/differential_gene_expression_dds.rda") # 

de_genes_shrink <- lfcShrink(dds_gene,
                             coef = "Condition_FNF_vs_CTL", format = "GRanges") |>
  plyranges::names_to_column("gene_id")



# Get the results from DESeq2
results_gene <- results(dds_gene)

# Extract the padj for SNRNP70
snrnp70_padj <- results_gene["ENSG00000104852.15", "padj"]

# Print the result
#cat("SNRNP70 adjusted p-value:", snrnp70_padj, "\n")



meta_cqtl <- fread("output/clu_fnf/meta_cqtl")
meta_cqtl_simple <- meta_cqtl %>%
  dplyr::filter(Condition %in% c("CTL","FNF")) %>%
  dplyr::filter(!str_detect(Donor, 'AM7352|AM7244'))

vsd_gene_expression <- fread("output/quant/normalized_vst_gene.txt")
snrnp70_gene_vst <- vsd_gene_expression |> filter(ENSG == "ENSG00000104852.15") |>
  pivot_longer(cols = -ENSG, names_to = "ID", values_to = "vst_expression")

vst_long_meta <- snrnp70_gene_vst %>%
  left_join(meta_cqtl_simple, by = "ID")

# Change the condition labels
vst_long_meta <- vst_long_meta %>%
  mutate(Condition = case_when(
    Condition == "CTL" ~ "PBS",
    Condition == "FNF" ~ "FN-f",
    TRUE ~ Condition
  ))

# Update the factor levels if needed
vst_long_meta$Condition <- factor(vst_long_meta$Condition, levels = c("PBS", "FN-f"))


# Check the actual direction of change (without absolute value)
#donor_diff_vst_direction <- vst_long_meta %>%
#  pivot_wider(id_cols = Donor, names_from = Condition, values_from = vst_expression) %>%
#  mutate(vst_diff = `FN-f` - PBS) %>%  # Remove abs() to see direction
#  select(Donor, vst_diff)
#donor_info <- donor_diff_vst_direction[donor_diff_vst_direction$vst_diff < 0 ,"Donor"]



snrnp70_expression_boxplot <- ggplot(vst_long_meta, aes(x = Condition, y = vst_expression, fill = Condition)) +
  #geom_hline(yintercept = mean(vst_long_meta$vst _expression), lty = 2, color = "grey25", linewidth = 0.25) +
  geom_jitter(width = 0.2, color = "grey40", size = 0.25) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.5, alpha = 0.4, width = 0.5, 
               color = c('#005587',"#FFB81C")) +
  stat_boxplot(geom = "errorbar", width = 0.5, color = c('#005587',"#FFB81C")) +
  scale_fill_manual(values = c('#005587',"#FFB81C")) +
  scale_y_continuous(name = "VST normalized SNRNP70 Expression",
                     limits = c(min(vst_long_meta$vst_expression) - 0.2, 
                                max(vst_long_meta$vst_expression) + 0.2)) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    strip.placement = "outside",
    axis.line.y = element_line(linewidth = 0.35),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length.y = unit(-0.1, "cm"),
    axis.title.y = element_markdown(size = 8, family = "Helvetica",
                                    margin = margin(r = -15)),
    text = element_text(family = "Helvetica"),
    axis.text.y = element_text(color = "black", size = 6),
    axis.text.x = element_text(size = 8, margin = margin(b = -1),
                               colour = c(darken("#005587", 0.3), darken('#FFB81C', 0.3))),
    strip.background = element_blank(),
    strip.text.x.bottom = element_markdown(size = 8, margin = margin(t = 1)),
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    panel.grid = element_blank(),
    legend.position = "none"
  )

print(snrnp70_expression_boxplot)

save(snrnp70_expression_boxplot, 
     file = "output/revision/plots/snrnp70_expression_boxplot.rda")

#-----------------------------------------------------------------------


# plot gardener to keep both plot at the same time------------------------------
pdf(file = "output/revision/plots/wetlab_plots.pdf",   # The directory you want to save the file in
    width = 6.5, # The width of the plot in inches
    height = 7)

pageCreate(width = 6.5, height =7 , default.units = "inches", showGuides = FALSE)
load("output/revision/plots/snrnp70_expression_boxplot.rda")
plotGG(snrnp70_expression_boxplot, x = 0.25, y = 0.25, width =2.5, height = 3)

load("output/revision/plots/qpcr_dotplot.rda")
plotGG(qpcr_dotplot, x = 3.25, y = 0.25, width = 3, height = 3)

load("output/revision/plots/westernblot_boxplot.rda")
plotGG(westernblot_dotplot, x = 0.25, y = 3.5, width = 2.5, height = 3)


dev.off()
