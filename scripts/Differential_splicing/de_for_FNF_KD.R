#Finding the how similar the test run vs chondro RNA-seq (Splicing) and gene expression
## Author: Seyoun Byun
## Date: 03.20.2024
## Edited:
#-------------------------------------------------------------------------------
setwd("/work/users/s/e/seyoun/CQTL_sQTL/")
library(limma)
library(magrittr)
library(data.table)
library(dplyr)
library(leafviz)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(RColorBrewer)
library(plotgardener)
library(grid)
library(ggvenn)
library(biomaRt)
library(ggtext)


diff_fnf <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/quant/sig_deGenes_pval01_l2fc15.csv")
diff_kd <- fread("/work/users/s/e/seyoun/crispr/02.test_seq/condition_de/sig_deGenes_pval_01_l2fc15.csv")


downsig_deGens_fnf <- diff_fnf |> 
  filter(log2FoldChange < 0) 

upsig_deGenes_fnf <- diff_fnf |> 
  filter(log2FoldChange > 0) 

downsig_deGens_kd <- diff_kd |> 
  filter(log2FoldChange > 0) 

upsig_deGenes_kd <- diff_kd |> 
  filter(log2FoldChange < 0) 

#Genes
downsigGenes_fnf <- downsig_deGens_fnf$gene_id
downsigGenes_KD <- downsig_deGens_kd$gene_id
upsigGenes_fnf <- upsig_deGenes_fnf$gene_id
upsigGenes_KD <- upsig_deGenes_kd$gene_id


down_list <- list("PBSvsFNF"=downsigGenes_fnf, "WDvsKD"=downsigGenes_KD)
up_list <- list("PBSvsFNF"=upsigGenes_fnf, "WDvsKD"=upsigGenes_KD)
down_venn <-ggvenn(
  down_list, 
  fill_color =   c("#FDCDAC", "#CBD5E8"),
  stroke_size = 0.5, set_name_size = 4,
  show_outside = "none",
)

up_venn <-ggvenn(
  up_list, 
  fill_color =   c("#FDCDAC", "#CBD5E8"),
  stroke_size = 0.5, set_name_size = 4,show_outside = "none",
)


# Assuming up_venn and down_venn are ggplot objects
up_venn <- up_venn + theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5, "points")) +
  ggtitle("Differentially up-regulated genes") +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2.5))

down_venn <- down_venn + theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5, "points")) +
  ggtitle("Differentially down-regulated genes") +
  theme(plot.title = element_text(hjust = 0.5, vjust =2.5))
pdf(file = "output/results_plots/deGens_venndiagram.pdf",   # The directory you want to save the file in
    width = 11, # The width of the plot in inches
    height = 6)
# Use grid.arrange to plot them together
grid.arrange(up_venn, down_venn, ncol = 2)


#------------------------------------------------------------------------------
#down- fingind hgnc
intersect_downGenes <- intersect(downsigGenes_fnf, downsigGenes_KD)
intersect_downGenes_noVer <- sapply(strsplit(intersect_downGenes, "\\."), `[`, 1)
#up
intersect_upGenes <- intersect(upsigGenes_fnf, upsigGenes_KD)
intersect_upGenes_noVer <- sapply(strsplit(intersect_upGenes, "\\."), `[`, 1)

# Get HGNC symbols
down_hgnc_mapping <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                      filters = 'ensembl_gene_id', 
                      values = intersect_downGenes_noVer, 
                      mart = ensembl)
up_hgnc_mapping <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                      filters = 'ensembl_gene_id', 
                      values = intersect_upGenes_noVer, 
                      mart = ensembl)



#------------------------------------------------------------------------------
#Now, try to understand truly because of SNRNP70? 
setdiff_down_kd <- setdiff(downsigGenes_KD,downsigGenes_fnf)
setdiff_down_kd_noVer <- sapply(strsplit(setdiff_down_kd, "\\."), `[`, 1)
setdiff_up_kd <- setdiff( upsigGenes_KD,upsigGenes_fnf)
setdiff_up_kd_noVer <- sapply(strsplit(setdiff_up_kd, "\\."), `[`, 1)



# Get HGNC symbols
down_kd_hgnc <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                           filters = 'ensembl_gene_id', 
                           values = setdiff_down_kd_noVer, 
                           mart = ensembl)|> write_csv("/work/users/s/e/seyoun/crispr/02.test_seq/condition_de/down_KD_only_gene.csv")
up_kd_hgnc <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                         filters = 'ensembl_gene_id', 
                         values = setdiff_up_kd_noVer, 
                         mart = ensembl) |> write_csv("/work/users/s/e/seyoun/crispr/02.test_seq/condition_de/up_KD_only_gene.csv")

#Pathway (Upregulated)
up_KD_pathway <- fread()
down_KD_pathway <- fread("/work/users/s/e/seyoun/crispr/02.test_seq/condition_de/pathway_go_results/down_kd_only_pathway.txt")

downSig_kd_pathway_data <- down_KD_pathway %>%
  filter(`q-value` <0.05) %>%
  distinct(pathway, .keep_all = TRUE) %>%
  mutate(`-log10(q-value)` = -log10(`q-value`)) %>%
  mutate(category ="Downregulated")

kegg_plotting <- downSig_kd_pathway_data %>% arrange(`-log10(q-value)`)

#Geneonoloty
down_KD_GO <- fread("/work/users/s/e/seyoun/crispr/02.test_seq/condition_de/pathway_go_results/down_kd_only_go.txt")


ggplot(kegg_plotting, aes(x = `-log10(q-value)`, y = pathway, fill = category)) +
  geom_vline(xintercept = 5, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 10, color = "grey75", alpha = 0.4) +
  geom_bar(stat = "identity") +
  scale_x_continuous(expand = c(0, 0), name = "-log~10~qval", limits = c(0, 10),
                     breaks = seq(0, 25, 5)) +
  scale_fill_manual(values = c("#FBBE67",  "#78A1Cd")) +
  facet_wrap(~category, ncol = 1, strip.position = "left", scales = "free_y") +
  geom_text(aes(x = 0, label = pathway), hjust = 0, family = "Helvetica") +
  theme(panel.background = element_blank(),
        text = element_text(family = "Helvetica"),
        legend.position = "None",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown(),
        axis.text.x = element_text(color = "black", size = 10),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_line(),
        strip.text = element_blank(),
        panel.spacing = unit(0, "mm"), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  ggtitle("KEGG Pathways")

