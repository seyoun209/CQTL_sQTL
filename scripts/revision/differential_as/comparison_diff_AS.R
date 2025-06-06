setwd("/work/users/s/e/seyoun/CQTL_sQTL/")
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(data.table)
library(httpgd)
library(tidyr)
library(tidyverse)
library(gridtext)
library(plotgardener)
library(ggtext)
library(ggrepel)
library(ggpmisc)
source("scripts/utils/utils.R")

#Load the data
# ratios_oa

pbs_oa <- load("output/clu_oa/PBSvsOA.Rdata")
for (name in pbs_oa) {
  # Construct the new name by prefixing with "fnf_"
  new_name <- paste0(name,"_oa")

  # Assign the object to the new name in the global environment
  assign(new_name, get(name))
  remove(name)
}

ratios_oa_qc <- read.table("output/clu_oa/ratio_oa.txt",sep='\t',header=T) |>
as.data.frame()
OA_deltapsiCalc_df <- calculate_delta_psi(ratios_oa_qc, "CTL", "OA")

OA_df <- OA_deltapsiCalc_df %>%
  mutate(junction_key = str_split_fixed(Junction, ":", 4)[,1:3] %>% 
           apply(1, paste, collapse = ":")) |>
            select(deltaPSI, junction_key) |> 
           mutate(group ="OA") 

#Ratios FNF

load("output/clu_fnf/introns_fnf_joinAll") #introns_fnf_pval_include

fnf_df <- introns_fnf_pval_include %>%
  filter(p.adjust < 0.05) %>%                     # keep only significant ones
  mutate(junction_key = str_split_fixed(phe_id, ":", 4)[,1:3] %>% 
           apply(1, paste, collapse = ":")) |> 
           select(deltapsi_batch, junction_key) |> 
           mutate(group ="FNF") |> 
           dplyr::rename(deltaPSI = "deltapsi_batch")

# Join the two datasets by junction_key
combined_deltaPSI <- inner_join(OA_df, fnf_df, by = "junction_key", suffix = c(".OA", ".FNF"))


#
combined_deltaPSI_color <- combined_deltaPSI %>%
  mutate(color_cat = case_when(
    deltaPSI.OA >= 0.15 & deltaPSI.FNF > 0 ~ "#37a7db",       # Quadrant 1: both high
    deltaPSI.OA <= -0.15 & deltaPSI.FNF < 0 ~ "#37a7db",     # Quadrant 3: both low
    (deltaPSI.OA > 0.15 & deltaPSI.FNF < 0) | (deltaPSI.OA < -0.15 & deltaPSI.FNF > 0) ~ "#7ecdbb",         # Quadrant 2: OA low, FN‑f high
    TRUE ~ "grey"                                               # Otherwise grey
  ))

fnf_df_allcol <- introns_fnf_pval_include %>%
  filter(p.adjust < 0.05) %>%                     # keep only significant ones
  mutate(junction_key = str_split_fixed(phe_id, ":", 4)[,1:3] %>% 
           apply(1, paste, collapse = ":")) |> 
           mutate(group ="FNF") |> 
           dplyr::rename(deltaPSI = "deltapsi_batch")


scatter_plot_PBS_OA <- ggplot(combined_deltaPSI_color, aes(x = deltaPSI.OA, y = deltaPSI.FNF)) +
  geom_hline(yintercept = 0, lty = 2, color = "grey25", linewidth = 0.25) +
  geom_vline(xintercept = c(0.15, -0.15), lty = 2, color = "grey25", linewidth = 0.25)+
  geom_point(aes(color = color_cat), size = 1) +
  scale_color_identity() +
  labs(x = "\u0394 PSI OA/PBS",
     y = "\u0394 PSI FNF/PBS")+
  stat_poly_eq(
    use_label(c( "R2", "P")),
    formula = y ~ x,
    parse = TRUE,
    size = 3, color = "black"
  ) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.2) +
  theme( strip.placement = "outside",
         axis.line.y = element_line(linewidth = 0.35),
         axis.line.x = element_line(linewidth = 0.35),
         axis.ticks.length.y = unit(-0.1, "cm"),
         axis.title.y = element_markdown(size = 8, family = "Helvetica",
                                         margin = margin(r = -15)),
         axis.title.x = element_markdown(size = 8, family = "Helvetica",
                                         margin = margin(t = +15)),
         text = element_text(family = "Helvetica"),
         axis.text.y = element_text(color = "black", size = 6),
         axis.text.x = element_text(color = "black", size = 6),
         strip.background = element_blank(),
         strip.text.x.bottom = element_markdown(size = 8, margin = margin(t = 1)),
         panel.background = element_rect(fill = "transparent", color = "transparent"),
         plot.background = element_rect(fill = "transparent", color = "transparent"),
         panel.grid = element_blank(),
         legend.position = "none")
# Save the plot
save(scatter_plot_PBS_OA, file = "output/revision/plots/scatter_plot_PBS_OA.rda")

save(combined_deltaPSI_color,OA_deltapsiCalc_df, introns_fnf_pval_include,
 file = "output/revision/data/combined_deltaPSI_color.Rdata")

load("output/revision/data/combined_deltaPSI_color.Rdata")


pdf(file = "output/revision/plots/pbs_fnf_scatterplot.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4)

pageCreate(width = 4, height =4 , default.units = "inches", 
showGuides = FALSE)

plotGG(scatter_plot_PBS_OA, x = 0.5, y = 0.5, width = 3.5, height = 3.5)

dev.off()


# Find the gene sets opposite

combined_deltaPSI_groups <- combined_deltaPSI %>%
  mutate(opposite = case_when(
    deltaPSI.OA >= 0.15 & deltaPSI.FNF > 0 ~"no",      
    deltaPSI.OA <= -0.15 & deltaPSI.FNF < 0 ~ "no",     
    (deltaPSI.OA >= 0.15 & deltaPSI.FNF < 0) |                 
    (deltaPSI.OA <= -0.15 & deltaPSI.FNF > 0) ~ "yes",        
    TRUE ~ "not_sig"                                               
  ))

junctions_key_opposite <- combined_deltaPSI_groups %>%
  filter(opposite == "yes") %>%
  select(junction_key)

opposite_ensg <- introns_fnf_pval_include %>%
  mutate(junction_key = str_split_fixed(phe_id, ":", 4)[,1:3] %>% 
           apply(1, paste, collapse = ":"),
         ensg = sub("\\..*$", "", ensemblID)) %>% 
  filter(junction_key %in% junctions_key_opposite$junction_key) %>% 
  select(ensg) %>% 
  distinct() 
write.table(opposite_ensg,
 file = "output/revision/data/opposite_genes.txt",sep='\t',quote=F,row.names=F,col.names=F)

#system("scripts/Differential_splicing/run_homer.sh /work/users/s/e/seyoun/CQTL_sQTL/output/revision/data/opposite_genes.txt /work/users/s/e/seyoun/CQTL_sQTL/output/clu_fnf/background_gene_set.txt /work/users/s/e/seyoun/CQTL_sQTL/output/revision/homer/homer_opposite_genes")



#-------------------------------------------------------------------------------------
# Pathway
reactome_data <- read_delim("output/revision/homer/homer_opposite_genes/reactome.txt") |>
  mutate(pval = exp(1)^logP) |>
  dplyr::filter(pval < 0.01) |>
  distinct(Term, .keep_all = TRUE) |>
  mutate(`-log10pval` = -log10(pval)) |>
  mutate(category = "Reactome Pathway")

kegg_data <- read_delim("output/revision/homer/homer_opposite_genes/kegg.txt") |>
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
write_csv(pathway_table, file = "output/revision/homer/table/pathway_table_sig.csv")

# Plot top 10 significant for each category
pathway_plotting <- pathway_table |>
  dplyr::filter(Term %in% c("Rho GTPase cycle",
    "Transport of inorganic cations/anions and amino acids/oligopeptides"
                            )) |>
  arrange(`-log10pval`)

pathway_plotting$Term <- factor(pathway_plotting$Term, levels = pathway_plotting$Term)

pathway_barplots <- ggplot(pathway_plotting, aes(x = `-log10pval`, y = Term)) +
  geom_vline(xintercept = 1, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 2, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 3, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 4, color = "grey75", alpha = 0.4) +
 # geom_vline(xintercept = 5, color = "grey75", alpha = 0.4) +
  geom_bar(stat = "identity", fill="#c1daf3") +
  scale_x_continuous(expand = c(0, 0), name = "-log~10~pval", limits = c(0, 4),
                     breaks = seq(0, 4, 1)) +
  #facet_wrap(~category, ncol = 1, strip.position = "left", scales = "free_y") +
  geom_text(aes(x = 0, label = Term), hjust = 0, family = "Helvetica",
            size = 2.5) +
  theme(panel.background = element_rect(fill = 'transparent', color = "transparent"),
        plot.background = element_rect(fill = 'transparent', color = "transparent"),
        text = element_text(family = "Helvetica"),
        legend.position = "None",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown(size = 6),
        axis.text.x = element_text(color = "black", size = 4),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_line(linewidth = 0.25),
        strip.text = element_blank(),
        panel.spacing = unit(0, "mm"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8)) +
  ggtitle("Pathways")
ggsave(filename = "output/revision/plot/pathway_barplot.pdf",
       plot = pathway_barplots, width = 5, height = 5, units = "in")
save(pathway_barplots, file = "output/revision/plot/pathway_barplots.rda")



go_data <- read_delim("output/revision/homer/homer_opposite_genes/biological_process.txt") |>
  mutate(pval = exp(1)^logP) |>
  dplyr::filter(pval < 0.01)
sig_go <- reduceGO(go_data,
                   category = "GO")

go_table <- sig_go |>
  dplyr::select(-`Entrez Gene IDs`, -pval, -logP, -size, -termUniqueness, -termUniquenessWithinCluster,
                -termDispensability, -category) |>
  relocate(`-log10pval`, .after = Enrichment) |>
  arrange(desc(`-log10pval`))

write_csv(go_table, file = "output/revision/homer/table/GO_sig.csv")


Siggo_plotting <- sig_go |>
  dplyr::filter(Term == parentTerm) |>
  dplyr::filter(parentTerm %in%  unique(go_table$parentTerm)[1:10]) |>
  arrange(`-log10pval`)

#------------------------------------------
#75 Gene effect size comparison

load("output/clu_fnf/introns_fnf_joinAll")
introns_fnf_sig <- introns_fnf_pval_include %>% dplyr::filter(p.adjust < 0.05)
introns_fnf_subset_sig <- introns_fnf_sig[abs(introns_fnf_sig$deltapsi_batch) > 0.15,]
introns_fnf_all_sig <- introns_fnf_sig %>% dplyr::filter(abs(deltapsi_batch) > 0.15)
load("./output/clu_oa/introns_oa_joinAll")
introns_oa_sig <- introns_oa_pval_include %>% dplyr::filter(p.adjust < 0.05)
introns_oa_subset_sig <- introns_oa_sig[abs(introns_oa_sig$deltapsi_batch) > 0.15,]


OA_df <- introns_oa_subset_sig %>%
  # Assuming 'Junction' contains "chr:start:end:..." and we need first 3 fields:
  mutate(junction_key = str_split_fixed(phe_id, ":", 4)[,1:3] %>% 
           apply(1, paste, collapse = ":")) %>%
  select(deltapsi,deltapsi_batch, junction_key, gene,ensemblID ) %>%   
  mutate(group = "OA")


fnf_df <- introns_fnf_subset_sig %>%
  # Build junction_key from the phe_id field (first 3 fields: chr, start, end)
  mutate(junction_key = str_split_fixed(phe_id, ":", 4)[,1:3] %>% 
           apply(1, paste, collapse = ":")) %>%
  select(deltapsi,deltapsi_batch, junction_key,gene,ensemblID ) %>% 
  mutate(group = "FNF")


combined_deltaPSI <- inner_join(OA_df, fnf_df, by = "junction_key", 
                                suffix = c(".OA", ".FNF"))

combined_deltaPSI_batch_color <- combined_deltaPSI %>%
  mutate(color_cat = case_when(
    deltapsi_batch.OA >= 0.15 & deltapsi_batch.FNF > 0  ~ "black",
    deltapsi_batch.OA <= -0.15 & deltapsi_batch.FNF < 0 ~ "black",
    (deltapsi_batch.OA >= 0.15 & deltapsi_batch.FNF < 0) | (deltapsi_batch.OA <= -0.15 & deltapsi_batch.FNF > 0) ~ "red",
    TRUE ~ "grey"
  ))



sGenes_sig_scatterPlot <- ggplot(combined_deltaPSI_batch_color, aes(x = deltapsi_batch.OA, y = deltapsi_batch.FNF)) +
  #geom_hline(yintercept = 0, lty = 2, color = "grey25", linewidth = 0.25) +
  #geom_vline(xintercept = 0, lty = 2, color = "grey25", linewidth = 0.25) +
  geom_point(aes(color = color_cat), size = 1) +
  scale_color_identity() +
labs(x = "\u0394 PSI OA/PBS",
     y = "\u0394 PSI FNF/PBS")+
  stat_poly_eq(
    use_label(c( "R2", "P")),
    formula = y ~ x,
    parse = TRUE,
    size = 3, color = "black"
  ) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.2) +
  theme(
    strip.placement = "outside",
    axis.line.y = element_line(linewidth = 0.35),
    axis.line.x = element_line(linewidth = 0.35),
    axis.ticks.length.y = unit(-0.1, "cm"),
    axis.title.y = element_markdown(size = 8, family = "Helvetica", margin = margin(r = -15)),
    axis.title.x = element_markdown(size = 8, family = "Helvetica", margin = margin(t = +15)),
    text = element_text(family = "Helvetica"),
    axis.text.y = element_text(color = "black", size = 6),
    axis.text.x = element_text(color = "black", size = 6),
    strip.background = element_blank(),
    strip.text.x.bottom = element_markdown(size = 8, margin = margin(t = 1)),
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    panel.grid = element_blank(),
    legend.position = "none"
  )


save(sGenes_sig_scatterPlot, file="output/revision/plots/scatterplot_sig_fnf_oa.rda")

pdf(file = "output/revision/plots/scatterplot_sig_fnf_oa.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4)

pageCreate(width = 4, height =4 , default.units = "inches", showGuides = FALSE)
plotGG(sGenes_sig_scatterPlot, x = 0.5, y = 0.5, width = 3.5, height = 3.5)

dev.off()