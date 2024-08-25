setwd("/work/users/s/e/seyoun/CQTL_sQTL")
library(plotgardener)
library(ggplot2)
library(RColorBrewer)
library(yaml)
library(ggvenn)
library(scales)
library(colorspace)
library(ggtext)
library(ggforce)
library(AnnotationDbi)
source("scripts/sQTL_rscripts/utils.R")

#Figure 3A -Venndiagrma --------------------------------------------------------

## open up data ####
response_pbs_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_pbs_resInv.rds")
response_fnf_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_fnf_resInv.rds")

# Subset significant reponse sQTL only
pbs_resQTL <- response_pbs_results %>%
  dplyr::filter(interaction_pval  < 0.05) |>
  arrange(interaction_pval)

fnf_resQTL <- response_fnf_results %>%
  dplyr::filter(interaction_pval  < 0.05) |>
  arrange(interaction_pval)

# Subset high confidence re-sQTls ####

pbs_highConf_resQtL <- response_pbs_results %>% 
  dplyr::filter(interaction_pval  < 0.05) %>%
  dplyr::filter(abs(delta_beta ) > 0.2) |>
  dplyr::filter(minor_alle_count >= 5) |>
  arrange(interaction_pval)

fnf_highConf_resQtL <- response_fnf_results %>% 
  dplyr::filter(interaction_pval  < 0.05) %>%
  dplyr::filter(abs(delta_beta ) > 0.2) |>
  dplyr::filter(minor_alle_count >= 5) |>
  arrange(interaction_pval)

pbs_specific_resQTL <- pbs_highConf_resQtL |> 
  dplyr::filter(!clusterID %in% fnf_highConf_resQtL$clusterID ) 

fnf_specific_resQTL <- fnf_highConf_resQtL |> 
  dplyr::filter(!clusterID %in% pbs_highConf_resQtL$clusterID )  

# pieChart  fro the significant response sQTL   

resqtl_tibble <- tibble(values = unique(c(pbs_resQTL$ensg,fnf_resQTL$ensg))) %>%
  mutate(PBS = values %in% pbs_resQTL$ensg,
         FNF = values %in% fnf_resQTL$ensg)

# Calculate the counts (adjust these based on your data)
total_PBS <- sum(resqtl_tibble$PBS)
total_FNF <- sum(resqtl_tibble$FNF)
overlap <- sum(resqtl_tibble$PBS & resqtl_tibble$FNF)
only_PBS <- total_PBS - overlap
only_FNF <- total_FNF - overlap

# Calculate the maximum count for scaling
max_count <- max(total_PBS, total_FNF)

resQTL_venn <- ggplot() +
  # PBS circle
  geom_circle(aes(x0 = 0, y0 = 0.6, r = sqrt(total_PBS/max_count)), 
              fill = "#BFDDFF", color = NA, alpha = 0.5) +
  # FNF circle
  geom_circle(aes(x0 = 0, y0 = -0.6, r = sqrt(total_FNF/max_count)), 
              fill = "#FFDDA2", color = NA, alpha = 0.5) +
  # Labels
  geom_text(aes(x = 0, y = 0.7, label = only_PBS), size = 3) +
  geom_text(aes(x = 0, y = -0.7, label = only_FNF), size = 3) +
  geom_text(aes(x = 0, y = -0.075, label = overlap), size = 3) +
  # Adjust the plot
  coord_fixed() +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", color = "transparent"),
        plot.margin = margin(t = 5, r = 10, b = 5, l = 5))

ggsave(filename = "output/results_plots/sqtl_plots/venn_re_sGenes.pdf", 
       plot =resQTL_venn,width = 4, height = 4, units = "in")
save(resQTL_venn, file = "output/results_plots/sqtl_plots/venn_resGenes.rda")

#Make the pie chart  for the high confidence------------------------------------

#PBS
# Calculate the counts
total_resQTL <- pbs_resQTL |> distinct(ensg) |> nrow()
highConf_resQTL <- pbs_specific_resQTL |> distinct(ensg) |> nrow()
other_resQTL <- total_resQTL - highConf_resQTL

# Create a data frame for the pie chart
pie_data <- data.frame(
  category = c("High Confidence", "Other"),
  count = c(highConf_resQTL, other_resQTL)
)


pie_data <- pie_data %>%
  mutate(label = case_when(
    category == "Other" ~ sprintf("\n%d", count),
    TRUE ~ sprintf("%s\n%s\n%d", "pbs-specific","high confidence" , count)
  ))


# Create the pie chart
highConf_pbs_pie <- ggplot(pie_data, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start =4.8) +
  scale_fill_manual(values = c("High Confidence" = lighten("#2057A7",amount = 0.3),"Other" =lighten("#BFDDFF",amount = 0.3))) +
  #geom_text(aes(label = label), position = position_stack(vjust = 0.5),size=3) +
  labs(subtitle = paste("Total PBS re-sGenes:", total_resQTL)) +
  theme_void() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background =element_blank(),
        legend.position = "none",
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        plot.margin = margin(t = 5, r = 10, b = 5, l = 5))

save(highConf_pbs_pie, file = "output/results_plots/sqtl_plots/highConf_pbs_pie.rda")
# fnf

# Calculate the counts
total_resQTL <- fnf_resQTL |> distinct(ensg) |> nrow()
highConf_resQTL <- fnf_specific_resQTL |> distinct(ensg) |> nrow()
other_resQTL <- total_resQTL - highConf_resQTL

# Create a data frame for the pie chart
pie_data <- data.frame(
  category = c("High Confidence", "Other"),
  count = c(highConf_resQTL, other_resQTL)
)

# Calculate percentages and create labels with both percentage and count
# pie_data <- pie_data %>%
#   mutate(percentage = count / sum(count) * 100,
#          label = case_when(
#            category == "Other" ~ sprintf("\nn=%d\n%.1f%%", count, percentage),
#            TRUE ~ sprintf("%s\nn=%d\n%.1f%%", category, count, percentage)
#          ))

pie_data <- pie_data %>%
  mutate(label = case_when(
    category == "Other" ~ sprintf("\n%d", count),
    TRUE ~ sprintf("%s\n%s\n%d", "FN-f-specific","high confidence" , count)
  ))



# Create the pie chart
highConf_fnf_pie <-ggplot(pie_data, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start =1) +
  scale_fill_manual(values = c("High Confidence" = lighten("#F2BC40",amount = 0.3), "Other" = lighten("#FFDDA2",0.3))) +
  #geom_text(aes(label = label), position = position_stack(vjust = 0.5),size=3) +
  labs(subtitle = paste("Total FN-f re-sGenes:", total_resQTL)) +
  theme_void() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background =element_blank(),
        legend.position = "none",
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        plot.margin = margin(t = 5, r = 10, b = 5, l = 5))

save(highConf_fnf_pie, file = "output/results_plots/sqtl_plots/highConf_fnf_pie.rda")


# Figure 3b for the distance plot  and the boxplots 
#Calculate distance  of snp and the intron junction either intron junction start or the intron junction end. 

# Calculate distances to nearest intron junction
resqtl_fnf_distance <- response_fnf_results %>%
  mutate(
    distance_to_start = abs(var_from - phe_from),
    distance_to_end = abs(var_from - phe_to),
    distance_to_nearest_junction = pmin(distance_to_start, distance_to_end),
    junction_side = ifelse(distance_to_start <= distance_to_end, "Start", "End")
  )

resqtl_pbs_distance <- response_pbs_results %>%
  mutate(
    distance_to_start = abs(var_from - phe_from),
    distance_to_end = abs(var_from - phe_to),
    distance_to_nearest_junction = pmin(distance_to_start, distance_to_end),
    junction_side = ifelse(distance_to_start <= distance_to_end, "Start", "End")
  )




pbs_qtl_distance <- resqtl_pbs_distance |> dplyr::select(distance_to_nearest_junction) |>
  mutate(group = "PBS") 

fnf_qtl_distance <- resqtl_fnf_distance |> dplyr::select(distance_to_nearest_junction) |>
  mutate(group = "FNF")

all_distance_qtl <- bind_rows(pbs_qtl_distance,fnf_qtl_distance) |> mutate(group= factor(group ,levels = c("PBS","FNF")))

summarize_sig_qtl_dis <- all_distance_qtl %>%
  group_by(group) %>%
  summarise(
    mean = mean(distance_to_nearest_junction), 
    median = median(distance_to_nearest_junction),
    total_snps = n(),
    snps_within_5kb = sum(distance_to_nearest_junction <= 5000),
    percent_within_5kb = (snps_within_5kb / total_snps) * 100
  )

# Calculate y-coordinates for mean points
summarize_dis_qtl_fi <- summarize_sig_qtl_dis %>%
  group_by(group) %>%
  mutate(y_at_mean = n() / (max(all_distance_qtl$distance_to_nearest_junction) / binwidth))
#line_plot
binwidth <- 5000

# Main plot
distance_linePlot <- ggplot(all_distance_qtl, aes(x = distance_to_nearest_junction, group=group, colour=group)) +
  stat_bin(geom = "line", binwidth = binwidth, size = 0.5) +
  scale_x_continuous(limits = c(0, 100000), name = "Distance to intron splice junction (kb)", expand = c(0,0),
                     labels = label_number(scale = .001, big.mark = "")) +
  scale_y_continuous(name = "Number of sSNPs", expand = c(0,0),
                     limits = c(0, NA)) +
  scale_color_manual(values = c("PBS" = "#0067B9", "FNF" = "#FCCE52"), labels = c("PBS" = "PBS", "FNF" = "FN-f")) +
  geom_vline(xintercept = 5000, linetype = "dashed", linewidth = 0.5) +
  # geom_text(data = subset(summarize_dis_qtl_fi, group == "PBS"), 
  #           aes(x = 10000, y = 2700, 
  #               label = sprintf("%s: %.1f%% of lead sQTLs are within \n5 Kb from the splice junction", 
  #                               group, percent_within_5kb)),
  #           hjust = 0, vjust = 2.5, size = 3, color = "#0067B9") +
  # geom_text(data = subset(summarize_dis_qtl_fi, group == "FNF"), 
  #           aes(x = 10000, y = 2400, 
  #               label = sprintf("%s: %.1f%% of lead sQTLs are within \n5 Kb from the splice junction", 
  #                               group, percent_within_5kb)),
  #           hjust = 0, vjust = 4, size = 3, color = "#FCCE52") +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica"),
    axis.line = element_line(linewidth = 0.25),
    axis.title.x = element_markdown(size = 6, margin = margin(r = 10)),
    axis.title.y = element_markdown(size = 6, margin = margin(r = 5)),
    axis.text = element_text(color = "black", size = 6),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.position = "none",
    plot.title =element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )
#facet_zoom(xlim = c(0, 20000), horizontal = FALSE, show.area = TRUE) +
#theme(strip.background = element_rect(fill = "grey85", colour = "grey85", size = 0.5)) 

save(distance_linePlot, file = "output/results_plots/sqtl_plots/distance_linePlot.rda")


# Assuming you have already created resqtl_fnf_distance and resqtl_pbs_distance

# Prepare the data
pbs_qtl_distance <- resqtl_pbs_distance |> 
  dplyr::select(distance_to_nearest_junction) |>
  mutate(group = "PBS") 
fnf_qtl_distance <- resqtl_fnf_distance |> 
  dplyr::select(distance_to_nearest_junction) |>
  mutate(group = "FNF")
all_distance_qtl <- bind_rows(pbs_qtl_distance, fnf_qtl_distance) |> 
  mutate(group = factor(group, levels = c("PBS", "FNF")))

binwidth <- 5000  # 5kb bins
max_distance <- 105000  # Set maximum distance to 100kb

all_distance_qtl_percentages <- all_distance_qtl %>%
  group_by(group) %>%
  mutate(bin = cut(pmin(distance_to_nearest_junction, max_distance), 
                   breaks = seq(0, max_distance, by = binwidth),
                   include.lowest = TRUE,
                   right = FALSE,
                   labels = FALSE)) %>%
  group_by(group, bin) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(percentage = count / sum(count) * 100,
         bin_start = (bin - 1) * binwidth)

# Main plot
distance_percentage_plot <- ggplot(all_distance_qtl_percentages, aes(x = bin_start, y = percentage, group = group, colour = group)) +
  geom_line(size = 0.5) +
  scale_x_continuous(
    limits = c(0, 100000),  # Set to exactly 100000
    name = "Distance to intron splice junction (kb)", 
    expand = c(0,0),
    labels = label_number(scale = .001, big.mark = "")
  ) +
  scale_y_continuous(
    name = "Percentage of sSNPs (%)", 
    expand = c(0,0),
    limits = c(0, 50),  # Set y-axis from 0 to 50%
    labels = function(x) paste0(x, "%")
  ) +
  scale_color_manual(values = c("PBS" = "#0067B9", "FNF" = "#FCCE52"), 
                     labels = c("PBS" = "PBS", "FNF" = "FN-f")) +
  geom_vline(xintercept = 5000, linetype = "dashed", linewidth = 0.5) +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica"),
    axis.line = element_line(linewidth = 0.25),
    axis.title.x = element_markdown(size = 6, margin = margin(t = 10)),
    axis.title.y = element_markdown(size = 6, margin = margin(r = 5)),
    axis.text = element_text(color = "black", size = 6),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.position = "none",
    plot.title = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )
save(distance_percentage_plot, file = "output/results_plots/sqtl_plots/distance_percentage_plot.rda")


#Density plot for sQTLs

density_distancePlot <- ggplot(all_distance_qtl, aes(x = distance_to_nearest_junction, fill = group, colour = group)) +
  geom_density(alpha = 0.5, na.rm = TRUE, color = NA) +
  labs(y = "Density", x = "Distance to intron splice junction (kb)") +
  #facet_wrap(~ group, scales = "free", ncol = 2) +
  scale_x_continuous(limits = c(0, 100000), name = "Distance to intron splice junction (kb)", expand = c(0,0),
                     labels = label_number(scale = .001, big.mark = ""))+
  scale_y_continuous(name = "Density", expand = c(0,0),
                     limits = c(0, NA)) +
  geom_vline(data = summarize_dis_qtl_fi, aes(xintercept = mean, color = group), 
             linetype = "dashed", size = 0.5) + 
  geom_text(data = summarize_dis_qtl_fi, 
            aes(x = -Inf, y = Inf, color = group, 
                label = sprintf("%s mean: %.2f kb", group, mean/1000)), 
            vjust = 2, hjust = -0.5, size = 3)+
  # geom_text(data = summarize_dis_qtl_fi, 
  #           aes(x = mean, y = Inf, color = group, 
  #               label = sprintf("Mean = %.2f kb", mean/1000)), 
  #           vjust = 5, hjust = -0.2, size = 3) +
  scale_fill_manual(values = c("PBS" = "#0067B9", "FNF" = "#FCCE52"), labels = c("PBS" = "PBS", "FNF" = "FN-f")) +
  scale_color_manual(values = c("PBS" = "#0067B9", "FNF" = "#FCCE52"), labels = c("PBS" = "PBS", "FNF" = "FN-f")) +
  theme_minimal()+
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
        axis.title.x = element_markdown(size = 6, family = "Helvetica", margin = margin(r = -15)),
        axis.ticks = element_blank(),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_markdown(size = 6, family = "Helvetica", margin = margin(r = -15)),
        text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.text.x =  element_text(color = "black", size = 6),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        #panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, family = "Helvetica", size = 10, margin = margin(b = 3)))

save(density_distancePlot, file = "output/results_plots/sqtl_plots/density_distancePlot.rda")


# Figure 3c Boxplot of where sQTLS are inside of the gene or outside of the gene. 


txdb <- loadDb("../crispr/02.test_seq/gencode.v45.annotation.TxDb")
txdb_genes <- genes(txdb)

remove_version <- function(ensg) {
  gsub("\\.\\d+$", "", ensg)
}

# Prepare txdb_genes data
txdb_genes_df <- as.data.frame(txdb_genes) %>%
  mutate(gene_id = remove_version(gene_id)) %>%
  dplyr::select(gene_id, seqnames, start, end)

# Function to check if var_id is within gene boundaries
check_within_gene <- function(df, gene_info) {
  df %>%
    left_join(gene_info, by = c("ensg" = "gene_id")) %>%
    mutate(within_gene = var_chr == seqnames & var_from >= start & var_from <= end) %>%
    dplyr::select(-seqnames, -start, -end)
}

# Apply the function to both datasets
response_pbs_results_checked <- check_within_gene(response_pbs_results, txdb_genes_df)
response_fnf_results_checked <- check_within_gene(response_fnf_results, txdb_genes_df)

# Function to calculate percentages
calculate_percentages <- function(df) {
  df %>%
    filter(!is.na(within_gene)) %>%  # Remove NA values
    group_by(within_gene) %>%
    summarise(count = n()) %>%
    mutate(
      total = sum(count),
      percentage = count / total * 100,
      category = if_else(within_gene, "Within Gene", "Outside Gene")
    ) %>%
    dplyr::select(category, count, percentage)
}

# Calculate percentages for both datasets
pbs_percentages <- calculate_percentages(response_pbs_results_checked)
fnf_percentages <- calculate_percentages(response_fnf_results_checked)

# Combine the results
combined_percentages <- bind_rows(
  mutate(pbs_percentages, dataset = "PBS"),
  mutate(fnf_percentages, dataset = "FNF")
)

combined_percentages <- combined_percentages |> mutate(dataset = factor(dataset, levels=c("PBS","FNF")))

# Create the bar plot
sqtl_locBarplot <- ggplot(combined_percentages, aes(x = category, y = percentage, fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 2) +
  scale_fill_manual(values = c("PBS" = "#BFDDFF", "FNF" = "#FFDDA2"), labels = c("PBS" = "PBS", "FNF" = "FN-f")) +
  labs(x = "", y = "% of sSNPs", fill = "Dataset") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(combined_percentages$percentage) * 1.1)) +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica"),
    axis.line = element_line(linewidth = 0.25),
    axis.title.x = element_blank(),
    axis.title.y = element_markdown(size = 6, margin = margin(r = 5)),
    axis.text = element_text(color = "black", size = 6),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.position = c(0.25,0.9),
    legend.title = element_blank(),
    legend.text = element_markdown(size = 6, margin = margin(r = 10)),
    legend.key.size = unit(4, 'mm'),
    plot.title =element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )
save(sqtl_locBarplot, file = "output/results_plots/sqtl_plots/sqtl_locBarplot.rda")

# Make box plot for the genes

plot_slc26a4 <- create_boxplot("SLC26A4", pbs_highConf_resQtL)
plot_MAPK8 <- create_boxplot("MAPK8", fnf_highConf_resQtL)
plot_RHBDF2 <- create_boxplot("RHBDF2", pbs_highConf_resQtL)



save(plot_slc26a4, file="output/results_plots/sqtl_plots/plot_slc26a4_boxplot.rda")
save(plot_MAPK8, file="output/results_plots/sqtl_plots/plot_MAPK8_boxplot.rda")
save(plot_RHBDF2, file="output/results_plots/sqtl_plots/plot_RHBDF2_boxplot.rda")


#Make the locus zoom plot for the  selected genes


#plottin------------------------------------------------------------------------

pdf(file = "output/results_plots/sqtl_plots/main_figure3.pdf",   # The directory you want to save the file in
    width = 9.55, # The width of the plot in inches
    height = 9)
pageCreate(width = 9.55, height = 9, showGuides = FALSE)
plotText("a", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

load("output/results_plots/sqtl_plots/venn_resGenes.rda")
load("output/results_plots/sqtl_plots/highConf_pbs_pie.rda")
load("output/results_plots/sqtl_plots/highConf_fnf_pie.rda")
plotText("Significant genetic and condition interaction effect",
         just = c("center", "top"), fontfamily = "Helvetica",fontsize = 8,x = 2, y = 0.35)
plotGG(resQTL_venn, x = 0.63, y = 0.5, width = 1.75, height = 2.5)

plotText("PBS \nre-sGenes",
         just = c("center","top"), fontfamily = "Helvetica",fontsize = 8,x = 0.5, y = 1.3, fontcolor = "#2057A7",fontface = "bold" )

plotText("FN-f \nre-sGenes",
         just =c("center","top"), fontfamily = "Helvetica",fontsize = 8,x = 0.5, y = 2.25, fontcolor ="#F2BC40",fontface = "bold" )

plotSegments(
  x0 = c(1.5,2.5), y0 = c(0.7,0.55), x1 = 1.5, y1 = 0.55,
  default.units = "inches",
  lwd = 1, lty = 2
)

plotSegments(
  x0 = c(1.5,2.5), y0 = c(2.75,2.9), x1 = 1.5, y1 = 2.9,
  default.units = "inches",
  lwd = 1, lty = 2
)

plotGG(highConf_pbs_pie, x = 2.5, y =0.27, width = 1.2, height = 1.2)
plotGG(highConf_fnf_pie, x = 2.5, y = 2, width = 1.2, height = 1.2)


plotText("477 \nhigh confidence \nPBS-specific",
         just = c("center","top"), fontfamily = "Helvetica",fontsize = 7,x = 3, y = 1.3, fontcolor = "#2057A7",fontface="bold")

plotText("200 \nnhigh confidence \nFN-f-specific",
         just =c("center","top"), fontfamily = "Helvetica",fontsize = 7,x = 3, y = 1.8, fontcolor ="#F2BC40",fontface="bold" )


plotText("557",
         just = c("center","top"), fontfamily = "Helvetica",fontsize = 7,x = 3.1, y = 0.65)

plotText("401",
         just =c("center","top"), fontfamily = "Helvetica",fontsize = 7,x = 3.1, y = 2.7)


# Distance plot 

plotText("b", x = 3.65, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

load("output/results_plots/sqtl_plots/distance_percentage_plot.rda")
plotGG(distance_percentage_plot, x = 3.7, y = 0.2, width = 3.7, height = 2.9)

# plotCircle(
#   x = 4.4, y = 2.0, r = 0.035, fill = "#2057A7",linecolor="#2057A7",
#   default.units = "inches"
# )
# 
# plotCircle(
#   x = 4.4, y = 2.1, r = 0.035, fill = "#F2BC40" ,linecolor="#F2BC40",
#   default.units = "inches"
# )

plotText("5kb",
         just = c("center","top"), fontfamily = "Helvetica",fontsize = 8,x = 4.4, y = 0.28)

plotSegments(
  x0 = 4.4, y0 = 2.0, x1 = 4.45, y1 = 2.0,
  default.units = "inches",
  lwd = 1, lty = 2,linecolor="#2057A7"
)


plotSegments(
  x0 = 4.45, y0 = 1.5, x1 = 4.45, y1 = 2.0,
  default.units = "inches",
  lwd = 1, lty = 2,linecolor="#2057A7"
)


plotSegments(
  x0 = 4.45, y0 = 1.5, x1 = 4.6, y1 = 1.5,
  default.units = "inches",
  lwd = 1, lty = 2,linecolor="#2057A7"
)


plotSegments(
  x0 = 4.4, y0 = 2.0, x1 = 4.5, y1 = 2.0,
  default.units = "inches",
  lwd = 1, lty = 2,linecolor="#F2BC40"
)

plotSegments(
  x0 = 4.5, y0 = 1.8, x1 = 4.5, y1 = 2.0,
  default.units = "inches",
  lwd = 1, lty = 2,linecolor="#F2BC40"
)

plotSegments(
  x0 = 4.5, y0 = 1.8, x1 = 4.6, y1 = 1.8,
  default.units = "inches",
  lwd = 1, lty = 2,linecolor="#F2BC40"
)

plotText("PBS: 43.3% of lead sSNPs are within \n5kb from the splice junction",
         just = c("left","top"), fontfamily = "Helvetica",fontsize = 7,x = 4.6, y = 1.4, fontcolor = "#2057A7")

plotText("FN-f: 44.8% of lead sSNPs are within \n5kb from the splice junction",
         just =c("left","top"), fontfamily = "Helvetica",fontsize = 7,x = 4.6, y = 1.75, fontcolor ="#F2BC40" )


# boxplot for C-----------------------------------------------------------------

plotText("C", x = 7.3, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
load("output/results_plots/sqtl_plots/sqtl_locBarplot.rda")
plotGG(sqtl_locBarplot, x = 7.4, y = 0.5, width = 2.2, height = 2.5)


# Locus zoom plot D-----------------------------------------------------------------

plotText("D", x = 0.1, y = 3.2, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

load("output/results_plots/sqtl_plots/plot_slc26a4_boxplot.rda")
load("output/results_plots/sqtl_plots/plot_MAPK8_boxplot.rda")
load("output/results_plots/sqtl_plots/plot_RHBDF2_boxplot.rda")
# shared gene  RHBDF2


plot_sQTL_manhattan("RHBDF2", "PBS", x_start = 0.6, y_start = 3.6, width = 2.2, height = 2,zoom_range = 30000)
plotGG(plot_RHBDF2,x = 0.3, y = 6.8, width = 2.8, height = 1.8)
plotText("chr17:76479941-76481375", x= 1.7,y=6.7, fontfamily = "Helvetica", just = "center", fontsize = 7)



# PBS specific geneSLC26A4

plotText("E", x = 3.25, y = 3.2, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")


plot_sQTL_manhattan("SLC26A4", "PBS", x_start = 3.75, y_start = 3.6, width = 2.2, height = 2,zoom_range = 100000)
plotGG(plot_slc26a4,x = 3.45, y = 6.8, width = 2.8, height = 1.8)
plotText("chr7:107694480-107694621", x= 4.85,y=6.7, fontfamily = "Helvetica", just = "center", fontsize = 7)

# FNF specific gene MAPK8

plotText("F", x = 6.4, y = 3.2, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

plot_sQTL_manhattan("MAPK8", "FNF", x_start = 6.9, y_start = 3.6, width = 2.2, height = 2,zoom_range = 200000)
plotGG(plot_MAPK8,x = 6.6, y = 6.8, width = 2.8, height = 1.8)
plotText("chr10:48363284-48401612", x=8,y= 6.7, fontfamily = "Helvetica", just = "center", fontsize = 7)

dev.off()



# Side works to get the gene PSI 

#SLC26A4
test_boxplotInfo <- pbs_highConf_resQtL |> dplyr::filter(SYMBOL == "SLC26A4") |> slice_min(PBS_p)
psi_matrix <- psi_ratio[rownames(psi_ratio) %in% test_boxplotInfo$phe_id, ]
psi_data <- data.frame(sampleID = names(psi_matrix), psi = as.double(psi_matrix))

# Process variant data
variantID <- test_boxplotInfo$var_id
minorAllele <- unique(test_boxplotInfo$minor_allele)
rsID <- unique(test_boxplotInfo$rsID)
alleles <- unlist(strsplit(variantID, ":"))
ref_allele <- alleles[3]
alt_allele <- alleles[4]
protective_allele <- ifelse(minorAllele == ref_allele, alt_allele, ref_allele)

# Process genotype data
geno_data <- all_geno_transpose_df[rownames(all_geno_transpose_df) %in% variantID, ] %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
colnames(geno_data)[2] <- "genotype"

# Combine metadata
meta_data_fi <- geno_data %>%
  left_join(meta_catl_simple, by = c("sampleID" = "ID")) %>%
  left_join(psi_data, by = c("sampleID" = "sampleID"))

# Prepare final dataset
meta_combined_all <- meta_data_fi %>%
  mutate(
    genotype = factor(genotype, levels = c("0", "1", "2")),
    Condition = ifelse(Condition == "CTL", "PBS", ifelse(Condition == "FNF", "FN-f", NA)),
    Condition = factor(Condition, levels = c("PBS", "FN-f")),
    genotype = case_when(
      genotype == "2" ~ paste(minorAllele, minorAllele, sep = "/"),
      genotype == "1" ~ paste(protective_allele, minorAllele, sep = "/"),
      genotype == "0" ~ paste(protective_allele, protective_allele, sep = "/")
    )
  ) %>%
  mutate(
    genotype = factor(genotype, levels = c(
      paste(protective_allele, protective_allele, sep = "/"),
      paste(protective_allele, minorAllele, sep = "/"),
      paste(minorAllele, minorAllele, sep = "/")
    ))
  )


summary_stats <- meta_combined_all %>%
  group_by(Condition, genotype) %>%
  summarise(
    mean_psi = mean(psi, na.rm = TRUE),
    sd_psi = sd(psi, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  )


# Filter psi_ratio for all rows matching the cluster ID
all_cluster_psi <- ratios_fnf[str_detect(rownames(ratios_fnf), test_boxplotInfo$clusterID), ]


# Create a long format dataframe of all PSI data for this cluster
psi_data_long <- as.data.frame(all_cluster_psi) %>%
  tibble::rownames_to_column("phe_id") %>%
  tidyr::pivot_longer(cols = -phe_id, names_to = "sampleID", values_to = "psi")


# Combine with existing meta_data
meta_combined_all_expanded <- meta_combined_all %>%
  dplyr::select(-psi) %>%  # Remove the original psi column
  left_join(psi_data_long, by = "sampleID")


# Calculate summary statistics
summary_stats <- meta_combined_all_expanded %>%
  group_by(Condition, genotype, phe_id) %>%
  summarise(
    mean_psi = mean(psi, na.rm = TRUE),
    sd_psi = sd(psi, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  )

print(summary_stats)

