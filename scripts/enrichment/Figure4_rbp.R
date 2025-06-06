setwd("/work/users/s/e/seyoun/CQTL_sQTL/")
source("./scripts/sQTL_rscripts/utils.R")
library(plotgardener)
library(dplyr)
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(ggrepel)
library(stringr)
library(httpgd)

#encode bigwig name needs to be change or find the right region of it. 

yl_gn_bu <- brewer.pal(n = 9, name = "YlGnBu")
"#FFFFD9" "#EDF8B1" "#C7E9B4" "#7FCDBB" "#41B6C4" "#1D91C0" "#225EA8" "#253494" "#081D58"

#Figure 4 RBP

# dot plot for the odd ratio





# ggplot(fnf_rbp_subset, aes(x = odd_med, y = -log10(epval))) +
#   geom_point(aes(color = (odd_med > 1 & epval < 0.05)), size = 2) +
#   geom_text_repel(
#     data = subset(fnf_rbp_subset, odd_med > 1 & epval < 0.05),
#     aes(label = rbp),
#     direction ="x",
#     #nudge_x = 0.15,
#     vjust=0.05,
#     size = 2, color =  "darkorange") +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey70") +
#   geom_vline(xintercept = 1, linetype = "dashed", color = "grey70") +
#   scale_color_manual(values = c("FALSE" = "grey70", "TRUE" =  "darkorange")) +
#   coord_cartesian(xlim = c(0, 2.5), ylim = c(0, 4)) +
#   labs(x = "log(odds ratio)", y = "-log10(p-value)", color = "Significant") +
#   theme_minimal() +
#   theme(
#     panel.grid = element_blank(),
#     strip.placement = "outside",
#     panel.background = element_rect(fill = 'transparent', color = "transparent"),
#     plot.background = element_rect(fill = 'transparent', color = "transparent"),
#     text = element_text(family = "Helvetica"),
#     axis.text.y = element_text(color = "black", size = 8),
#     axis.title.y = element_text(size = 8),
#     axis.title.x = element_text(size = 8),
#     axis.text.x = element_text(color = "black", size = 4),
#     strip.text = element_text(face = "bold", hjust = 0),
#     axis.line.x = element_line(linewidth = 0.25),
#     axis.line.y = element_line(linewidth = 0.25),
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
#     strip.background = element_blank(),
#     axis.ticks = element_blank(),
#     legend.position = "none"
#   )
# 

# I will subset the too big  odd ratios. 

load("external_data/encode_rbp/rbp_prep_rbpOnly/pbs_rbp_elip.Rdata") #pbs_rbp_eclip 
load("external_data/encode_rbp/rbp_prep_rbpOnly/fnf_rbp_elip.Rdata") #fnf_rbp_eclip

pbs_rbp_eclip$ratio_dnv_med <- pbs_rbp_eclip$odd_dnv / pbs_rbp_eclip$odd_med
fnf_rbp_eclip$ratio_dnv_med <- fnf_rbp_eclip$odd_dnv / fnf_rbp_eclip$odd_med

pbs_rbp_eclip$ratio_upv_med <- pbs_rbp_eclip$odd_upv / pbs_rbp_eclip$odd_med
fnf_rbp_eclip$ratio_upv_med <- fnf_rbp_eclip$odd_upv / fnf_rbp_eclip$odd_med

# Set thresholds (you can adjust these based on your data and requirements)
threshold_dnv <- 1.5  # For example, odd_dnv shouldn't be more than 2.5 times odd_med
threshold_upv <- 1.5  # For example, odd_upv shouldn't be more than 1.5 times odd_med

pbs_rbp_eclip_filtered <- pbs_rbp_eclip[which(pbs_rbp_eclip$ratio_dnv_med <= threshold_dnv & pbs_rbp_eclip$ratio_upv_med <= threshold_upv), ]
fnf_rbp_eclip_filtered <- fnf_rbp_eclip[which(fnf_rbp_eclip$ratio_dnv_med <= threshold_dnv & fnf_rbp_eclip$ratio_upv_med <= threshold_upv), ]



pbs_rbp_subset <- pbs_rbp_eclip_filtered |>  group_by(rbp) %>%
  dplyr::filter(observed > 5) |>
  ungroup()  |>
  arrange(desc(odd_med))

fnf_rbp_subset <- fnf_rbp_eclip_filtered |>  group_by(rbp) %>%
  dplyr::filter(observed > 5) |>
  ungroup()  |>
  arrange(desc(odd_med))

ggplot(pbs_rbp_subset, aes(x = odd_med, y = -log10(epval))) +
  geom_point(aes(color = (odd_med > 1 & epval < 0.05)), size = 2) +
  geom_text_repel(
    data = subset(pbs_rbp_subset, odd_med > 1 & epval < 0.05),
    aes(label = rbp),
    direction = "y",
    nudge_x = 0.06,
    hjust = 0,
    segment.size = 0.15,
    segment.color = "#005587",
    size = 2.5,
    color = "#005587",
    box.padding = 0.25,
    point.padding = 0.1,
    force = 1,
    max.overlaps = Inf
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey70") +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#005587")) +
  coord_cartesian(xlim = c(0, 2.5), ylim = c(0, 4)) +
  labs(x = "log(odds ratio)", y = "-log10(p-value)", color = "Significant") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    strip.placement = "outside",
    panel.background = element_rect(fill = 'transparent', color = "transparent"),
    plot.background = element_rect(fill = 'transparent', color = "transparent"),
    text = element_text(family = "Helvetica"),
    axis.text.y = element_text(color = "black", size = 8),
    axis.title.y = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    axis.text.x = element_text(color = "black", size = 4),
    strip.text = element_text(face = "bold", hjust = 0),
    axis.line.x = element_line(linewidth = 0.25),
    axis.line.y = element_line(linewidth = 0.25),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
    strip.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

edited_boxplot_plot <- function(data, x, y, width, height, default.units = "inches", just = c("left", "top"), 
                                significant_color = "#005587", non_significant_color = "grey70",segment_length_factor = 0.5, label_spread_factor = 0.8) {
  
  # Create the base plot
  p <- ggplot(data, aes(x = odd_med, y = -log10(epval))) +
    geom_point(aes(color = (odd_med > 1 & epval < 0.05)), size = 1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey70") +
    scale_color_manual(values = c("FALSE" = non_significant_color, "TRUE" = significant_color)) +
    coord_cartesian(xlim = c(0, 2.2), ylim = c(0, 3)) +
    labs(x = "Odd ratio", y = "-log10(p-value)", color = "Significant") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      strip.placement = "outside",
      panel.background = element_rect(fill = 'transparent', color = "transparent"),
      plot.background = element_rect(fill = 'transparent', color = "transparent"),
      text = element_text(family = "Helvetica"),
      axis.text.y = element_text(color = "black", size = 6),
      axis.title.y = element_text(size = 6),
      axis.ticks.length = unit(-.1, "cm"),
      axis.title.x = element_text(size = 6),
      axis.text.x = element_text(color = "black", size = 6),
      strip.text = element_text(face = "bold", hjust = 0),
      axis.line.x = element_line(linewidth = 0.25),
      axis.line.y = element_line(linewidth = 0.25),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
      strip.background = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    )
  
  # Get significant points and sort by p-value
  sig_points <- subset(data, odd_med > 1 & epval < 0.05)
  sig_points <- sig_points[order(sig_points$epval), ]
  
  # Create a grid for labels with adjustable spread
  label_count <- nrow(sig_points)
  label_mid <- 0.5  # middle of the plot
  label_half_range <- (label_count - 1) / (2 * label_count) * label_spread_factor
  label_positions <- seq(label_mid + label_half_range, 
                         label_mid - label_half_range, 
                         length.out = label_count)
  
  # Calculate the unit length in the x and y directions
  x_range <- diff(layer_scales(p)$x$range$range)
  y_range <- diff(layer_scales(p)$y$range$range)
  unit_length <- sqrt((x_range/2.5)^2 + (y_range/4)^2)
 
  # Add segments and labels
  for (i in 1:nrow(sig_points)) {
    start_x <- sig_points$odd_med[i]
    start_y <- -log10(sig_points$epval[i])
    end_y <- label_positions[i] * 3.2
    segment_length <- segment_length_factor * unit_length
    # Calculate end_x to make segment length 1
    dx <- sqrt(segment_length^2 / (1 + ((end_y - start_y) / (2.5 - start_x))^2))
    end_x <- start_x + dx
    
    p <- p + 
      geom_segment(
        x = start_x,
        y = start_y,
        xend = end_x,
        yend = end_y,
        color = significant_color,
        linewidth = 0.05,
        linetype = "dashed"
      ) +
      annotate(
        "text",
        x = end_x,
        y = end_y,
        label = sig_points$rbp[i],
        hjust = 0,
        vjust = 0.1,
        size = 2.5,
        color = significant_color
      )
  }
  
  return(p)
}

pbs_rbp_eclip_subset <- pbs_rbp_eclip |> group_by(rbp) %>%
  dplyr::filter(observed > 5) |>
  ungroup()  |>
  arrange(desc(odd_med))

fnf_rbp_eclip_subset <- fnf_rbp_eclip |> group_by(rbp) %>%
  dplyr::filter(observed > 5) |>
  ungroup()  |>
  arrange(desc(odd_med))

rbp_pbs <- edited_boxplot_plot(pbs_rbp_eclip_subset, x = 0, y = 0, width = 5, height = 5, 
                    significant_color = "#005587", non_significant_color = "grey70",segment_length_factor = 0.28,
                    label_spread_factor = 0.5)

rbp_fnf <- edited_boxplot_plot(fnf_rbp_eclip_subset, x = 0, y = 0, width = 5, height = 5, 
                               significant_color = "darkorange", non_significant_color = "grey70",segment_length_factor = 0.25,
                               label_spread_factor = 0.4)

save(rbp_pbs, file = "output/results_plots/rbp/rbp_pbs_oddratio.rda")
save(rbp_fnf, file = "output/results_plots/rbp/rbp_fnf_oddratio.rda")


# Finding the MGRN1 for the locus zoom and the gene  
#mgrn1 chr16: 4624826 4690972
#qtl
response_pbs_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_pbs_resInv.rds")
response_fnf_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_fnf_resInv.rds")

pbs_sQtl_bed <- fread("output/Enrichment/rbp/data_prep/significant_pbs_rank0.bed")
colnames(pbs_sQtl_bed) <- c("var_chr", "var_start", "var_end", "var_id", "phe_id", "phe_strd")
fnf_sQtl_bed <- fread("output/Enrichment/rbp/data_prep/significant_fnf_rank0.bed")
colnames(fnf_sQtl_bed) <-  c("var_chr", "var_start", "var_end", "var_id", "phe_id", "phe_strd")

# Convert sQTL data to GRanges objects
pbs_gr <- with(pbs_sQtl_bed, GRanges(var_chr, IRanges(var_start, var_end)))
fnf_gr <- with(fnf_sQtl_bed, GRanges(var_chr, IRanges(var_start, var_end)))


load("external_data/encode_rbp/rbp_prep_rbpOnly/pbs_rbp_elip.Rdata")
load("external_data/encode_rbp/rbp_prep_rbpOnly/fnf_rbp_elip.Rdata")
combined_data_rbp_eclip <- bind_rows(pbs_rbp_eclip, fnf_rbp_eclip)

combined_data_rbp_eclip_subset <- combined_data_rbp_eclip %>%
  group_by(rbp) %>%
  dplyr::filter(observed >5) |>
  #dplyr::filter(!all(epval > 0.05)) %>%
  ungroup()  |>
  arrange(desc(odd_med))


all_subset_enriched_rbps  <- unique(combined_data_rbp_eclip_subset$rbp)
#rbp_bed

rbp_subset.list <- list()
for(i in all_subset_enriched_rbps){
  print(i)
  rbp_bed <- fread(paste0("external_data/encode_rbp/rbp_prep_rbpOnly/",i,".bed"))
  colnames(rbp_bed) <- c("chr","start","end","strand","rbp_nm","cell_line","rep_n","-log2FC","-log10pval")
  rbp_subset.list[[i]] <- rbp_bed
}


# Function to find overlaps and return a data frame
find_overlaps <- function(sqtl_gr, rbp_bed) {
  rbp_gr <- with(rbp_bed, GRanges(chr, IRanges(start, end)))
  overlaps <- findOverlaps(sqtl_gr, rbp_gr)
  data.frame(
    sqtl_index = queryHits(overlaps),
    rbp_index = subjectHits(overlaps)
  )
}


rbp_all.df <- fread("external_data/encode_rbp/rbp_all_dataframe.txt")

#pbs_overlaps <- lapply(rbp_sig.list, function(rbp) find_overlaps(pbs_gr, rbp))
#fnf_overlaps <- lapply(rbp_sig.list, function(rbp) find_overlaps(fnf_gr, rbp))



process_gene_region <- function(gene_gr, gene_name, pbs_sQtl_bed, fnf_sQtl_bed, rbp_all.df, response_pbs_results, response_fnf_results, vsd_geneExp) {
  
  # Find overlaps with RBP data
  find_overlaps_gene <- find_overlaps(gene_gr, rbp_all.df)
  rbp_unique_gene <- rbp_all.df[find_overlaps_gene$rbp_index,] %>% distinct(rbp_nm)
  
  # Find overlaps with sQTL data
  find_overlaps_qtl <- function(snp_gr, qtl_bed) {
    sqtl_gr <- with(qtl_bed, GRanges(var_chr, IRanges(var_start, var_end)))
    overlaps <- findOverlaps(snp_gr, sqtl_gr)
    data.frame(
      snp_index = queryHits(overlaps),
      sqtl_index = subjectHits(overlaps)
    )
  }
  
  pbs_sqtl_gene <- find_overlaps_qtl(gene_gr, pbs_sQtl_bed)
  fnf_sqtl_gene <- find_overlaps_qtl(gene_gr, fnf_sQtl_bed)
  
  fnf_sqtl_overlaps_gene <- fnf_sQtl_bed[fnf_sqtl_gene$sqtl_index,] 
  pbs_sqtl_overlaps_gene <- pbs_sQtl_bed[pbs_sqtl_gene$sqtl_index,] 
  
  # Get unique phe_ids
  gene_phe_id <- c(pbs_sqtl_overlaps_gene$phe_id, fnf_sqtl_overlaps_gene$phe_id) %>% unique()
  
  # Find most significant intron junction
  fnf_mostSig_gene <- response_fnf_results %>% 
    dplyr::filter(phe_id %in% gene_phe_id) %>% 
    slice_min(FNF_p) %>% 
    pull(phe_id)
  
  # Filter RBP expression data
  gene_rbp_exp <- vsd_geneExp %>% dplyr::filter(SYMBOL %in% rbp_unique_gene$rbp_nm)
  
  # Return results as a list
  return(list(
    rbp_exp = gene_rbp_exp,
    phe_id = gene_phe_id,
    mostSig_phe_id = fnf_mostSig_gene
  ))
}

#response_pbs_results |> dplyr::filter(SYMBOL %in% c('CLCN6','PPP1CB','RAF1','LMNA','YWHAB','RAP1A','FAM114A2','BRAF',
#'PFKL','H6PD','PFKP','RPE','PGM2')) |> arrange(PBS_p)

#pbs_results |> dplyr::filter(SYMBOL.y %in% c('CLCN6','PPP1CB','RAF1','LMNA','YWHAB','RAP1A','FAM114A2','BRAF',
#                                           'PFKL','H6PD','PFKP','RPE','PGM2' )) |> dplyr::filter( minor_alle_count >= 5, dist_phe_var < 1000) |> data.frame()


#make gr for the test 
snhg29_gr <- GRanges(seqnames = "chr17", IRanges(16439518,16439519))
pcyt1a_gr <- GRanges(seqnames = "chr3", IRanges(196287494,196287495))
mica_gr <- GRanges(seqnames = "chr4", IRanges(37855720,37855721))
# calu_gr <- GRanges(seqnames = "chr7", IRanges(128756527,128756528))


# For SNHG29
snhg29_results <- process_gene_region(snhg29_gr, "SNHG29", pbs_sQtl_bed, fnf_sQtl_bed, rbp_all.df, response_pbs_results, response_fnf_results, vsd_geneExp)
pcyt1a_results <- process_gene_region(pcyt1a_gr, "PCYT1A", pbs_sQtl_bed, fnf_sQtl_bed, rbp_all.df, response_pbs_results, response_fnf_results, vsd_geneExp)
MICA_results <- process_gene_region(mica_gr, "MICA", pbs_sQtl_bed, fnf_sQtl_bed, rbp_all.df, response_pbs_results, response_fnf_results, vsd_geneExp)
# calu_results <- process_gene_region(calu_gr, "CALU", pbs_sQtl_bed, fnf_sQtl_bed, rbp_all.df, response_pbs_results, response_fnf_results, vsd_geneExp)

# For SNHG29
snhg29_rbp_exp <- snhg29_results$rbp_exp
snhg29_phe_id <- snhg29_results$phe_id
fnf_mostSig_snhg29 <- snhg29_results$mostSig_phe_id



# Signal data 
dir_merged <- "/work/users/s/e/seyoun/CQTL_sQTL/output/signals/merged_norm/"
ctl_signal <- paste0(dir_merged,"CTL_norm.bw")
fnf_signal <- paste0(dir_merged,"FNF_norm.bw")


# Find the right RBP signals ID
encode_meta_signals <- fread("external_data/encode_rbp/signal/metadata.tsv")

RBP_name <- "AATF"
encode_rbp_AATF <- encode_meta_signals |> dplyr::filter(str_detect(`Experiment target` ,RBP_name), `File assembly` == "GRCh38", 
                                                            `Output type` == "plus strand signal of unique reads", `Biological replicate(s)` == 1 ) |> 
  mutate(filename  = paste0("/work/users/s/e/seyoun/CQTL_sQTL/external_data/encode_rbp/signal/",`File accession`, ".bigWig")) |> pull(filename) 


RBP_name <- "EFTUD2"
encode_rbp_eftud2 <- encode_meta_signals |> dplyr::filter(str_detect(`Experiment target` ,RBP_name), `File assembly` == "GRCh38", 
                                                        `Output type` == "minus strand signal of unique reads", `Biological replicate(s)` == 1 ) |> 
  mutate(filename  = paste0("/work/users/s/e/seyoun/CQTL_sQTL/external_data/encode_rbp/signal/",`File accession`, ".bigWig")) |> pull(filename) 

# RBP_name <- c("AATF|CSTF2|TRA2A")
# encode_rbp_MICA <- encode_meta_signals |> dplyr::filter(str_detect(`Experiment target` ,RBP_name), `File assembly` == "GRCh38", 
#                                                         `Output type` == "plus strand signal of unique reads", `Biological replicate(s)` == 1, !str_detect(`Experiment target`, "CSTF2T") ) |> 
#   mutate(filename  = paste0("/work/users/s/e/seyoun/CQTL_sQTL/external_data/encode_rbp/signal/",`File accession`, ".bigWig")) |> pull(filename) 
# 


# plot gardener to keep both plot at the same time------------------------------
pdf(file = "output/results_plots/rbp/main_figure4.pdf",   # The directory you want to save the file in
    width = 9.5, # The width of the plot in inches
    height = 6)

pageCreate(width = 9.5, height =6 , default.units = "inches", showGuides = TRUE)


# scatter plots of PBS and FNF--------------------------------------------------
plotText("A", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

#save(rbp_pbs, file = "output/results_plots/rbp/rbp_pbs_oddratio.rda")
#save(rbp_fnf, file = "output/results_plots/rbp/rbp_fnf_oddratio.rda")
load("output/results_plots/rbp/rbp_pbs_oddratio.rda")
load("output/results_plots/rbp/rbp_fnf_oddratio.rda")
plotGG(rbp_pbs, x = 0.5, y= 0.65, width = 2.6, height = 2.5, just = c("left","top"))
plotGG(rbp_fnf, x = 0.5, y= 3.3, width = 2.6, height = 2.5, just = c("left","top"))


plotText("Enrichment of sSNPs in RNA binding protein sites", x = 0.5, y = 0.4, just = "left", fontfamily = "Helvetica",
         fontsize = 8, fontcolor ="black", fontface="bold" )

plotText("PBS", x = 1.8, y = 0.65, just = "center", fontfamily = "Helvetica",
         fontsize = 8, fontcolor ="#005587" )

plotText("FN-f", x = 1.8, y = 3.3, just = "center", fontfamily = "Helvetica",
         fontsize = 8, fontcolor ="darkorange" )


# visualize locus zoom plot

plotText("B", x = 3.4, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

plot_sQTL_manhattan("SNHG29", "FNF", x_start = 4.0, y_start = 0.4, width = 2.2, height = 2,zoom_range = 100000)

#SNHG29
region_pg <- pgParams(assembly = "hg38",chrom = "chr17",
                      chromstart = 16439519-5000,
                      chromend = 16439519+5000,
                      x = 4.0, width = 2.2)

plotgenes <- plotGenes(params = region_pg,
                       y = 3.4,x=4.0,
                       height = 0.25,
                       geneHighlights = data.frame("gene" = 'SNHG29',
                                                   "color" = "#37a7db"),fontsize = 6,geneOrder="SNHG29",strandLabels = FALSE,
                       geneBackground = "transparent")

RNA_signals <- plotSignal(encode_rbp_AATF,
                               params = region_pg, x=4.0,
                               y = 3.7, height = 0.25 ,linecolor = yl_gn_bu[9], 
                               fill= yl_gn_bu[9],
                               default.units = "inches")

# Extract the y-range
y_range <- RNA_signals$range
# Format as a label with 2-digit precision
rna_scale_label <- sprintf("[%.1f - %.1f]", y_range[1], y_range[2])

plotText(label = rna_scale_label,
         x = 4.0, 
         y = 3.7 + 0.25 + 0.05,  # place it slightly above the signal track
         just = c("left", "bottom"),
         fontsize = 6,
         rot = 90,            # rotate the text 90 degrees
         default.units = "inches")


annoGenomeLabel(plot = RNA_signals, params = region_pg, fontsize = 6, 
                y = 4.0)

annoHighlight(
  plot = RNA_signals,
  chrom = "chr17",
  chromstart = 16439518-25, chromend = 16439519+25,
  y =3.45, height = 0.5 ,just = c("left", "top"),
    default.units = "inches", fill = "red",
)

plotText("AATF", x=4.0, y= 3.8, fontfamily = "Helvetica", just = "left", fontsize = 6, fontcolor = yl_gn_bu[9])
plotText("rs11871968", x= 5.1,y=3.4, fontfamily = "Helvetica", just = "center", fontsize = 6,fontcolor = "red")

plotRect(
  x = 3.9, y = 3.32, width = 2.4, height = 0.84,
  just = c("left", "top"), default.units = "inches",
  lwd = 1, fill = "transparent",lty = 2 ,linecolor = "grey40"
)

load("output/results_plots/rbp/rbp_snhg29_Genoboxplot.rda")
load("output/results_plots/rbp/snhg29_rbpCorrPlot.rda")

plotGG(snhg29_Genoboxplot,x = 3.65, y = 4.4, width = 1.4, height = 1.4)
plotGG(snhg29_rbpCorrPlot,x = 5, y = 4.4, width = 1.4, height = 1.4)
plotText("chr17:16439414-16439528", x= 5.1,y=4.3, fontfamily = "Helvetica", just = "center", fontsize = 7)

plotText("C", x = 6.4, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

# PCYT1A
plot_sQTL_manhattan("PCYT1A", "FNF", x_start =7.0, y_start = 0.4, width = 2.2, height = 2,zoom_range = 100000)

region_pg <- pgParams(assembly = "hg38",chrom = "chr3",
                      chromstart =196287495-5000,
                      chromend = 196287495+5000 ,
                      x = 7.0, width = 2.2)


plotgenes <- plotGenes(params = region_pg,
                       y = 3.4,x=7.0,
                       height = 0.25,
                       geneHighlights = data.frame("gene" = 'PCYT1A',
                                                   "color" = "#37a7db"),fontsize = 6,geneOrder="PCYT1A",
                       strandLabels = FALSE,geneBackground = "transparent")

eftud2_pcy1a <- readBigwig(encode_rbp_eftud2[2],
  chrom = "chr3",
  chromstart =196287495-5000,
  chromend = 196287495+5000 
) |> mutate(score =abs(score))


RNA_signals <- plotSignal(eftud2_pcy1a,
                          params = region_pg, x=7.0,
                          y = 3.7, height = 0.25 ,linecolor = yl_gn_bu[9], 
                          fill= yl_gn_bu[9],
                          default.units = "inches")


# Extract the y-range
y_range <- RNA_signals$range
# Format as a label with 2-digit precision
rna_scale_label <- sprintf("[%.1f - %.1f]", y_range[1], y_range[2])

plotText(label = rna_scale_label,
         x = 7.0, 
         y = 3.7 + 0.25 + 0.05,  # place it slightly above the signal track
         just = c("left", "bottom"),
         fontsize = 6,
         rot = 90,            # rotate the text 90 degrees
         default.units = "inches")



annoGenomeLabel(plot = RNA_signals, params = region_pg, fontsize = 6, 
                y = 4.0)

annoHighlight(
  plot = plotgenes,
  chrom = "chr3",
  chromstart = 196287495-25, chromend = 196287495+25,
  y =3.45, height = 0.5 ,just = c("left", "top"),
  default.units = "inches", fill = "red",
)


plotText("EFTUD2", x=7.0, y= 3.8, fontfamily = "Helvetica", just = "left", fontsize = 6, fontcolor = yl_gn_bu[9] )

plotText("rs6809764", x= 8.1,y=3.4, fontfamily = "Helvetica", just = "center", fontsize = 6,fontcolor = "red")

plotRect(
  x = 6.9, y = 3.32, width = 2.4, height = 0.84,
  just = c("left", "top"), default.units = "inches",
  lwd = 1, fill = "transparent",lty = 2 ,linecolor = "grey40"
)


load("output/results_plots/rbp/PCYT1A_Genoboxplot.rda")
load("output/results_plots/rbp/PCYT1A_rbpCorrPlot.rda")
plotGG(PCYT1A_Genoboxplot,x = 6.65, y = 4.4, width = 1.4, height = 1.4)
plotGG(PCYT1A_rbpCorrPlot,x = 8.0, y = 4.4, width = 1.4, height = 1.4)
plotText(pcyt1a_results$phe_id[2], x= 8.1,y=4.3, fontfamily = "Helvetica", just = "center", fontsize = 7)


dev.off()


#Find the y-axis AATF AND eftud2

region_pg <- pgParams(assembly = "hg38",chrom = "chr17",
                      chromstart = 16439519-5000,
                      chromend = 16439519+5000,
                      x = 4.0, width = 2.2)

encode_aatf <- readBigwig(encode_rbp_AATF, params = region_pg)


region_pg <- pgParams(assembly = "hg38",chrom = "chr3",
                      chromstart =196287495-5000,
                      chromend = 196287495+5000 ,
                      x = 7.0, width = 2.2)

encode_eftud2 <- readBigwig(encode_rbp_eftud2[2])
