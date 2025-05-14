# rMATs data
setwd("/work/users/s/e/seyoun/CQTL_sQTL")  

library(maser)
library(rtracklayer)
library(dplyr)
library(ggrepel)
library(stringr)
library(ggplot2)
#-------------------------------------------------------------------------------
# Comparison wd vs hetKD
edited_path <- "output/rmats_edited"
kd_wd <- maser(edited_path, c("kd", "wd"), ftype = "JCEC")
geneEvents(kd_wd, "SNRNP70")
#head(summary(KD_WD, type = "SE")[, 1:8])
kd_filt <- filterByCoverage(kd_wd, avg_reads = 5)
sig_kd <- topEvents(kd_filt, fdr = 0.05, deltaPSI = 0.1)
summary(sig_kd,'SE') |> dim()
geneEvents(sig_kd, "SNRNP70")
rmats_kd_genes <- c(summary(sig_kd,'SE')[,"geneSymbol"],
                    summary(sig_kd,'RI')[,"geneSymbol"],
                    summary(sig_kd,'A3SS')[,"geneSymbol"],
                    summary(sig_kd,'A5SS')[,"geneSymbol"],
                    summary(sig_kd,'MXE')[,"geneSymbol"]) |> unique()

save(kd_wd, sig_kd, kd_filt, rmats_kd_genes,
 file = "output/rmats_edited/sig_kd_wd.RData")
#------------------------------------------------------------------------------
# Visualize the MA? plot looking? 

# Function to process PSI values and calculate mean PSI
process_psi <- function(psi_string) {
  psi_vals <- as.numeric(str_split(psi_string, ",")[[1]])
  mean(psi_vals, na.rm = TRUE)
}

# Function to combine all event types with proper junction formatting
combine_events <- function(sig_kd) {
  # Skipped Exon (SE)
  se <- as.data.frame(summary(sig_kd, "SE")) %>%
    mutate(
      Junction = paste0(Chr, ":", exon_target, ":SE"),
      EventType = "SE"
    )
  
  # Retained Intron (RI)
  ri <- as.data.frame(summary(sig_kd, "RI")) %>%
    mutate(
      Junction = paste0(Chr, ":", exon_ir, ":RI"),
      EventType = "RI"
    )
  
  # Alternative 3' Splice Site (A3SS)
  a3ss <- as.data.frame(summary(sig_kd, "A3SS")) %>%
    mutate(
      Junction = paste0(Chr, ":", exon_long, ":A3SS"),
      EventType = "A3SS"
    )
  
  # Alternative 5' Splice Site (A5SS)
  a5ss <- as.data.frame(summary(sig_kd, "A5SS")) %>%
    mutate(
      Junction = paste0(Chr, ":", exon_long, ":A5SS"),
      EventType = "A5SS"
    )
  
  # Mutually Exclusive Exons (MXE)
  mxe <- as.data.frame(summary(sig_kd, "MXE")) %>%
    mutate(
      Junction = paste0(Chr, ":", exon_1, "-", exon_2, ":MXE"),
      EventType = "MXE"
    )
  
  bind_rows(se, ri, a3ss, a5ss, mxe)
}

# Get all significant events
all_events <- combine_events(wd_kd)

# Process the data
plot_data <- all_events %>%
  mutate(
    # Calculate mean PSI for each condition
    PSI_1_mean = sapply(PSI_1, process_psi),
    PSI_2_mean = sapply(PSI_2, process_psi),
    # Correct IncLevelDifference to PSI_2 - PSI_1 (KD - WD)
    IncLevelDifference_corrected = PSI_2_mean - PSI_1_mean,
    # Calculate overall mean PSI across all samples
    Mean_PSI = (PSI_1_mean + PSI_2_mean) / 2
  ) %>%
  dplyr::select(Junction, Mean_PSI, IncLevelDifference_corrected, FDR, geneSymbol, EventType)

# Define significance thresholds
up_significant <- plot_data %>%
  filter(IncLevelDifference_corrected > 0.1, FDR < 0.05)

down_significant <- plot_data %>%
  filter(IncLevelDifference_corrected < -0.1, FDR < 0.05)

# SNRNP70 specific data
snrnp70_data <- plot_data %>% filter(geneSymbol == "SNRNP70") |> filter(abs(IncLevelDifference_corrected) > 0.1, FDR < 0.05)

# Create the MA plot
ggplot(plot_data, aes(x = Mean_PSI, y = IncLevelDifference_corrected)) +
  # Non-significant points
  geom_point(data = plot_data %>% filter(!Junction %in% c(up_significant$Junction, down_significant$Junction)),
             color = "gray", alpha = 0.5, size = 1) +
  # Up-regulated points
  geom_point(data = up_significant,
             color = "#e07653", alpha = 0.5, size = 1) +
  # Down-regulated points
  geom_point(data = down_significant,
             color = "#1e87a5", alpha = 0.5, size = 1) +
  # Significant points with outline
  geom_point(data = up_significant,
             color = "black", fill = "#e07653", size = 2, shape = 21, stroke = 0.5) +
  geom_point(data = down_significant,
             color = "black", fill = "#1e87a5", size = 2, shape = 21, stroke = 0.5) +
  # Labels for significant up points
  # ggrepel::geom_text_repel(data = up_significant,
  #                          aes(label = geneSymbol),
  #                          color = "black",
  #                          size = 2.5,
  #                          segment.color = "black",
  #                          segment.size = 0.2,
  #                          box.padding = 0.5,
  #                          point.padding = 0.5,
  #                          max.overlaps = Inf) +
  # # Labels for significant down points
  # ggrepel::geom_text_repel(data = down_significant,
  #                          aes(label = geneSymbol),
  #                          color = "black",
  #                          size = 2.5,
  #                          segment.color = "black",
  #                          segment.size = 0.2,
  #                          box.padding = 0.5,
  #                          point.padding = 0.5,
  #                          max.overlaps = Inf) +
  # SNRNP70 points and labels
  geom_point(data = snrnp70_data,
             color = "black", fill = "purple", size = 2.5, shape = 21, stroke = 0.5) +
  ggrepel::geom_text_repel(data = snrnp70_data,
                           aes(label = geneSymbol),
                           color = "black",
                           size = 2.5,
                           segment.color = "black",
                           segment.size = 0.2,
                           box.padding = 0.5,
                           point.padding = 0.5,
                           max.overlaps = Inf) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = 'transparent', color = "transparent"),
    plot.background = element_rect(fill = 'transparent', color = "transparent"),
    text = element_text(family = "Helvetica"),
    legend.position = "none",
    axis.text = element_text(color = "black", size = 6),
    axis.title = element_text(size = 6),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
    axis.line = element_line(linewidth = 0.25)
  ) +
  labs(
    x = "Mean PSI",
    y = "Delta PSI (KD/WD)",
    title = "MA Plot: All Splicing Events"
  )

#-------------------------------------------------------------------------------
# Compare with the sGene only?
load("/work/users/s/e/seyoun/CQTL_sQTL/output/clu_fnf/introns_fnf_joinAll")
load("/work/users/s/e/seyoun/CQTL_sQTL/output/clu_oa/introns_oa_joinAll")
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

  
  return(introns_processed)
}

# Process the data
fnf_data <- process_introns_data(introns_fnf_pval_include)
oa_data <- process_introns_data(introns_oa_pval_include)


#------------------------------------------------------------------------------
# Try another way? SNRNP70 gene expression to see the correlation with the PSI. 







