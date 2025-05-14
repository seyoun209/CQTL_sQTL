# SNRNP70 edited splicing change (Test)
setwd("/work/users/s/e/seyoun/crispr/02.test_seq")
library(leafcutter)
library(stringr)
library(magrittr)
library(ggplot2)
library(dplyr)
library(ggrepel)

#Leafcutter

#Functions
normalize_column <- function(x) {
  x / sum(x, na.rm = TRUE)
}
#-------------------------------------------------------------------------------

wd_kd <- load("clu_wd_kd/WD_KD.Rdata")

#-------------------------------------------------------------------------------
#Step1 Leafcutter Calculate the PSI ratio

process_count_data <- function(count_data, output_suffix = "_processed", remove_original = FALSE) {
  # Store original name
  original_name <- deparse(substitute(count_data))
  
  # Construct new name with suffix
  new_name <- paste0(original_name, output_suffix)
  
  # Process the data
  result <- count_data %>%
    mutate(clu = str_split_fixed(rownames(.), ":", 4)[,4]) %>%
    group_by(clu) %>%
    mutate_all(normalize_column) %>%
    ungroup() %>%
    as.data.frame() %>%
    set_rownames(rownames(count_data)) %>%
    dplyr::select(-clu)
  
  # QC filtering (optional - only if you want to include the NA handling)
  result_qc <- result[rowMeans(is.na(result)) <= 0.4,,drop=F]
  
  # Handle NAs with row means
  row_means <- rowMeans(result_qc, na.rm = TRUE)
  row_means_outer <- outer(row_means, rep(1, ncol(result_qc)))
  result_qc[is.na(result_qc)] <- row_means_outer[is.na(result_qc)]
  
  # Add Junction column
  result_qc <- cbind(rownames(result_qc), result_qc)
  colnames(result_qc)[1] <- "Junction"
  
  # Assign to global environment
  assign(new_name, result_qc, envir = .GlobalEnv)
  
  # Remove original if specified
  if(remove_original && exists(original_name, envir = .GlobalEnv)) {
    remove(list = original_name, envir = .GlobalEnv)
  }
  
  # Return the processed data
  return(result_qc)
}

ratios_edited <- process_count_data(counts, "_edited_v1", FALSE)

#MA plot -----------------------------------------------------------------------
# Calculate mean normalized counts across all samples for each junction
counts_mean <- rowMeans(counts_edited_v1[,2:5], na.rm = TRUE)  # Assuming columns 2-5 are your sample columns

# Create a data frame with junction IDs and mean counts
plot_data <- data.frame(
  Junction = rownames(counts_edited_v1),
  Mean_Counts = counts_mean
)

# Merge with introns_edited_V1 data
# First, create a Junction column in introns_edited_V1 that matches counts_edited_v1 format
introns <- introns %>%
  mutate(Junction = paste0(chr, ":", start, ":", end, ":", clusterID))

# Prepare the data
# Prepare the initial plot data
plot_data <- data.frame(
  Junction = rownames(counts_edited_v1),
  Mean_Counts = counts_mean
) %>%
  inner_join(introns %>% 
               mutate(Junction = paste0(chr, ":", start, ":", end, ":", clusterID)) %>%
               dplyr::select(Junction, deltapsi, gene, clusterID), 
             by = "Junction")

# Get max absolute deltaPSI per cluster
plot_data_max <- plot_data %>%
  group_by(clusterID) %>%
  mutate(abs_deltapsi = abs(deltapsi)) %>%
  slice_max(abs_deltapsi, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(-abs_deltapsi)


plot_data_max <- plot_data_max %>%
  left_join(dplyr::select(clusters, clusterID, FDR), by = "clusterID")

# Define significance thresholds based on FDR and deltaPSI
up_significant <- plot_data_max %>%
  filter(deltapsi > 0.1, FDR < 0.05)

down_significant <- plot_data_max %>%
  filter(deltapsi < -0.1, FDR < 0.05)

# SNRNP70 specific data
snrnp70_data <- plot_data_max %>% filter(gene == "SNRNP70")

ggplot(plot_data_max, aes(x = Mean_Counts, y = deltapsi)) +
  # Non-significant points
  geom_point(data = plot_data_max %>% filter(!Junction %in% c(up_significant$Junction, down_significant$Junction)),
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
  ggrepel::geom_text_repel(data = up_significant,
                           aes(label = gene),
                           color = "black",
                           size = 2.5,
                           segment.color = "black",
                           segment.size = 0.2,
                           box.padding = 0.5,
                           point.padding = 0.5,
                           max.overlaps = Inf) +
  # Labels for significant down points
  ggrepel::geom_text_repel(data = down_significant,
                           aes(label = gene),
                           color = "black",
                           size = 2.5,
                           segment.color = "black",
                           segment.size = 0.2,
                           box.padding = 0.5,
                           point.padding = 0.5,
                           max.overlaps = Inf) +
  # SNRNP70 points and labels
  geom_point(data = snrnp70_data,
             color = "black", fill = "purple", size = 2.5, shape = 21, stroke = 0.5) +
  ggrepel::geom_text_repel(data = snrnp70_data,
                           aes(label = gene),
                           color = "black",
                           size = 3,
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
    y = "Delta PSI (Max per Cluster)",
    title = ""
  )

#-------------------------------------------------------------------------------
