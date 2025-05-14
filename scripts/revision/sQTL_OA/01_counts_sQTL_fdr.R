# Calculate the counts
#Libraries
library(ggplot2)
library(dplyr)
library(tidyr)

#-------------------------------------------------------------------------------
# Define the possible arguments files, significant p-values etc
adjusted_beta_p <- 0.1
pbs_perm_dir <- "output/03.qtltools_re_oa/perm_pbs"
oa_perm_dir <- "output/03.qtltools_re_oa/perm_oa"


# Function for reading and processing files --------------------------------
process_files <- function(filepath_pattern, pattern_suffix) {
  filepaths <- list.files(filepath_pattern, pattern = pattern_suffix, full.names = TRUE)
  
  # Sort filepaths numerically based on the PC number
  sorted_filepaths <- filepaths[order(as.numeric(gsub(".*pc(\\d+)_.*", "\\1", filepaths)))]
  
  # Read files and process them
  list_processed <- lapply(sorted_filepaths, function(filepath) {
    df <- fread(filepath, header = FALSE, stringsAsFactors = FALSE)
    
    # Filter out rows with "chrX" or "chrY" from the specified column
    df_filtered <- df[!(df[, 2] %in% c("chrX", "chrY")), ] |> filter(!is.na(V20))
    
    # Add adjusted p-values using FDR
    df_filtered$p_adjusted <- p.adjust(df_filtered$V20, method = "fdr")
    
    # Filter based on significance threshold
    df_sig <- df_filtered[df_filtered$p_adjusted <= adjusted_beta_p, ]
    
    return(df_sig)
  })
  
  return(list_processed)
}

# Process PBS and OA files
pbs_list_processed <- process_files(pbs_perm_dir, "_allchr.pbs.perm")
oa_list_processed <- process_files(oa_perm_dir, "_allchr.oa.perm")

# Count unique values in column 1 of each list element (intron clusters)
pbs_unique_counts <- lapply(pbs_list_processed, function(x) {
  length(unique(x$V1))   
})

oa_unique_counts <- lapply(oa_list_processed, function(x) {
  length(unique(x$V1)) 
})

# Create a data frame with the list elements and names
df <- data.frame(
  pc = paste0("pc", 1:20), 
  pbs_counts = unlist(pbs_unique_counts), 
  oa_counts = unlist(oa_unique_counts)
)

# Format for plotting
df$pc <- factor(df$pc, levels = paste0("pc", 1:20))

df_long <- df %>%
  pivot_longer(
    cols = c(pbs_counts, oa_counts),
    names_to = "type",
    values_to = "counts"
  ) %>%
  mutate(type = dplyr::recode(type, pbs_counts = "PBS", oa_counts = "OA"))

df_long$type <- factor(df_long$type, levels = c("PBS", "OA"))

# Calculate the max for each group
max_counts <- df_long %>%
  dplyr::group_by(type) %>%
  dplyr::summarize(
    max_count = max(counts), 
    pc_max = pc[which.max(counts)]
  ) %>%
  ungroup()

# Set colors for PBS and OA
max_counts$color <- ifelse(max_counts$type == "PBS", "#18417c", "#aa5e39")

# Create the dot plot
counts_dotPlot <- ggplot(df_long, aes(x = pc, y = counts, color = type)) +
  geom_point(size = 2) +
  labs(x = "PC", y = "Number of significant introns", color = "Type") +
  scale_color_manual(values = c("PBS" = "#2057A7", "OA" = "#d6663c")) +
  # Adjust y-axis limits as needed based on your data
  # scale_y_continuous(limits = c(2500, 6000)) +
  theme_classic() +
  geom_text(
    data = max_counts, 
    aes(x = pc_max, y = max_count + 100, label = paste0(pc_max, ": ", max_count)),
    color = max_counts$color, 
    vjust = 0, 
    size = 3.5
  ) +
  theme(
    text = element_text(family = "Helvetica"),
    panel.background = element_rect(fill = 'transparent', color = "transparent"),
    plot.background = element_rect(fill = 'transparent', color = "transparent"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    legend.background = element_blank(),
    legend.position = c(0.9, 0.9),
    axis.title = element_text(size = 9),
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.ticks.length.y = unit(-0.1, "cm"),
    axis.ticks.y = element_line(color = "black", linewidth = 0.25),
    axis.text = element_text(color = "black", size = 8),
    axis.ticks.x = element_line(color = "black", linewidth = 0.25),
    axis.text.x = element_text(color = "black", size = 8, angle = 45),
    plot.margin = margin(t = 5, r = 10, b = 5, l = 5)
  )

# Save the plot
ggsave(
  "output/results_plots/PBS_OA_Significant_introns.pdf",
  counts_dotPlot,
  width = 11,
  height = 6.5,
  dpi = 300,
  bg = "transparent"
)

# Print the plot
print(counts_dotPlot)

# Output the PC with maximum counts for each group
print(max_counts)