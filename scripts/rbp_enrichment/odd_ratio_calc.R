library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(ggplot2)
library(ggtext)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(grid)
library(colorspace)

setwd("/work/users/s/e/seyoun/CQTL_sQTL/")



#Fenrich analysis odd ratio

# Function to process a single file
process_file <- function(file_path) {
  # Extract tissue, RBP, and condition from filename
  #file_info <- str_match(basename(file_path), "(.+)_(.+)_(.+)_sQTL_enrichment.txt")
  file_info <- str_match(basename(file_path), "(.+)_(.+)_sQTL_enrichment.txt")
  #tissue <- file_info[2]
  rbp <- file_info[2]
  condition <- file_info[3]
  
  # Read the file
  data <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)
  
  # Check if the file is empty or has unexpected content
  if (nrow(data) == 0 || ncol(data) != 8) {
    warning(paste("Skipping file due to unexpected format:", file_path))
    return(NULL)
  }
  
  # Rename columns
  colnames(data) <- c("observed", "total_qtls", "expected_mean", "expected_sd", "epval", "odd_ratio_dnv", "odd_ratio_med", "odd_ratio_upv")
  
  # Calculate fold change
  fold_change <- data$observed / data$expected_mean
  
  # Calculate odds ratio, confidence interval, and p-value
  fisher_result <- tryCatch({
    fisher.test(matrix(c(data$observed, data$total_qtls, round(data$expected_mean), data$total_qtls), ncol=2))
  }, error = function(e) {
    warning(paste("Error in Fisher's test for file:", file_path, "\nError:", e$message))
    return(list(conf.int = c(NA, NA), estimate = NA, p.value = NA))
  })

  
  # Create a data frame with all required information
  result <- data.frame(
    observed = data$observed,
    total_qtls = data$total_qtls,
    expected_mean = data$expected_mean,
    expected_sd = data$expected_sd,
    #tissue = tissue,
    rbp = rbp,
    condition = condition,
    fold_change = fold_change,
    epval = data$epval,
    odd_dnv = data$odd_ratio_dnv,
    odd_med = data$odd_ratio_med,
    odd_upv= data$odd_ratio_upv,
    ci_lower = fisher_result$conf.int[1],
    ci_upper = fisher_result$conf.int[2],
    odds_ratio = fisher_result$estimate,
    fpval = fisher_result$p.value
  )
  
  return(result)
}

# Process all files in PBS directory
#pbs_files <- list.files("output/Enrichment/rbp/PBS", full.names = TRUE, pattern = "*.txt")
pbs_files <- list.files("output/Enrichment/rbp_all/PBS", full.names = TRUE, pattern = "*.txt")
pbs_rbp_results <- map_df(pbs_files, process_file)

# Process all files in FNF directory
fnf_files <- list.files("output/Enrichment/rbp_all/FNF", full.names = TRUE, pattern = "*.txt")
fnf_rbp_results <- map_df(fnf_files, process_file)

pbs_rbp_results$ratio_dnv_med <- pbs_rbp_results$odd_dnv / pbs_rbp_results$odd_med
fnf_rbp_results$ratio_dnv_med <- fnf_rbp_results$odd_dnv / fnf_rbp_results$odd_med

pbs_rbp_results$ratio_upv_med <- pbs_rbp_results$odd_upv / pbs_rbp_results$odd_med
fnf_rbp_results$ratio_upv_med <- fnf_rbp_results$odd_upv / fnf_rbp_results$odd_med

# Set thresholds (you can adjust these based on your data and requirements)
threshold_dnv <- 1.5  # For example, odd_dnv shouldn't be more than 2.5 times odd_med
threshold_upv <- 1.5  # For example, odd_upv shouldn't be more than 1.5 times odd_med

pbs_rbp_results_filtered <- pbs_rbp_results[which(pbs_rbp_results$ratio_dnv_med <= threshold_dnv & pbs_rbp_results$ratio_upv_med <= threshold_upv), ]
fnf_rbp_results_filtered <- fnf_rbp_results[which(fnf_rbp_results$ratio_dnv_med <= threshold_dnv & fnf_rbp_results$ratio_upv_med <= threshold_upv), ]


combined_data_rbp_eclip <- bind_rows(pbs_rbp_results_filtered, fnf_rbp_results_filtered)

combined_data_rbp_sig <- combined_data_rbp_eclip %>%
  group_by(rbp) %>%
  dplyr::filter(observed >5) |>
  dplyr::filter(!all(epval > 0.05)) %>%
  ungroup()  |>
  arrange(desc(odd_med))


# Create a position dodge object
dodge <- position_dodge(width = 0.9)

# Create the combined plot
ggplot(combined_data_rbp_sig, aes(x = odd_med, y = rbp, color = condition, group = condition)) +
  geom_vline(aes(xintercept = 1), linewidth = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmin = odd_dnv, xmax = odd_upv), 
                 size = 0.2, height = 0.4, position = position_dodge(width = 0.7), color = "grey50") +
  geom_point(size = 3, aes(alpha = ifelse(epval < 0.05, "< 0.05", "≥ 0.05")), 
             position = position_dodge(width = 0.7)) +
  scale_x_continuous(
    breaks = c(0, 0.5, 1, 1.5, 2),
    labels = c(0, 0.5, 1, 1.5, 2))+
  scale_color_manual(values = c("fnf" = "#FFB81C", "pbs" = "#005587"), 
                     labels = c("pbs" = "PBS", "fnf" = "FN-f"),
                     name = "Group",
                     breaks = c("pbs", "fnf")) +
  scale_alpha_manual(values = c("< 0.05" = 1, "≥ 0.05" = 0.3), name = "p-value") +
  coord_cartesian(xlim = c(0, 2.2)) +
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
    plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
    strip.background = element_blank(),
    axis.ticks = element_blank(),
    legend.direction = "horizontal",
    legend.title = element_text(size=6),
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "bottom",
    legend.box.spacing = unit(0.1,"mm")
  ) +
  labs(y = "", x = "log(odds ratio)", 
       title = "Enrichment of sSNPs in RNA binding protein sites",
       color = "Group",
       alpha = "p-value") +
  guides(color = guide_legend(order = 1),
         alpha = guide_legend(order = 2))





#-------------------------------------------------------------------------------

# Function to process a single file

process_chromhmm_file <- function(file_path) {
  # Extract state and condition from filename
  file_info <- str_match(basename(file_path), "(.+)_(.+)_sQTL_enrichment.txt")
  state <- file_info[2]
  condition <- file_info[3]
  
  # Read the file
  data <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)
  
  # Check if the file is empty or has unexpected content
  if (nrow(data) == 0 || ncol(data) != 8) {
    warning(paste("Skipping file due to unexpected format:", file_path))
    return(NULL)
  }
  
  # Rename columns
  colnames(data) <- c("observed", "total_qtls", "expected_mean", "expected_sd", "epval", "odd_ratio_dnv", "odd_ratio_med", "odd_ratio_upv")
  
  # Calculate fold change
  fold_change <- data$observed / data$expected_mean
  
  # Calculate odds ratio, confidence interval, and p-value
  fisher_result <- tryCatch({
    fisher.test(matrix(c(data$observed, data$total_qtls, round(data$expected_mean), data$total_qtls), ncol=2))
  }, error = function(e) {
    warning(paste("Error in Fisher's test for file:", file_path, "\nError:", e$message))
    return(list(conf.int = c(NA, NA), estimate = NA, p.value = NA))
  })
  
  # Create a data frame with all required information
  result <- data.frame(
    observed = data$observed,
    total_qtls = data$total_qtls,
    expected_mean = data$expected_mean,
    expected_sd = data$expected_sd,
    state = state,
    condition = condition,
    fold_change = fold_change,
    epval = data$epval,
    odd_dnv = data$odd_ratio_dnv,
    odd_med = data$odd_ratio_med,
    odd_upv= data$odd_ratio_upv,
    ci_lower = fisher_result$conf.int[1],
    ci_upper = fisher_result$conf.int[2],
    odds_ratio = fisher_result$estimate,
    fpval = fisher_result$p.value
  )
  
  return(result)
}
# Set the base directory
base_dir <- "/work/users/s/e/seyoun/CQTL_sQTL/output/Enrichment/chromhmm"

# Process all files in PBS directory
pbs_files_chm <- list.files(file.path(base_dir, "PBS"), full.names = TRUE, pattern = "*.txt")
pbs_chromhmm_results <- map_df(pbs_files_chm, process_chromhmm_file, .id = "file_id")

# Process all files in FNF directory
fnf_files_chm <- list.files(file.path(base_dir, "FNF"), full.names = TRUE, pattern = "*.txt")
fnf_chromhmm_results <- map_df(fnf_files_chm, process_chromhmm_file, .id = "file_id")

pbs_chromhmm_results_addGroup <- pbs_chromhmm_results %>%
  mutate(group = case_when(
    state %in% c("1_TssA", "2_TssAFlnk", "10_TssBiv", "11_BivFlnk") ~ "Promoter",
    state %in% c("3_TxFlnk", "4_Tx", "5_TxWk") ~ "Transcription",
    state %in% c("6_EnhG", "7_Enh", "12_EnhBiv") ~ "Enhancer",
    state %in% c("8_ZNF/Rpts", "9_Het", "13_ReprPC", "14_ReprPCWk", "15_Quies","8_ZNF_Rpts") ~ "Inactive chromatin",
    TRUE ~ "Other"
  ))


fnf_chromhmm_results_addGroup <- fnf_chromhmm_results %>%
  mutate(group = case_when(
    state %in% c("1_TssA", "2_TssAFlnk", "10_TssBiv", "11_BivFlnk") ~ "Promoter",
    state %in% c("3_TxFlnk", "4_Tx", "5_TxWk") ~ "Transcription",
    state %in% c("6_EnhG", "7_Enh", "12_EnhBiv") ~ "Enhancer",
    state %in% c("8_ZNF/Rpts", "9_Het", "13_ReprPC", "14_ReprPCWk", "15_Quies","8_ZNF_Rpts") ~ "Inactive chromatin",
    TRUE ~ "Other"
  ))



# Define the state descriptions and group order
state_descriptions <- c(
  "1_TssA" = "Active TSS",
  "2_TssAFlnk" = "Flanking Active TSS",
  "3_TxFlnk" = "Transcr. at gene 5' and 3'",
  "4_Tx" = "Strong transcription",
  "5_TxWk" = "Weak transcription",
  "6_EnhG" = "Genic enhancers",
  "7_Enh" = "Enhancers",
  "8_ZNF_Rpts" = "ZNF genes & repeats",
  "9_Het" = "Heterochromatin",
  "10_TssBiv" = "Bivalent/Poised TSS",
  "11_BivFlnk" = "Flanking Bivalent TSS/Enh",
  "12_EnhBiv" = "Bivalent Enhancer",
  "13_ReprPC" = "Repressed PolyComb",
  "14_ReprPCWk" = "Weak Repressed PolyComb",
  "15_Quies" = "Quiescent/Low"
)

custom_log_trans <- function(base = exp(1)) {
  trans <- function(x) ifelse(x > 0, log(x, base), x)
  inv <- function(x) ifelse(x > 0, base^x, x)
  scales::trans_new("custom_log", trans, inv)
}

# Preprocess the data to handle infinity
group_order <- c( "Transcription","Promoter", "Enhancer", "Inactive chromatin")


fnf_chromhmm_results_sorted <- fnf_chromhmm_results_addGroup %>%
  mutate(group = factor(group, levels = group_order),
         state = factor(state_descriptions[state], levels = state_descriptions)) %>%
  arrange(group, odds_ratio) %>%
  mutate(state = factor(state, levels = unique(state)))

pbs_chromhmm_results_sorted <- pbs_chromhmm_results_addGroup %>%
  mutate(group = factor(group, levels = group_order),
         state = factor(state_descriptions[state], levels = state_descriptions)) %>%
  arrange(group, odds_ratio) %>%
  mutate(state = factor(state, levels = fnf_chromhmm_results_sorted$state))



ggplot(pbs_chromhmm_results_sorted, aes(x= odd_med, y=pmin(-log10(epval))))


# Create a common theme
common_theme <- theme_minimal() +
  theme(
    panel.grid = element_blank(),
    strip.placement = "outside",
    panel.background = element_rect(fill = 'transparent', color = "transparent"),
    plot.background = element_rect(fill = 'transparent', color = "transparent"),
    text = element_text(family = "Helvetica"),
    axis.text.y = element_text(color = "black", size = 4),
    axis.title.y = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    axis.text.x = element_text(color = "black", size = 4),
    strip.text = element_text(face = "bold", hjust = 0),
    axis.line.x = element_line(linewidth = 0.25),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
    strip.background = element_blank(),
    axis.ticks = element_blank(),
    legend.direction = "horizontal",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "bottom"
  )

common_color_scale <- scale_color_gradient(
  low = "grey",
  high = "orange",
  name = "-log10(pval)",
  limits = c(0, 3),  # Set the limits from 0 to 3
  breaks = seq(0, 3, by = 1)  # Set breaks at 0, 1, 2, 3
)




gpbs_plot <- ggplot(pbs_chromhmm_results_sorted, aes(x = odds_ratio, y = state)) +
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = ci_upper, xmin = ci_lower), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, aes(color = pmin(-log10(epval), 3))) +  # Cap at 3
  facet_grid(group ~ ., scale = "free_y", switch = "y") +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4), labels = c(0, 1, 2, 3, 4)) +
  common_color_scale +
  common_theme +
  theme(strip.text.y.left = element_text(angle = 90, face = "bold", size = 8, vjust = 0.5, hjust = 0.5))+
  ylab("") +
  xlab("log(odds ratio)") +
  ggtitle("PBS")

# Create the FNF plot
fnf_plot <- ggplot(fnf_chromhmm_results_sorted, aes(x = odds_ratio, y = state)) +
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = ci_upper, xmin = ci_lower), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, aes(color = pmin(-log10(epval), 3))) +  # Cap at 3
  facet_grid(group ~ ., scale = "free_y", switch = "y") +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4), labels = c(0, 1, 2, 3, 4)) +
  common_color_scale +
  common_theme +
  theme(strip.text = element_blank())+
  theme(axis.text.y = element_blank(), strip.text.y = element_blank()) +
  ylab("") +
  xlab("log(odds ratio)") +
  ggtitle("FN-f")

# function to extract legend from plot 
get_only_legend <- function(plot) { 
  plot_table <- ggplot_gtable(ggplot_build(plot)) 
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box") 
  legend <- plot_table$grobs[[legend_plot]] 
  return(legend) 
} 

# extract legend from plot1 using above function 
legend <- get_only_legend(pbs_plot)


# Combine the plots
combined_plot <- grid.arrange(
  pbs_plot, fnf_plot, 
  ncol = 2, 
  widths = c(1.2, 1),
  top = textGrob("Enrichement sSNPs within functional genomic annotations", 
                 gp = gpar(fontface = "bold", fontsize = 10))
)

grid.arrange(combined_plot, legend, nrow = 2, heights = c(10, 1))



# pvalue
# Define colors
sig_color <- "orange"
nonsig_color <- "grey"

pbs_plot <- ggplot(pbs_chromhmm_results_sorted, aes(x = odds_ratio, y = state)) +
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = ci_upper, xmin = ci_lower), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, aes(color = epval < 0.05)) +
  facet_grid(group ~ ., scale = "free_y", switch = "y") +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4), labels = c(0, 1, 2, 3, 4)) +
  scale_color_manual(values = c("TRUE" = sig_color, "FALSE" = nonsig_color), 
                     name = "p-value < 0.05",
                     labels = c("TRUE" = "Yes", "FALSE" = "No")) +
  common_theme +
  theme(strip.text.y.left = element_text(angle = 90, face = "bold", size = 8, vjust = 0.5, hjust = 0.5)) +
  ylab("") +
  xlab("log(odds ratio)") +
  ggtitle("PBS")







#Try again




# Function to handle Inf and large values
handle_extreme_values <- function(x, max_value = 100) {
  ifelse(is.infinite(x) | abs(x) > max_value, NA, x)
}

# Prepare and combine the data
fnf_chromhmm_results_sorted <- fnf_chromhmm_results_addGroup %>%
  mutate(state = factor(state_descriptions[state], levels = state_descriptions),
         dataset = "FNF",
         odd_dnv = handle_extreme_values(as.numeric(odd_dnv)),
         odd_upv = handle_extreme_values(as.numeric(odd_upv)),
         odd_med = handle_extreme_values(as.numeric(odd_med)))

pbs_chromhmm_results_sorted <- pbs_chromhmm_results_addGroup %>%
  mutate(state = factor(state_descriptions[state], levels = state_descriptions),
         dataset = "PBS",
         odd_dnv = handle_extreme_values(as.numeric(odd_dnv)),
         odd_upv = handle_extreme_values(as.numeric(odd_upv)),
         odd_med = handle_extreme_values(as.numeric(odd_med)))

combined_data <- bind_rows(fnf_chromhmm_results_sorted, pbs_chromhmm_results_sorted) %>%
  mutate(state = factor(state, levels = rev(state_descriptions)))

# Filter out cases where observed < 5 and print them
low_observed <- combined_data %>%
  filter(observed < 5)
print("Cases where observed < 5 (these will be removed from the plot):")
print(low_observed[, c("state", "dataset", "observed", "odd_med", "odd_dnv", "odd_upv")])

# Remove low observed cases from the data for plotting
combined_data_filtered <- combined_data %>%
  filter(observed >= 5)

# Create a position dodge object
dodge <- position_dodge(width = 0.9)

# Create the combined plot
chromHMM_ORplot <- ggplot(combined_data_filtered, aes(x = odd_med, y = state, color = dataset, group = dataset)) +
  geom_vline(aes(xintercept = 1), linewidth = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmin = odd_dnv, xmax = odd_upv), 
                 size = 0.2, height = 0.4, position = dodge, color = "grey50") +
  geom_point(size = 2, position = dodge) +
  scale_x_continuous(breaks = c(0, 0.5,1, 1.5,2), labels = c(0, 0.5,1, 1.5,2)) +
  scale_color_manual(values = c("FNF" = "#FFB81C", "PBS" = "#005587"), name = "Dataset") +
  coord_cartesian(xlim = c(0,2))+
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    strip.placement = "outside",
    panel.background = element_rect(fill = 'transparent', color = "transparent"),
    plot.background = element_rect(fill = 'transparent', color = "transparent"),
    text = element_text(family = "Helvetica"),
    axis.text.y = element_text(color = "black", size = 4),
    axis.title.y = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    axis.text.x = element_text(color = "black", size = 4),
    strip.text = element_text(face = "bold", hjust = 0),
    axis.line.x = element_line(linewidth = 0.25),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
    strip.background = element_blank(),
    axis.ticks = element_blank(),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.5, "cm"),
    legend.position = c(0.9,0.05),
  ) +
  ylab("") +
  xlab("log(odds ratio)") +
  ggtitle("Chromatin State Enrichment Analysis")



#-------------------------------------------------------------------------------
# This is only for the eclip data from the encode

# Function to process a single file
process_file <- function(file_path) {
  # Extract tissue, RBP, and condition from filename
  #This section can be changed depends on if I have tissue or not. 
  file_info <- str_match(basename(file_path), "(.+)_(.+)_(.+)_sQTL_enrichment.txt")
  #file_info <- str_match(basename(file_path), "(.+)_(.+)_sQTL_enrichment.txt")
  tissue <- file_info[2]
  rbp <- file_info[3]
  condition <- file_info[4]
  
  # Read the file
  data <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)
  
  # Check if the file is empty or has unexpected content
  if (nrow(data) == 0 || ncol(data) != 8) {
    warning(paste("Skipping file due to unexpected format:", file_path))
    return(NULL)
  }
  
  # Rename columns
  colnames(data) <- c("observed", "total_qtls", "expected_mean", "expected_sd", "epval", "odd_ratio_dnv", "odd_ratio_med", "odd_ratio_upv")
  
  # Calculate fold change
  fold_change <- data$observed / data$expected_mean
  
  # Calculate odds ratio, confidence interval, and p-value
  fisher_result <- tryCatch({
    fisher.test(matrix(c(data$observed, data$total_qtls, round(data$expected_mean), data$total_qtls), ncol=2))
  }, error = function(e) {
    warning(paste("Error in Fisher's test for file:", file_path, "\nError:", e$message))
    return(list(conf.int = c(NA, NA), estimate = NA, p.value = NA))
  })
  
  
  # Create a data frame with all required information
  result <- data.frame(
    observed = data$observed,
    total_qtls = data$total_qtls,
    expected_mean = data$expected_mean,
    expected_sd = data$expected_sd,
    tissue = tissue,
    rbp = rbp,
    condition = condition,
    fold_change = fold_change,
    epval = data$epval,
    odd_dnv = data$odd_ratio_dnv,
    odd_med = data$odd_ratio_med,
    odd_upv= data$odd_ratio_upv,
    ci_lower = fisher_result$conf.int[1],
    ci_upper = fisher_result$conf.int[2],
    odds_ratio = fisher_result$estimate,
    fpval = fisher_result$p.value
  )
  
  return(result)
}

# Process all files in PBS directory
#pbs_files <- list.files("output/Enrichment/rbp/PBS", full.names = TRUE, pattern = "*.txt")
pbs_files <- list.files("output/Enrichment/eclip_tissue_sep/PBS", full.names = TRUE, pattern = "*.txt")
pbs_files_hepg2 <- pbs_files[grepl("HepG2", pbs_files)]
pbs_files_k562 <- pbs_files[grepl("K562", pbs_files)]
pbs_rbp_eclip_hepg2 <- map_df(pbs_files_hepg2, process_file)
pbs_rbp_eclip_k562 <- map_df(pbs_files_k562, process_file)
# Process all files in FNF directory
fnf_files <- list.files("output/Enrichment/eclip_tissue_sep/FNF", full.names = TRUE, pattern = "*.txt")
fnf_files_hepg2 <- fnf_files[grepl("HepG2", fnf_files)]
fnf_files_k562 <- fnf_files[grepl("K562", fnf_files)]
fnf_rbp_eclip_hepg2 <- map_df(fnf_files_hepg2, process_file)
fnf_rbp_eclip_k562 <- map_df(fnf_files_k562, process_file)



combined_data_hepg2 <- bind_rows(pbs_rbp_eclip_hepg2, fnf_rbp_eclip_hepg2)

combined_data_k562 <- bind_rows(pbs_rbp_eclip_k562, fnf_rbp_eclip_k562)

# Filter out cases where observed < 5 and print them
# low_observed <- combined_data %>%
#   filter(observed < 5)
# print("Cases where observed < 5 (these will be removed from the plot):")
# print(low_observed[, c("condition","observed", "odd_med", "odd_dnv", "odd_upv")])

# Remove low observed cases from the data for plotting
combined_data_hepg2_sig <- combined_data_hepg2 %>%
  group_by(rbp) %>%
  dplyr::filter(!all(epval > 0.05)) %>%
  ungroup()  |>
  arrange(desc(odd_med))


#combined_data_hepg2_sigFilltered <- combined_data_hepg2_sig |> filter(observed > 5) %>%
#  arrange(odd_med, epval)

combined_data_k562_sig <- combined_data_k562 %>%
  group_by(rbp) |>
  dplyr::filter(observed >5) |>
  dplyr::filter(!all(epval > 0.05)) |>
  ungroup() |>
  arrange(odd_med) 




#HEPG2 -Liver cancer cell line 
#K562 immortalized myelogenous leukemia cell line (Hematopoiesis, leukemia, cellular response) blood cell line


ggplot(combined_data_hepg2_sig, aes(x = odd_med, y = rbp, color = condition, group = condition)) +
  geom_vline(aes(xintercept = 1), linewidth = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmin = odd_dnv, xmax = odd_upv), 
                 size = 0.2, height = 0.4, position = position_dodge(width = 0.7), color = "grey50") +
  geom_point(size = 3, aes(alpha = ifelse(epval < 0.05, "< 0.05", "≥ 0.05")), 
             position = position_dodge(width = 0.7)) +
  scale_x_continuous(breaks = c(0,  1, 2, 3), labels = c(0,  1, 2, 3)) +
  scale_color_manual(values = c("fnf" = "#FFB81C", "pbs" = "#005587"), 
                     labels = c("pbs" = "PBS", "fnf" = "FN-f"),
                     name = "Group",
                     breaks = c("pbs", "fnf")) +
  scale_alpha_manual(values = c("< 0.05" = 1, "≥ 0.05" = 0.5), name = "p-value") +
  coord_cartesian(xlim = c(0, 3)) +
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
    plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
    strip.background = element_blank(),
    axis.ticks = element_blank(),
    legend.direction = "horizontal",
    legend.title = element_text(size=6),
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "bottom",
    legend.box.spacing = unit(0.1,"mm")
  ) +
  labs(y = "", x = "log(odds ratio)", 
       title = "RBP Enrichment Analysis-HePG2",
       color = "Group",
       alpha = "p-value") +
  guides(color = guide_legend(order = 1),
         alpha = guide_legend(order = 2))




# For K562
ggplot(combined_data_k562_sig, aes(x = odd_med, y = rbp, color = condition, group = condition)) +
  geom_vline(aes(xintercept = 1), linewidth = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmin = odd_dnv, xmax = odd_upv), 
                 size = 0.2, height = 0.4, position = position_dodge(width = 0.7), color = "grey50") +
  geom_point(size = 3, aes(alpha = ifelse(epval < 0.05, "< 0.05", "≥ 0.05")), 
             position = position_dodge(width = 0.7)) +
  scale_x_continuous(breaks = seq(0,5,1), labels = seq(0,5,1)) +
  scale_color_manual(values = c("fnf" = "#FFB81C", "pbs" = "#005587"), 
                     labels = c("pbs" = "PBS", "fnf" = "FN-f"),
                     name = "Group",
                     breaks = c("pbs", "fnf")) +
  scale_alpha_manual(values = c("< 0.05" = 1, "≥ 0.05" = 0.5), name = "p-value") +
  coord_cartesian(xlim = c(0, 5)) +
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
    plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
    strip.background = element_blank(),
    axis.ticks = element_blank(),
    legend.direction = "horizontal",
    legend.title = element_text(size=6),
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "bottom",
    legend.box.spacing = unit(0.1,"mm")
  ) +
  labs(y = "", x = "log(odds ratio)", 
       title = "RBP Enrichment Analysis- K562",
       color = "Group",
       alpha = "p-value") +
  guides(color = guide_legend(order = 1),
         alpha = guide_legend(order = 2))


# If I want to plot just all RBP not consider the cell-lines


# Function to process a single file
process_file <- function(file_path) {
  # Extract tissue, RBP, and condition from filename
  #This section can be changed depends on if I have tissue or not. 
  #file_info <- str_match(basename(file_path), "(.+)_(.+)_(.+)_sQTL_enrichment.txt")
  file_info <- str_match(basename(file_path), "(.+)_(.+)_sQTL_enrichment.txt")
  #tissue <- file_info[2]
  rbp <- file_info[2]
  condition <- file_info[3]
  
  # Read the file
  data <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)
  
  # Check if the file is empty or has unexpected content
  if (nrow(data) == 0 || ncol(data) != 8) {
    warning(paste("Skipping file due to unexpected format:", file_path))
    return(NULL)
  }
  
  # Rename columns
  colnames(data) <- c("observed", "total_qtls", "expected_mean", "expected_sd", "epval", "odd_ratio_dnv", "odd_ratio_med", "odd_ratio_upv")
  
  # Calculate fold change
  fold_change <- data$observed / data$expected_mean
  
  # Calculate odds ratio, confidence interval, and p-value
  fisher_result <- tryCatch({
    fisher.test(matrix(c(data$observed, data$total_qtls, round(data$expected_mean), data$total_qtls), ncol=2))
  }, error = function(e) {
    warning(paste("Error in Fisher's test for file:", file_path, "\nError:", e$message))
    return(list(conf.int = c(NA, NA), estimate = NA, p.value = NA))
  })
  
  
  # Create a data frame with all required information
  result <- data.frame(
    observed = data$observed,
    total_qtls = data$total_qtls,
    expected_mean = data$expected_mean,
    expected_sd = data$expected_sd,
    #tissue = tissue,
    rbp = rbp,
    condition = condition,
    fold_change = fold_change,
    epval = data$epval,
    odd_dnv = data$odd_ratio_dnv,
    odd_med = data$odd_ratio_med,
    odd_upv= data$odd_ratio_upv,
    ci_lower = fisher_result$conf.int[1],
    ci_upper = fisher_result$conf.int[2],
    odds_ratio = fisher_result$estimate,
    fpval = fisher_result$p.value
  )
  
  return(result)
}

# Process all files in PBS directory
#pbs_files <- list.files("output/Enrichment/rbp/PBS", full.names = TRUE, pattern = "*.txt")
pbs_files <- list.files("output/Enrichment/eclip_rbp/PBS", full.names = TRUE, pattern = "*.txt")
#pbs_files_hepg2 <- pbs_files[grepl("HepG2", pbs_files)]
#pbs_files <- pbs_files[grepl("K562|HepG2", pbs_files)]
pbs_rbp_eclip <- map_df(pbs_files, process_file)

# Process all files in FNF directory
fnf_files <- list.files("output/Enrichment/eclip_rbp/FNF", full.names = TRUE, pattern = "*.txt")
#fnf_files <- fnf_files[grepl("K562|HepG2", fnf_files)]
fnf_rbp_eclip <- map_df(fnf_files, process_file)

save(pbs_rbp_eclip, file="external_data/encode_rbp/rbp_prep_rbpOnly/pbs_rbp_elip.Rdata")
save(fnf_rbp_eclip, file="external_data/encode_rbp/rbp_prep_rbpOnly/fnf_rbp_elip.Rdata")

pbs_rbp_eclip$ratio_dnv_med <- pbs_rbp_eclip$odd_dnv / pbs_rbp_eclip$odd_med
fnf_rbp_eclip$ratio_dnv_med <- fnf_rbp_eclip$odd_dnv / fnf_rbp_eclip$odd_med

pbs_rbp_eclip$ratio_upv_med <- pbs_rbp_eclip$odd_upv / pbs_rbp_eclip$odd_med
fnf_rbp_eclip$ratio_upv_med <- fnf_rbp_eclip$odd_upv / fnf_rbp_eclip$odd_med

# Set thresholds (you can adjust these based on your data and requirements)
#threshold_dnv <- 1.5  # For example, odd_dnv shouldn't be more than 2.5 times odd_med
#threshold_upv <- 1.5  # For example, odd_upv shouldn't be more than 1.5 times odd_med

pbs_rbp_eclip_filtered <- pbs_rbp_eclip[which(pbs_rbp_eclip$ratio_dnv_med <= threshold_dnv & pbs_rbp_eclip$ratio_upv_med <= threshold_upv), ]
fnf_rbp_eclip_filtered <- fnf_rbp_eclip[which(fnf_rbp_eclip$ratio_dnv_med <= threshold_dnv & fnf_rbp_eclip$ratio_upv_med <= threshold_upv), ]


combined_data_rbp_eclip <- bind_rows(pbs_rbp_eclip_filtered, fnf_rbp_eclip_filtered)
combined_data_rbp_eclip <- bind_rows(pbs_rbp_eclip, fnf_rbp_eclip)

combined_data_rbp_eclip_sig <- combined_data_rbp_eclip %>%
  group_by(rbp) %>%
  dplyr::filter(observed >5) |>
  dplyr::filter(!all(epval > 0.05)) %>%
  ungroup()  |>
  arrange(desc(odd_med))


ggplot(combined_data_rbp_eclip_sig, aes(x = odd_med, y = rbp, color = condition, group = condition)) +
  geom_vline(aes(xintercept = 1), linewidth = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmin = odd_dnv, xmax = odd_upv), 
                 linewidth = 0.2, height = 0.4, position = position_dodge(width = 0.7), color = "grey50") +
  geom_point(size = 3, aes(alpha = ifelse(epval < 0.05, "< 0.05", "≥ 0.05")), 
             position = position_dodge(width = 0.7)) +
  scale_x_continuous(
    breaks = c(0, 0.5, 1, 1.5, 2),
    labels = c(0, 0.5, 1, 1.5, 2))+
  scale_color_manual(values = c("fnf" = "#FFB81C", "pbs" = "#005587"), 
                     labels = c("pbs" = "PBS", "fnf" = "FN-f"),
                     name = "Group",
                     breaks = c("pbs", "fnf")) +
  scale_alpha_manual(values = c("< 0.05" = 1, "≥ 0.05" = 0.3), name = "p-value") +
  coord_cartesian(xlim = c(0, 2)) +
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
    plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
    strip.background = element_blank(),
    axis.ticks = element_blank(),
    legend.direction = "horizontal",
    legend.title = element_text(size=6),
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "bottom",
    legend.box.spacing = unit(0.1,"mm")
  ) +
  labs(y = "", x = "log(odds ratio)", 
       title = "Enrichment of sSNPs in RNA binding protein sites",
       color = "Group",
       alpha = "p-value") +
  guides(color = guide_legend(order = 1),
         alpha = guide_legend(order = 2))





