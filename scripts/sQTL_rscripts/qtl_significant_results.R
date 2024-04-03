# Finding the total counts and find the significant QTL
## Author: Seyoun Byun
## Date: 03.08.2024
## Edited:
#-------------------------------------------------------------------------------
setwd("/work/users/s/e/seyoun/CQTL_sQTL")
library(tidyr)
library(dplyr)
library(ggplot2)


# Define the possible arguments files, significant p-values etc
adjusted_beta_p <- 0.05
pbs_perm_dir <- "output/01.qtltools_re/perm_pbs/"
fnf_perm_dir <- "output/01.qtltools_re/perm_fnf/"
pbs_wasp_perm_dir <- "output/01.qtltools_re/perm_pbs_wasp/"
fnf_wasp_perm_dir <- "output/01.qtltools_re/perm_fnf_wasp/"


# Function for reading and processing PBS files --------------------------------
process_pbs_files <- function(filepath_pattern, pattern_suffix) {
filepaths <- list.files(filepath_pattern, pattern = pattern_suffix, full.names = TRUE)   # List files based on the pattern

# Sort filepaths numerically based on the PC number
sorted_filepaths <- filepaths[order(as.numeric(gsub(".*pc(\\d+)_.*", "\\1", filepaths)))]

# Read files and process them
list_processed <- lapply(sorted_filepaths, function(filepath) {
  df <- read.table(filepath, header = FALSE, stringsAsFactors = FALSE)
  
# Filter out rows with "chrX" or "chrY" from the specified column
df_filtered <- df[!(df[, 2] %in% c("chrX", "chrY")), ]
    
# Add adjusted p-values using FDR
df_filtered$p_adjusted <- p.adjust(df_filtered[, 20], method = "fdr")
    
# Filter based on significance threshold
df_sig <- df_filtered[df_filtered$p_adjusted <= adjusted_beta_p, ]

return(df_sig)
})
  
  return(list_processed)
}


pbs_list_processed <- process_pbs_files(pbs_perm_dir, "_allchr.pbs.perm")
fnf_list_processed <- process_pbs_files(fnf_perm_dir, "_allchr.fnf.perm")

pbs_wasp_list_processed <- process_pbs_files(pbs_wasp_perm_dir, "allchr.pbs.wasp.perm")
fnf_wasp_list_processed <- process_pbs_files(fnf_wasp_perm_dir, "allchr.fnf.perm")



#dotplot to see the counts of snps only for the ---------------------------------

# Count unique values in column 1 of each list element
pbs_unique_counts <- lapply(pbs_list_processed, function(x) {
  length(unique(x$V1))
})

fnf_unique_counts <- lapply(fnf_list_processed, function(x) {
  length(unique(x$V1))
})


# Create a data frame with the list elements and names
df <- data.frame(pc = paste0("pc", 1:20), pbs_counts = unlist(pbs_unique_counts), fnf_counts=unlist(fnf_unique_counts))
df$pc <- factor(df$pc, levels = paste0("pc", 1:20))
df_long <- df %>%
  pivot_longer(cols = c(pbs_counts, fnf_counts),
               names_to = "type",
               values_to = "counts") %>%
  mutate(type = recode(type, pbs_counts = "PBS", fnf_counts = "FNF"))

df_long$type <- factor(df_long$type, levels = c("PBS","FNF"))

# Calculate the max for each group
max_counts <- df_long %>%
  group_by(type) %>%
  summarize(max_count = max(counts), pc_max = pc[which.max(counts)]) %>%
  ungroup()

max_counts$color <- ifelse(max_counts$type == "PBS", "#2ea7e8", "#f77943")


pdf(file = "output/results_plots/Significant_count_introns.pdf",   # The directory you want to save the file in
    width = 11, # The width of the plot in inches
    height = 6.5) # The height of the plot in inches

# Create the dot plot
ggplot(df_long, aes(x = pc, y = counts, color = type)) +
  geom_point(size = 3) +
  labs(x = "PC", y = "Number of significant introns", color = "Type") +
  scale_color_manual(values = c("PBS" = "#9FCCE4", "FNF" = "#FAB394")) +
  scale_y_continuous(limits = c(2500, 5000)) +
  theme_classic() +
  # geom_segment(data = max_counts, aes(x = "pc6", xend = pc_max, y = max_count +150, yend = max_count), 
  #              arrow = arrow(type = "closed", length = unit(0.2, "cm")), color = "black") +
  geom_text(data = max_counts, aes(x = pc_max, y = max_count + 100, label = max_count), 
            color = max_counts$color, vjust = 0, size = 3.5)
dev.off()


#dotplot to see WASP counts------------------- ---------------------------------

# Count unique values in column 1 of each list element
pbs_wasp_counts <- lapply(pbs_wasp_list_processed, function(x) {
  length(unique(x$V1))
})

fnf_wasp_counts <- lapply(fnf_wasp_list_processed, function(x) {
  length(unique(x$V1))
})


# Create a data frame with the list elements and names
df_wasp <- data.frame(pc = paste0("pc", 1:20), pbs_counts = unlist(pbs_wasp_counts), fnf_counts=unlist(fnf_wasp_counts))
df_wasp$pc <- factor(df_wasp$pc, levels = paste0("pc", 1:20))
df_wasp_long <- df_wasp %>%
  pivot_longer(cols = c(pbs_counts, fnf_counts),
               names_to = "type",
               values_to = "counts") %>%
  mutate(type = recode(type, pbs_counts = "PBS", fnf_counts = "FNF"))

df_wasp_long$type <- factor(df_wasp_long$type, levels = c("PBS","FNF"))

# Calculate the max for each group
max_counts_wasp <-df_wasp_long %>%
  group_by(type) %>%
  summarize(max_count = max(counts), pc_max = pc[which.max(counts)]) %>%
  ungroup()

max_counts_wasp$color <- ifelse(max_counts_wasp$type == "PBS", "#2ea7e8", "#f77943")

pdf(file = "output/results_plots/wasp_Significant_count_introns.pdf",   # The directory you want to save the file in
    width = 11, # The width of the plot in inches
    height = 6.5)
# Create the dot plot
ggplot(df_wasp_long, aes(x = pc, y = counts, color = type)) +
  geom_point(size = 3) +
  labs(x = "PC", y = "Number of significant introns", color = "Type") +
  scale_color_manual(values = c("PBS" = "#9FCCE4", "FNF" = "#FAB394")) +
  scale_y_continuous(limits = c(300, 700)) +
  theme_classic() +
  # geom_segment(data = max_counts, aes(x = "pc6", xend = pc_max, y = max_count +150, yend = max_count), 
  #              arrow = arrow(type = "closed", length = unit(0.2, "cm")), color = "black") +
  geom_text(data = max_counts_wasp, aes(x = pc_max, y = max_count + 10, label = max_count), 
            color = max_counts_wasp$color, vjust = 0, size = 3.5)


dev.off()


#check pvalue are okay ---------------------------------------------------------
perm_pbs_pc5 <- read.table(paste0(pbs_perm_dir,"pc5_allchr.pbs.perm"), header = FALSE, stringsAsFactors = FALSE)
perm_fnf_pc4 <- read.table(paste0(fnf_perm_dir,"pc4_allchr.fnf.perm"), header = FALSE, stringsAsFactors = FALSE)

## check first beta approximated P-values are okay

pdf(file = "output/results_plots/beta_approximation_p_value_check_Both_cond.pdf",   # The directory you want to save the file in
    width = 11, # The width of the plot in inches
    height = 6.5) # The height of the plot in inches

par(mfrow=c(1,2))
plot(perm_pbs_pc5$V19, perm_pbs_pc5$V20, xlab="Direct method", ylab="Beta approximation", main="PBS-Beta approximation check (PC5)")
abline(0, 1, col="red")

plot(perm_fnf_pc4$V19, perm_fnf_pc4$V20, 
     xlab="Direct method", 
     ylab="Beta approximation", 
     main="FNF-Beta approximation check (PC4)")
abline(0, 1, col="red") # Add a red line to this plot as well
par(mfrow=c(1,1))

dev.off()


#-------------------------------------------------------------------------------
#Finding the significant 


