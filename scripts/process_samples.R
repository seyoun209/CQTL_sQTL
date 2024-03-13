#Make the updated samplesheet with the cofounder.
#run code 
#Rscript scripts/process_samples.R ./aligned_samplesheet.txt "./donor_samples.txt" "./rna_extraction.txt" "CTL FNF" FragmentBatch RNAextractionKitBatch RIN 
#Rscript scripts/process_samples.R ./aligned_samplesheet.txt "./donor_samples.txt" "./rna_extraction.txt" "CTL FNF" FragmentBatch RNAextractionKitBatch RIN wasp
#Rscript scripts/process_samples.R ./aligned_samplesheet.txt "./donor_samples.txt" "./rna_extraction.txt" "CTL OA" FragmentBatch RNAextractionKitBatch RIN 
#Rscript scripts/process_samples.R ./aligned_samplesheet.txt "./donor_samples.txt" "./rna_extraction.txt" "CTL OA" FragmentBatch RNAextractionKitBatch RIN wasp
# process_samples.R
args <- commandArgs(trailingOnly = TRUE)
aligned_samplesheet_path <- args[1]
donor_samples_path <- args[2]
rna_extraction_path <- args[3]
conditions_to_include <- unlist(strsplit(args[4], " ")) # Split the first argument into condition values
use_wasp_id <- tail(args, 1) == "wasp" # Check if the last argument is "wasp"

# Correctly adjust additional_confounds based on the presence of "wasp"
if(use_wasp_id) {
  additional_confounds <- args[5:(length(args) - 1)] # Exclude last argument if it's "wasp"
} else {
  additional_confounds <- args[5:length(args)] # Include all remaining arguments otherwise
}

library(dplyr)
library(yaml)
library(readr)
library(data.table)

# Load configuration from YAML file
config <- yaml::read_yaml("./config/rna_prcoess.yaml")
aligned_samples <- fread(aligned_samplesheet_path)
donor_samples <- fread(donor_samples_path)
rna_extraction <- fread(rna_extraction_path)

# Identify common columns between aligned_samples and donor_samples, excluding 'Donor'
common_cols_donor <- setdiff(intersect(colnames(aligned_samples), colnames(donor_samples)), "Donor")

# Remove these common columns from donor_samples before the join
donor_samples_clean <- dplyr::select(donor_samples, -all_of(common_cols_donor))

# Perform the first left join
combined_data <- aligned_samples %>%
  left_join(donor_samples_clean, by = "Donor")

common_cols_rna <- setdiff(intersect(colnames(combined_data), colnames(rna_extraction)), "Read2")
rna_extraction_clean <- dplyr::select(rna_extraction, -all_of(common_cols_rna))
combined_data <- combined_data %>%
  left_join(rna_extraction_clean, by = "Read2")

# Process ID column with and without "_wasp"
combined_data <- combined_data %>%
  mutate(ID = paste(Donor, Condition, Tech_Rep, Sex, sep = "_"),
         ID_wasp = paste(Donor, Condition, Tech_Rep, Sex, "wasp", sep = "_"))

# Conditionally adjust FragmentBatch for 'OA' condition
combined_data <- combined_data %>%
  mutate(FragmentBatch = ifelse(Condition == "OA", 0, FragmentBatch))

# Omit samples and filter rows based on configuration and conditions
combined_data <- combined_data[!combined_data$Donor %in% config$samples_to_omit, ] %>%
  filter(Condition %in% conditions_to_include)

# Select and rename columns dynamically, based on the use of wasp ID
id_column <- if(use_wasp_id) "ID_wasp" else "ID"
selected_columns <- c(id_column, "Condition", additional_confounds)
final_data <- combined_data %>% dplyr::select(all_of(selected_columns))

# Writing the output
output_prefix <- if(use_wasp_id) "wasp_" else ""
output_file <- sprintf("output/%s%s_group.txt", output_prefix, paste(conditions_to_include, collapse="vs"))

write.table(final_data, output_file, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
