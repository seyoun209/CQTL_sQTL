library(dplyr)
library(tidyr)
#library(formattable)
library(flextable)
library(officer)
library(openxlsx)

# This is the making Extended Data table 1 (Both chondrocytes and OA (16 samples) )

# Table 1 Characteristic 
load("output/combined_meta_data.RData") # It is loading name is combined_data
selected_columns <- c("Donor","ID","Condition","Sex", "Age","OAGradeAvg") #select column needed it based
meta_data <- combined_data %>% dplyr::select(all_of(selected_columns))
#Ancestry 
ancestry_df <- fread("/proj/phanstiel_lab/Data/processed/CQTL/geno/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_COA6_COA7/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_COA6_COA7_predictedAncestry.csv")
ancestry_df$Donor <- sub("^.*_(AM[0-9]+)_.*$", "\\1", ancestry_df$Donor)

ancestry_OA_df <- fread("/proj/phanstiel_lab/Data/processed/CQTL/geno/COA8_OA/ancestry/CQTL_COA8_predictedAncestry.csv")
ancestry_OA_df_filtered <- ancestry_OA_df %>%
  filter(grepl("^OA\\d+", Donor))
ancestry_OA_df_filtered$Donor <- gsub("_r2","",ancestry_OA_df_filtered$Donor)

ancestry_cqtl <-rbind(ancestry_df,ancestry_OA_df_filtered)
meta_cqtl <- merge(meta_data,ancestry_cqtl,by="Donor",all.x=TRUE)
unique_data <- meta_cqtl %>% distinct(Donor, Condition, .keep_all = TRUE)
## Summary statistics-----------------------------------------------------------
# Define age groups
unique_data <- unique_data %>%
  mutate(
    Age_Group = case_when(
      Age < 40 ~ "<40",
      Age >= 40 & Age < 60 ~ "40-59",
      Age >= 60 & Age < 80 ~ "60-79",
      Age >= 80 ~ "80+"
    )
  )

# Order the age groups manually
unique_data$Age_Group <- factor(unique_data$Age_Group, levels = c("<40", "40-59", "60-79", "80+"))

# Function to create summary tables for each condition
create_summary_tables <- function(data, condition) {
  data_condition <- data %>% filter(Condition == condition)
  
  # Age summary
  age_table <- data_condition %>%
    group_by(Age_Group) %>%
    summarise(
      n = n(),
      Age_mean = round(mean(Age), 1),
      Age_sd = round(sd(Age), 1)
    ) %>%
    mutate(
      Percentage = round(n / sum(n) * 100, 1),
      Age = paste0(Age_mean, " ± ", Age_sd)
    ) %>%
    select(Age_Group, n, Percentage, Age)
  
  # Sex summary
  sex_table <- data_condition %>%
    group_by(Sex) %>%
    summarise(
      n = n()
    ) %>%
    mutate(
      Percentage = round(n / sum(n) * 100, 1)
    )
  
  # Predicted Ancestry summary
  ancestry_table <- data_condition %>%
    group_by(Predicted_Ancestry) %>%
    summarise(
      n = n()
    ) %>%
    complete(Predicted_Ancestry = c("AMR", "EUR", "AFR", "EAS"), fill = list(n = 0)) %>%
    mutate(
      Percentage = round(n / sum(n) * 100, 1)
    )
  
  list(
    "Age (years)" = age_table,
    "Sex (Male/Female)" = sex_table,
    "Predicted Ancestry" = ancestry_table
  )
}

# Create summary tables for each condition
conditions <- unique(unique_data$Condition)
summary_tables <- lapply(conditions, function(cond) {
  create_summary_tables(unique_data, cond)
})

# Print the summary tables
for (i in seq_along(conditions)) {
  cat("\nCondition:", conditions[i], "\n")
  for (section in names(summary_tables[[i]])) {
    cat("\n", section, "\n")
    print(summary_tables[[i]][[section]])
  }
}





# Functions
### Function to create the flextable for a given condition----------------------

#This is for the flextable
# create_condition_flextable <- function(condition_data, condition_name, n_total) {
#   # Age table
#   age_data <- condition_data[["Age (years)"]] %>%
#     mutate(Characteristic = ifelse(row_number() == 1, "Age, n(%), Mean ± SD", ""),
#            Group = paste0("    ", Age_Group),
#            Value = paste0(n, " (", Percentage, "%)", ", ", Age)) %>%
#     select(Characteristic, Group, Value)
#   
#   # Sex table
#   sex_data <- condition_data[["Sex (Male/Female)"]] %>%
#     arrange(desc(n)) %>%
#     mutate(Sex = ifelse(Sex == "F", "Female", ifelse(Sex == "M", "Male", NA))) %>%
#     mutate(Characteristic = ifelse(row_number() == 1, "Sex, n(%)", ""),
#            Group = paste0("    ", Sex),
#            Value = paste0(n, " (", Percentage, "%)")) %>%
#     select(Characteristic, Group, Value)
#   
#   # Ancestry table
#   ancestry_data <- condition_data[["Predicted Ancestry"]] %>%
#     arrange(desc(n)) %>%
#     mutate(Characteristic = ifelse(row_number() == 1, "Predicted Ancestry, n(%)", ""),
#            Group = paste0("    ", Predicted_Ancestry),
#            Value = paste0(n, " (", Percentage, "%)")) %>%
#     select(Characteristic, Group, Value)
#   
#   # Combine all data and add spacing rows
#   all_data <- bind_rows(
#     age_data,
#     data.frame(Characteristic = "", Group = "", Value = ""),
#     sex_data,
#     data.frame(Characteristic = "", Group = "", Value = ""), # Spacing row
#     ancestry_data
#   )
#   
#   # Create flextable
#   ft <- flextable(all_data) %>%
#     set_header_labels(
#       Characteristic = paste(condition_name, "  (n =", n_total,")"),
#       Group = "",
#       Value = "Value"
#     ) %>%
#     merge_h(part = "header") %>%
#     align(align = "left", part = "all") %>%
#     align(j = 3, align = "center", part = "all") %>%
#     align(align = "center", part = "header") %>%
#     bold(j = 1, i = ~ Characteristic != "", bold = TRUE, part = "body") %>%
#     bold(bold = TRUE, part = "header") %>%
#     fontsize(size = 10, part = "all") %>%
#     padding(padding = 2, part = "all") %>%
#     border_remove() %>%
#     hline(i = 1, border = fp_border(color="black", width = 1), part = "header") %>%
#     height(height = 0.5, unit = "cm", part = "body") %>%
#     autofit()
#   
#   return(ft)
# }

# # Usage example:
# ctl_flextable <- create_condition_flextable(summary_tables[[1]], "Normal chondrocytes", 101)
# oa_flextable <- create_condition_flextable(summary_tables[[2]], "OA chondrocytes", 16)

create_condition_table <- function(condition_data, condition_name, n_total) {
  # Age table
  age_data <- condition_data[["Age (years)"]] %>%
    mutate(Characteristic = ifelse(row_number() == 1, "Age, n(%), Mean ± SD", ""),
           Group = paste0("    ", Age_Group),
           Value = paste0(n, " (", Percentage, "%)", ", ", Age)) %>%
    select(Characteristic, Group, Value)
  
  # Sex table
  sex_data <- condition_data[["Sex (Male/Female)"]] %>%
    arrange(desc(n)) %>%
    mutate(Sex = ifelse(Sex == "F", "Female", ifelse(Sex == "M", "Male", NA))) %>%
    mutate(Characteristic = ifelse(row_number() == 1, "Sex, n(%)", ""),
           Group = paste0("    ", Sex),
           Value = paste0(n, " (", Percentage, "%)")) %>%
    select(Characteristic, Group, Value)
  
  # Ancestry table
  ancestry_data <- condition_data[["Predicted Ancestry"]] %>%
    arrange(desc(n)) %>%
    mutate(Characteristic = ifelse(row_number() == 1, "Predicted Ancestry, n(%)", ""),
           Group = paste0("    ", Predicted_Ancestry),
           Value = paste0(n, " (", Percentage, "%)")) %>%
    select(Characteristic, Group, Value)
  
  # Combine all data and add spacing rows
  all_data <- bind_rows(
    data.frame(Characteristic = paste(condition_name, "n =", n_total), Group = "", Value = ""),
    data.frame(Characteristic = "", Group = "", Value = ""),
    age_data,
    data.frame(Characteristic = "", Group = "", Value = ""),
    sex_data,
    data.frame(Characteristic = "", Group = "", Value = ""),
    ancestry_data
  )
  
  return(all_data)
}

# Create tables for both conditions
ctl_table <- create_condition_table(summary_tables[[1]], "Normal chondrocytes", 101)
oa_table <- create_condition_table(summary_tables[[3]], "OA chondrocytes", 16)

# Create a new workbook
wb <- createWorkbook()

# Add worksheets
addWorksheet(wb, "Normal_Donor_characteristics")
addWorksheet(wb, "OA_Donor_characteristics")

# Write data to worksheets
writeData(wb, "Normal_Donor_characteristics", ctl_table, colNames = FALSE)
writeData(wb, "OA_Donor_characteristics", oa_table, colNames = FALSE)


# Apply some basic formatting
for (sheet in c("Normal_Donor_characteristics", "OA_Donor_characteristics")) {
  # Bold the title and category rows
  boldStyle <- createStyle(textDecoration = "bold")
  addStyle(wb, sheet, boldStyle, rows = c(1, 3, 8, 12), cols = 1)
  
  # Center the "Value" column
  centerStyle <- createStyle(halign = "center")
  addStyle(wb, sheet, centerStyle, rows = 1:nrow(ctl_table), cols = 3, gridExpand = TRUE)
  
  # Adjust column widths
  setColWidths(wb, sheet, cols = 1:3, widths = c(30, 20, 30))
}

# Save the workbook
saveWorkbook(wb, "output/results_plots/Supp_Data/Supplementary_Data1.xlsx", overwrite = TRUE)

























