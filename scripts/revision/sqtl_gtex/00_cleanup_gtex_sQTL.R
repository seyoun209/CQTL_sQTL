setwd("/work/users/s/e/seyoun/CQTL_sQTL/")
# Load libraries
library(dplyr)
library(arrow)
library(data.table)
library(ggplot2)
library(httpgd)


gtex_dir <- "external_data/GTEx_v10_sQTL/GTEx_Analysis_v10_sQTL_updated"
#------------------------------------------------------------------------------------
# Load the datas
## load the sQTL results
response_pbs_results <- readRDS("/work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/conditional_pbs/response_pbs_resInv.rds")
response_fnf_results <- readRDS("/work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/conditional_pbs/response_fnf_resInv.rds")

#------------------------------------------------------------------------------------

## Gtex data
# ThiS is the passing from the nominal threshold 
tissue_files <- list.files(path = gtex_dir, pattern = "\\.sQTLs\\.signif_pairs\\.parquet$", full.names = TRUE)
tissues <- gsub("\\.v10\\.sQTLs\\.signif_pairs\\.parquet$", "", basename(tissue_files))

#read_parquet(tissue_files[1]) %>%
#  mutate(tissue = tissues[1]) %>%
#  as.matrix() |> head()


# GTEx sQTLs for the sGenes

#-------------------------------------------------------------------------------------
# Make a function to process GTEx sQTL files and the sQTL LD files. 


process_gtex_sqtl_files <- function(gtex_dir, tissues = NULL) {
  # List all GTEx sQTL files, or use provided tissue list
  if (is.null(tissues)) {
    files <- list.files(gtex_dir, pattern = "*.sGenes.txt.gz", full.names = TRUE)
    tissues <- gsub(".v10.sGenes.txt.gz", "", basename(files))
  } else {
    files <- file.path(gtex_dir, paste0(tissues, ".v10.sGenes.txt.gz"))
  }
  
  # Initialize list to store processed data frames
  tissue_sqtls <- list()
  
  # Process each tissue file
  for (i in seq_along(files)) {
    print(i)
    tissue <- tissues[i]
    file <- files[i]
    
    # Check if file exists
    if (!file.exists(file)) {
      warning(paste("File not found:", file))
      next
    }
    
    # Read GTEx sQTL data
    gtex_data <- fread(file, sep = "\t", header = TRUE)
    
    # Extract relevant columns and standardize format
    gtex_processed <- gtex_data %>%
      select(phenotype_id, gene_id, variant_id, pval_nominal, qval) %>%
      # Parse the variant_id to match your format
      mutate(
        # Convert variant_id from "chr1_17005_A_G_b38" to "chr1:17005:A:G"
        var_id = gsub("_b38$", "", variant_id),
        var_id = gsub("_", ":", var_id),
        # Extract ENSG without version number
        ensg = gsub("\\.\\d+$", "", gene_id),
        # Extract intron/cluster info from phenotype_id 
        phe_id = sub("^(chr[^:]+:[^:]+:[^:]+).*", "\\1", phenotype_id)
      )
    
    tissue_sqtls[[tissue]] <- gtex_processed
  }
  
  return(tissue_sqtls)
}

sQTL_gtex_list <- process_gtex_sqtl_files(gtex_dir, tissues = NULL)

#-------------------------------------------------------------------------------------
# Comparison sQTLs and GTEx sQTLs
#-------------------------------------------------------------------------------------

compare_ld_variants_with_gtex <- function(chon_sqtls, gtex_sqtls, ld_dir, ld_threshold = 0.7) {
  # Initialize results tracking with additional columns
  results <- data.frame(
    var_id = character(),
    phe_id = character(),    # Added chondrocyte phenotype ID
    ensg = character(),      # Added chondrocyte ENSG
    symbol = character(),    # Added gene symbol if available
    tissue = character(),
    ld_variant = character(),
    gtex_match = logical(),
    gtex_gene = character(),
    gtex_phe_id = character(), # Added GTEx phenotype ID
    gtex_pval = numeric(),
    ld_r2 = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Process each sQTL in your data
  for (i in seq_len(nrow(chon_sqtls))) {
    if (i %% 100 == 0) cat("Processing sQTL", i, "of", nrow(chon_sqtls), "\n")
    
    # Extract key information
    var_id <- chon_sqtls$var_id[i]
    phe_id <- chon_sqtls$genomicLoc[i]        # Use genomicLoc directly
    ensg <- chon_sqtls$ensg[i]            # Use ensg directly
    symbol <- if ("SYMBOL" %in% colnames(chon_sqtls)) chon_sqtls$SYMBOL[i] else NA_character_
    
    # Extract chromosome from var_id
    chr_parts <- strsplit(var_id, ":")[[1]]
    if (length(chr_parts) < 1) {
      warning("Invalid var_id format for row ", i, ": ", var_id)
      next
    }
    chr_sub <- chr_parts[1]
    ld_file <- file.path(ld_dir, chr_sub, paste0(var_id, ".ld"))
    
    # Check if the LD file exists
    if (!file.exists(ld_file)) {
      cat("No LD file found for variant:", var_id, "\n")
      next
    }
    
    # Load and filter LD data
    ld_data <- fread(ld_file)
    ld_hi <- ld_data[R2 > ld_threshold]
    
    if (nrow(ld_hi) == 0) {
      cat("No variants in LD (R² >", ld_threshold, ") found for:", var_id, "\n")
      next
    }
    
    # Get unique LD variants
    ld_variants <- unique(ld_hi$SNP_B)
    
    # Check each GTEx tissue
    for (tissue_name in names(gtex_sqtls)) {
      gtex_tissue_data <- gtex_sqtls[[tissue_name]]
      
      # Get all variants in this tissue
      gtex_variants <- unique(gtex_tissue_data$var_id)
      
      # Find intersection between LD variants and GTEx variants
      common_vars <- intersect(ld_variants, gtex_variants)
      
      # If matches found, record details
      if (length(common_vars) > 0) {
        for (ld_var in common_vars) {
          # Get LD value
          ld_value <- ld_hi[SNP_B == ld_var, R2][1]
          
          # Get GTEx information for this variant
          gtex_matches <- gtex_tissue_data[gtex_tissue_data$var_id == ld_var, ]
          
          for (j in seq_len(nrow(gtex_matches))) {
            results <- rbind(results, data.frame(
              var_id = var_id,
              phe_id = phe_id,       # Include chondrocyte phenotype ID
              ensg = ensg,           # Include chondrocyte ENSG
              symbol = symbol,       # Include gene symbol
              tissue = tissue_name,
              ld_variant = ld_var,
              gtex_match = TRUE,
              gtex_gene = gtex_matches$ensg[j],
              gtex_phe_id = gtex_matches$phe_id[j], # Include GTEx phenotype ID
              gtex_pval = gtex_matches$pval_nominal[j],
              ld_r2 = ld_value,
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
  }
  
  # Summarize results
  if (nrow(results) > 0) {
    # Count unique matches by tissue
    tissue_summary <- results %>%
      group_by(tissue) %>%
      summarise(
        unique_sQTLs = n_distinct(var_id),
        unique_LD_variants = n_distinct(ld_variant),
        total_matches = n()
      )
    
    # Count variants that match in multiple tissues
    variant_tissue_counts <- results %>%
      group_by(var_id, phe_id, ensg) %>%  # Group by all identifiers
      summarise(
        num_tissues = n_distinct(tissue),
        tissues = paste(sort(unique(tissue)), collapse=",")
      )
    
    # Print summary statistics
    cat("\nSummary of LD variant matches in GTEx tissues:\n")
    print(tissue_summary)
    
    cat("\nVariants with matches in multiple tissues:\n")
    print(table(variant_tissue_counts$num_tissues))
    
    return(list(
      all_matches = results,
      tissue_summary = tissue_summary,
      variant_tissue_counts = variant_tissue_counts
    ))
  } else {
    cat("No LD variant matches found in any GTEx tissue.\n")
    return(NULL)
  }
}
#your_sqtls, gtex_sqtls, ld_dir, ld_threshold = 0.7

pbs_comparison_gtex <- compare_ld_variants_with_gtex(chon_sqtls =response_pbs_results, 
gtex_sqtls = sQTL_gtex_list,
ld_dir = "/work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld",
ld_threshold = 0.7)

fnf_comparison_gtex <- compare_ld_variants_with_gtex(chon_sqtls =response_fnf_results, 
gtex_sqtls = sQTL_gtex_list,
ld_dir = "/work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld",
ld_threshold = 0.7)


save(pbs_comparison_gtex,fnf_comparison_gtex, sQTL_gtex_list, 
     file = "/work/users/s/e/seyoun/CQTL_sQTL/output/revision/data/comparison_gtex.rdata")

load("/work/users/s/e/seyoun/CQTL_sQTL/output/revision/data/comparison_gtex.rdata")
#-------------------------------------------------------------------------------------
#Summary Table sQTLs and sQTL GTEx
#-------------------------------------------------------------------------------------
# Makeing a table
pbs_gtex_summary_filtered <- pbs_comparison_gtex$all_matches %>%
  filter(phe_id == gtex_phe_id) %>% 
  # Group by chondrocyte sQTL (var_id, phe_id) and by each tissue separately:
  group_by(var_id, phe_id, tissue) %>%
  summarize(
    ld_info = paste(unique(paste0(ld_variant, " (R2=", ld_r2, ")")), collapse = "; "),
    .groups = "drop"
  ) %>%
  # Now, for each sQTL group, aggregate the tissues and their ld_info together.
  group_by(var_id, phe_id) %>%
  summarize(
    tissues = paste(sort(unique(tissue)), collapse = ", "),
    tissue_count = n_distinct(tissue),
    ld_info = paste(paste0(tissue, ": ", ld_info), collapse = " | "),
    .groups = "drop"
  ) %>%
  # Rename phe_id to genomicLoc for the join
  rename(genomicLoc = phe_id)
# Merge the summarized GTEx info back to your chondrocyte sQTL results.
pbs_final_table <- response_pbs_results %>%
  left_join(pbs_gtex_summary_filtered, by = c("var_id", "genomicLoc"))

gtex_exist_pbs <- pbs_final_table |> filter(!is.na(tissue_count)) 


fnf_gtex_summary_filtered <- fnf_comparison_gtex$all_matches %>%
  filter(phe_id == gtex_phe_id) %>% 
  # Group by chondrocyte sQTL (var_id, phe_id) and by each tissue separately:
  group_by(var_id, phe_id, tissue) %>%
  summarize(
    ld_info = paste(unique(paste0(ld_variant, " (R2=", ld_r2, ")")), collapse = "; "),
    .groups = "drop"
  ) %>%
  # Now, for each sQTL group, aggregate the tissues and their ld_info together.
  group_by(var_id, phe_id) %>%
  summarize(
    tissues = paste(sort(unique(tissue)), collapse = ", "),
    tissue_count = n_distinct(tissue),
    ld_info = paste(paste0(tissue, ": ", ld_info), collapse = " | "),
    .groups = "drop"
  ) %>%
  # Rename phe_id to genomicLoc for the join
  rename(genomicLoc = phe_id)
# Merge the summarized GTEx info back to your chondrocyte sQTL results.
fnf_final_table <- response_fnf_results %>%
  left_join(fnf_gtex_summary_filtered, by = c("var_id", "genomicLoc"))

gtex_exist_fnf <- fnf_final_table |> filter(!is.na(tissue_count)) 

#-------------------------------------------------------------------------------------
library(tidyr)
library(ComplexHeatmap)

# Create a unique identifier for each chondrocyte sQTL (using var_id and genomicLoc)
heatmap_df <- pbs_comparison_gtex$all_matches %>%
  filter(phe_id == gtex_phe_id, ld_r2 > 0.7) %>% 
  mutate(sqtl_id = paste(var_id, phe_id, sep = "|")) %>%
  select(sqtl_id, tissue) %>%
  distinct() %>%
  mutate(presence = 1)

# Pivot the data to create a matrix where rows = unique sQTL IDs and columns = tissues
heatmap_matrix <- heatmap_df %>%
  pivot_wider(names_from = tissue, values_from = presence, values_fill = list(presence = 0)) %>%
  as.data.frame()

# Set row names and drop the sqtl_id column
rownames(heatmap_matrix) <- heatmap_matrix$sqtl_id
heatmap_matrix$sqtl_id <- NULL

# Convert to matrix
heatmap_matrix <- as.matrix(heatmap_matrix)

# Optionally limit the number of tissues (columns) to 50 if too many exist:
if(ncol(heatmap_matrix) > 50) {
  heatmap_matrix <- heatmap_matrix[, 1:50]
}

# Plot the heatmap with ComplexHeatmap
Heatmap(heatmap_matrix,
        name = "GTEx Hit",
        col = c("0" = "white", "1" = "darkblue"),
        show_row_names = TRUE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_title = "Chondrocyte sQTL (var_id|genomicLoc)",
        column_title = "GTEx Tissues",
        heatmap_legend_param = list(title = "Presence"))
