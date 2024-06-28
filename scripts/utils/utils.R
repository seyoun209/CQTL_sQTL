library(dplyr)
library(data.table)
library(rrvgo)
# Function to reduce similar GO terms
reduceGO <- function(goterms, category, ont = "BP", threshold = 0.8){
  
  # Calculate similarity matrix for GO terms based on ontology
  simMatrix <- calculateSimMatrix(goterms$TermID,
                                  orgdb = "org.Hs.eg.db",
                                  ont = ont,
                                  method = "Rel")
  
  # Create named vector of scores
  # Higher is better -> -log10 transforming p-values
  scores <- setNames(-log10(goterms$pval), goterms$TermID)
  
  # Group GO terms based on similarity threshold
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold = threshold,
                                  orgdb = "org.Hs.eg.db") |>
    dplyr::rename(`-log10pval` = score)
  
  # Join grouped terms with original data
  joined_reducedTerms <- left_join(goterms, reducedTerms,
                                   by = join_by("TermID" == "go", "Term" == "term")) |>
    dplyr::select(-cluster) |>
    relocate(parent, .after = Term) |>
    relocate(parentTerm, .after = parent) |>
    dplyr::rename(parentTermID = parent) |>
    mutate("category" = category)
  
  return(joined_reducedTerms)
}


#This is the normalize the counts to psi
normalize_column <- function(x) {
  x / sum(x, na.rm = TRUE)
}



#Prep for the PCA plot
pca_prep <- function(data_matrix, metadata,remove_vars) {
  # Transpose the data matrix
  data_transposed <- t(data_matrix)
  
  # Convert to data frame and numeric
  data_df <- as.data.frame(data_transposed)
  data_df[] <- lapply(data_df, as.numeric)
  
  # Merge with metadata
  merged_data <- merge(data_df, metadata, by.x = 'row.names', by.y = "ID", all = TRUE)
  
  # Remove unnecessary columns
  merged_data_num <- merged_data[, !(colnames(merged_data) %in% remove_vars)]
  zero_variance_columns <- which(apply(merged_data_num, 2, var) == 0)
  if( length(zero_variance_columns) != 0 ){
    merged_data_num <- merged_data_num[, -zero_variance_columns]
  }

  # Perform PCA
  pca_obj <- prcomp(merged_data_num, scale. = TRUE)
  pca_df <- as.data.frame(pca_obj$x)
  pca_df$Condition <- merged_data$Condition
  pca_df$Donor <- merged_data$Donor
  
  # Calculate variance explained
  variance_explained <- summary(pca_obj)$importance[2, ]
  
  return(list(pca_data = pca_df, variance_explained = variance_explained))
}

# Calculate deltaPSI from the condition
calculate_delta_psi <- function(data, pattern_control, pattern_case) {
  data$deltaPSI <- apply(data, 1, function(row) {
    ctl_cols <- grep(pattern_control, colnames(data))
    case_cols <- grep(pattern_case, colnames(data))
    
    ctl_mean <- mean(as.numeric(row[ctl_cols]), na.rm = TRUE)
    case_mean <- mean(as.numeric(row[case_cols]), na.rm = TRUE)
    
    if (is.na(ctl_mean) || is.na(case_mean)) {
      NA
    } else {
      case_mean - ctl_mean
    }
  })
  return(data)
  }



join_introns_deltapsi_fdr <- function(limma_data, introns_data, cluster_sig_file) {
  psi_df_subset <- limma_data %>%
    rownames_to_column("phe_id") %>%
    dplyr::select(phe_id, deltaPSI)
  
  introns_data$phe_id <- paste(introns_data$chr,
                               introns_data$start,
                               introns_data$end,
                               introns_data$clusterID,
                               sep = ":")
  
  introns_filtered_unknown <- introns_data %>%
    dplyr::filter(verdict != "unknown_strand")
  
  introns_filtered_subset <- semi_join(introns_filtered_unknown, psi_df_subset, by = "phe_id")
  introns_new <- left_join(introns_filtered_subset, psi_df_subset, by = "phe_id")
  introns_new <- introns_new %>%
    dplyr::rename(deltapsi_batch = deltaPSI)
  
  cluster_sig <- fread(cluster_sig_file)
  cluster_sig_ID <- cluster_sig %>%
    dplyr::mutate(chr = str_split(cluster_sig$cluster, ":", simplify = TRUE)[, 1],
                  clusterID = str_split(cluster_sig$cluster, ":", simplify = TRUE)[, 2]) %>%
    dplyr::select(clusterID, p, p.adjust, loglr, genes)
  
  introns_pval_include <- left_join(introns_new, cluster_sig_ID, by = "clusterID")
  
  return(introns_pval_include)
}




#gene log2FC
get_sample_l2fc <- function(gene, countMatrix){
  
  # Extract row of gene from countMatrix
  gene_counts <- countMatrix[gene,]
  
  # Convert to dataframe and extract donors/conditions into separate columns
  donor_gene_counts <- data.frame(gene_counts) %>%
    rownames_to_column(var = "Sample") %>%
    separate_wider_delim("Sample", 
                         delim = "_", 
                         names = c("Donor", "Condition", NA, "Sex")) %>%
    # Group by each donor and calculate l2FC 
    group_by(Donor) %>%
    reframe(log2FC = 
                log2(gene_counts[Condition == "FNF"]/gene_counts[Condition == "CTL"])) %>%
    ungroup() %>%
    mutate(ENSEMBL = gene)
  
  return(donor_gene_counts)
}


get_gene_condition_Counts <- function(gene, dds){
  geneCounts <- DESeq2::plotCounts(dds, gene = gene,
                                   intgroup = "Condition",
                                   normalized = TRUE,
                                   returnData = TRUE) |>
    remove_rownames() |>
    mutate(gene_id = gene)
  
  return(geneCounts)
}


#Fraction calculate -Coloc

calculate_case_fraction <- function(case_control_sizes) {
  # Check if the required columns exist
  if (!all(c("Max_Cases", "Max_Controls") %in% colnames(case_control_sizes))) {
    stop("The data frame must contain 'Max_Cases' and 'Max_Controls' columns.")
  }
  
  # Check if there's at least one row of data
  if (nrow(case_control_sizes) == 0) {
    warning("The data frame is empty.")
    return(NA)
  }
  
  # Calculate the fraction
  cases <- case_control_sizes$Max_Cases[1]
  controls <- case_control_sizes$Max_Controls[1]
  total <- cases + controls
  
  if (total == 0) {
    warning("Total number of cases and controls is zero.")
    return(NA)
  }
  
  fraction <- cases / total
  
  # Return a named list with additional information
  return(list(
    fraction_cases = fraction,
    total_samples = total,
    number_cases = cases,
    number_controls = controls
  ))
}



# Function to process var_id
process_var_id <- function(var_id) {
  # Remove "chr" prefix
  var_id_no_chr <- str_remove(var_id, "^chr")
  
  # Split by ":" and keep only first two elements
  split_elements <- str_split(var_id_no_chr, ":", n = 3, simplify = TRUE)
  
  # Combine the first two elements with a colon
  processed_id <- paste(split_elements[,1], split_elements[,2], sep = ":")
  
  return(processed_id)
}
