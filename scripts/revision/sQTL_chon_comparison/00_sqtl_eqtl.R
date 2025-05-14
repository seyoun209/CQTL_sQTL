setwd("/work/users/s/e/seyoun/CQTL_sQTL/output")

# Loading all the packages 
library(jsonlite)
library(stringr)
library(data.table)
library(dplyr)
library(GenomicRanges)
library(GenomicFeatures)
library(stringr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggvenn)
library(scales)
library(colorspace)
library(ggtext)
library(ggforce)
library(tidyverse)
library(httpgd)
library(plotgardener)
library(colorspace)
library(lme4)
library(lmerTest)

#---------------------------------------------------------------
# Load the data
#---------------------------------------------------------------
## load the sQTL results
response_pbs_results <- readRDS("/work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/conditional_pbs/response_pbs_resInv.rds")
response_fnf_results <- readRDS("/work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/conditional_pbs/response_fnf_resInv.rds")

  # Subset significant reponse sQTL only
  #pbs_sQTL <- response_pbs_results %>%
  #  dplyr::filter(rank  == 0)

  #fnf_sQTL <- response_fnf_results %>%
  #  dplyr::filter(rank  == 0)

  ## load the eQTL results
  eqtl_dir <- "/proj/phanstiel_lab/Data/processed/CQTL/eqtl/"

  #This data set have the signal information
  pbs_eqtl <- fread(paste0(eqtl_dir,"CTL_PEER_k20_genoPC_cond1Mb_topSignals_rsID_signalRanges_distal.csv"))
  fnf_eqtl <- fread(paste0(eqtl_dir,"FNF_PEER_k22_genoPC_cond1Mb_topSignals_rsID_signalRanges_distal.csv"))

  #---------------------------------------------------------------
  # Find the LD with the 0.7 sQTL
  #---------------------------------------------------------------

  annotate_eqtl_ld <- function(sqtl, eqtl, ld_dir, ld_threshold = 0.7) {
    # sqtl: data.table with sQTL info (must include column 'var_id')
    # eqtl: data.table with eQTL info (must include columns 'variantID' & 'gene_symbol')
    # ld_dir: directory where LD files are stored in chromosome-specific subdirectories
    # ld_threshold: R2 threshold for high LD variants
      if(!data.table::is.data.table(sqtl)) {
      sqtl <- as.data.table(sqtl)
      }
    
    # Add new columns for combined annotation and for match count
    sqtl[, eqtl_annotation := NA_character_]
    sqtl[, eqtl_match_count := NA_integer_]
    
    # Iterate over each sQTL
    for(i in seq_len(nrow(sqtl))) {
      var_id <- sqtl$var_id[i]
      # Extract chromosome from var_id (e.g., "chr10" from "chr10:1066481:A:G")
      chr_sub <- strsplit(var_id, ":")[[1]][1]
      ld_file <- file.path(ld_dir, chr_sub, paste0(var_id, ".ld"))
      
      # Check if the LD file exists
      if(file.exists(ld_file)) {
        ld_data <- fread(ld_file)
        # Filter LD data to get variants with R2 > ld_threshold
        ld_hi <- ld_data[R2 > ld_threshold]
        
        if(nrow(ld_hi) > 0) {
          # Get unique variants from LD file (assumed column name "SNP_B")
          ld_variants <- unique(ld_hi$SNP_B)
          # Find intersection between LD variants and eQTL variantID
          common_vars <- intersect(ld_variants, eqtl$variantID)
          
          if(length(common_vars) > 0) {
            # Get matching eQTL rows
            eqtl_matches <- eqtl[variantID %in% common_vars]
            # Merge with LD data to bring in the R2 values
            common_pairs <- merge(eqtl_matches, ld_hi, by.x = "variantID", by.y = "SNP_B")
            
            if(nrow(common_pairs) > 0) {
              # Build annotation strings for each match: "variantID:gene_symbol (R2=0.XXX)"
              pair_strings <- apply(common_pairs, 1, function(r) {
                paste0(r["variantID"], ":", r["gene_id"], ":", r["gene_symbol"], ":", r["signal"], " (R2=", round(as.numeric(r["R2"]), 3), ")")
              })
              annotation_str <- paste(pair_strings, collapse = ";")
              sqtl$eqtl_annotation[i] <- annotation_str
              sqtl$eqtl_match_count[i] <- nrow(common_pairs)
            }
          }
        }
      }
    }
    
    return(sqtl)
  }


  pbs_sqtl_eqtlLD_07 <- annotate_eqtl_ld(response_pbs_results, pbs_eqtl, ld_dir = "/work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld")
  fnf_sqtl_eqtlLD_07 <- annotate_eqtl_ld(response_fnf_results, fnf_eqtl, ld_dir = "/work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld")

  pbs_sqtl_eqtlLD_99 <- annotate_eqtl_ld(response_pbs_results, pbs_eqtl, 
  ld_dir = "/work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld",
  ld_threshold = 0.99)
  fnf_sqtl_eqtlLD_99 <- annotate_eqtl_ld(response_fnf_results, fnf_eqtl, 
  ld_dir = "/work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld",
  ld_threshold = 0.99)

  save(pbs_sqtl_eqtlLD_07, fnf_sqtl_eqtlLD_07, file = "/work/users/s/e/seyoun/CQTL_sQTL/output/revision/data/sqtl_eqtl_LD07.RData")
  save(pbs_sqtl_eqtlLD_99, fnf_sqtl_eqtlLD_99, file = "/work/users/s/e/seyoun/CQTL_sQTL/output/revision/data/sqtl_eqtl_LD99.RData")

  load("/work/users/s/e/seyoun/CQTL_sQTL/output/revision/data/sqtl_eqtl_LD07.RData")

  #---------------------------------------------------------------
  # eQTL and sQTL comparsion whether eQTL and sQTL are independent from gene expression
  #---------------------------------------------------------------
  # Step 1 yes or no of the eqtl and sqtl overlaps.
  # For PBS: add an annotation column ("Yes"/"No")
  pbs_sqtl_eqtlLD_07[, annot := ifelse(is.na(eqtl_match_count), "No", "Yes")]
  # For FNF: add the same annotation column
  fnf_sqtl_eqtlLD_07[, annot := ifelse(is.na(eqtl_match_count), "No", "Yes")]

  # Step 2: Check if the gene symbol is the same.
  #pbs_sqtl_eqtl_ld07_annot <- pbs_sqtl_eqtlLD_07 |> filter(annot == "Yes") 
  #fnf_sqtl_eqtl_ld07_annot <- fnf_sqtl_eqtlLD_07 |> filter(annot == "Yes") 

  pbs_sqtl_eqtlLD_07[, same_gene := sapply(seq_len(.N), function(i) {
    annot <- eqtl_annotation[i]
    gene  <- ensg[i]
    
    # Check if either annotation or gene are missing or empty
    if (is.na(annot) || is.na(gene) || gene == "") return("No")
    
    # Extract all ENSG IDs from the annotation string
    ids <- unlist(str_extract_all(annot, "ENSG\\d+"))
    
    if (length(ids) == 0) {
      return("No")
    } else if (gene %in% ids) {
      return("Yes")
    } else {
      return("No")
    }
  })]

  fnf_sqtl_eqtlLD_07[, same_gene := sapply(seq_len(.N), function(i) {
    annot <- eqtl_annotation[i]
    gene  <- ensg[i]
    
    # Check if either annotation or gene are missing or empty
    if (is.na(annot) || is.na(gene) || gene == "") return("No")
    
    # Extract all ENSG IDs from the annotation string
    ids <- unlist(str_extract_all(annot, "ENSG\\d+"))
    
    if (length(ids) == 0) {
      return("No")
    } else if (gene %in% ids) {
      return("Yes")
    } else {
      return("No")
    }
  })]

  #----------------------------------------------------------------
  # Make a plot allelic heterogeneity
  #----------------------------------------------------------------
  # First plot is for the Upset plot for allelic heterogeneity
  pbs_sQTL_eqtl_summary_by_rank <- pbs_sqtl_eqtlLD_07 %>%
    mutate(
      rank_category = case_when(
        rank == 0 ~ "Rank0",
        rank == 1 ~ "Rank1",
        rank == 2 ~ "Rank2",
        rank == 3 ~ "Rank3"
      ),
      overlap_category = case_when(
        is.na(eqtl_match_count) ~ "No_Overlap",
        eqtl_match_count == 1 & same_gene == "Yes" ~ "Variant_Overlap_Exact_Gene_Match",
        eqtl_match_count == 1 & same_gene == "No" ~ "Variant_Overlap_Gene_Mismatch",
        eqtl_match_count > 1  ~ "Variant_Overlap_Multiple_Gene_Matches"
      )
    ) %>%
    group_by(rank_category, overlap_category) %>%
    summarise(count = n(), .groups = "drop")


  fnf_sQTL_eqtl_summary_by_rank <- fnf_sqtl_eqtlLD_07 %>%
    mutate(
      rank_category = case_when(
        rank == 0 ~ "Rank0",
        rank == 1 ~ "Rank1",
        rank == 2 ~ "Rank2"
      ),
      overlap_category = case_when(
        is.na(eqtl_match_count) ~ "No_Overlap",
        eqtl_match_count == 1 & same_gene == "Yes" ~ "Variant_Overlap_Exact_Gene_Match",
        eqtl_match_count == 1 & same_gene == "No" ~ "Variant_Overlap_Gene_Mismatch",
        eqtl_match_count > 1  ~ "Variant_Overlap_Multiple_Gene_Matches"
      )
    ) %>%
    group_by(rank_category, overlap_category) %>%
    summarise(count = n(), .groups = "drop")



  pbs_sQTL_eqtl_by_rank_two_category <- pbs_sqtl_eqtlLD_07 %>%
    mutate(
      overlap_category = case_when(
      !is.na(eqtl_match_count) & same_gene == "Yes" ~ "Variant_Overlap_with_Gene",
    TRUE ~ "No_Overlap"
      )
    )

  fnf_sQTL_eqtl_by_rank_two_category <- fnf_sqtl_eqtlLD_07 %>%
    mutate(
      overlap_category = case_when(
      !is.na(eqtl_match_count) & same_gene == "Yes" ~ "Variant_Overlap_with_Gene",
    TRUE ~ "No_Overlap"
      )
    )

  #Create the plot for PBS and FNF how the sQTL and eQTL are related

  pbs_sQTL_eqtl_summary_rank0 <- pbs_sQTL_eqtl_summary_by_rank %>%
    filter(rank_category == "Rank0")

  fnf_sQTL_eqtl_summary_rank0 <- fnf_sQTL_eqtl_summary_by_rank %>%
    filter(rank_category == "Rank0")

  # Combine the summaries and add a group column
  combined_summary_rank0 <- bind_rows(
    pbs_sQTL_eqtl_summary_rank0 %>% mutate(group = "PBS"),
    fnf_sQTL_eqtl_summary_rank0 %>% mutate(group = "FN-f")
  ) %>%
    mutate(group = factor(group, levels = c("PBS", "FN-f")))

  #----------------------------------------------------------------
  #Make summary only for the two groups

    pbs_sQTL_eqtl_summary <- pbs_sqtl_eqtlLD_07 %>%
    mutate(
      overlap_category = case_when(
      !is.na(eqtl_match_count) & same_gene == "Yes" ~ "Variant_Overlap_with_Gene",
    TRUE ~ "No_Overlap"
      )
    )%>%
    group_by(overlap_category) %>%
    summarise(count = n(), .groups = "drop")
  
  fnf_sQTL_eqtl_summary <- fnf_sqtl_eqtlLD_07 %>%
    mutate(
      overlap_category = case_when(
      !is.na(eqtl_match_count) & same_gene == "Yes" ~ "Variant_Overlap_with_Gene",
    TRUE ~ "No_Overlap"
      )
    )%>%
    group_by(overlap_category) %>%
    summarise(count = n(), .groups = "drop")


  #Create two categories for the overlap to make a plot
    combined_summary <- bind_rows(
    pbs_sQTL_eqtl_summary %>% mutate(group = "PBS"),
    fnf_sQTL_eqtl_summary %>% mutate(group = "FN-f")
  ) %>%
    mutate(group = factor(group, levels = c("PBS", "FN-f")))


  # Create the grouped barplot
  eqtl_sqtl_overlaps_Barplot <- ggplot(combined_summary, aes(x = overlap_category, y = count, fill = group)) +
    geom_col(position = position_dodge(width = 0.9), width = 0.8) +
    geom_text(aes(label = count),
              position = position_dodge(width = 0.9),
              vjust = -0.5,
              size = 2) +
    scale_fill_manual(values = c("PBS" = "#BFDDFF", "FN-f" = "#FFDDA2")) +
    scale_y_continuous(
      name = "Counts of sQTL-eQTL Overlap",
      trans = scales::pseudo_log_trans(base = 2),
      breaks = c(0, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096),
      labels = scales::comma,
      expand = expansion(mult = c(0, 0.1))
    ) +
    scale_x_discrete(name = "Overlap Category") +
    coord_cartesian(clip = "off") +
    theme(
      strip.placement = "outside",
      axis.line = element_line(linewidth = 0.25),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(color = "black", linewidth = 0.25),
      axis.ticks.length.y = unit(-0.1, "cm"),
      axis.title.x = element_markdown(size = 8, family = "Helvetica", margin = margin(t = 5)),
      axis.title.y = element_markdown(size = 8, family = "Helvetica", margin = margin(r = 5)),
      text = element_text(family = "Helvetica"),
      axis.text.y = element_text(color = "black", size = 6),
      axis.text.x = element_text(color = "black", size = 6, margin = margin(t = 5)),
      strip.background = element_blank(),
      strip.text.x.top = element_text(size = 8, margin = margin(b = 5)),
      panel.background = element_rect(fill = "transparent", color = "transparent"),
      plot.background = element_rect(fill = "transparent", color = "transparent"),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.position.inside = c(0.9, 0.9),
      panel.spacing.y = unit(0.5, "cm")
    )

save(eqtl_sqtl_overlaps_Barplot, 
file = "/work/users/s/e/seyoun/CQTL_sQTL/output/revision/plots/eqtl_sqtl_overlaps_Barplot.rda")
  #----------------------------------------------------------------
  #This is for the exact match of the lead variant
  #----------------------------------------------------------------

  pbs_sqtl_eqtlLD_99[, annot := ifelse(is.na(eqtl_match_count), "No", "Yes")]
  fnf_sqtl_eqtlLD_99[, annot := ifelse(is.na(eqtl_match_count), "No", "Yes")]

  pbs_sqtl_eqtlLD_99[, same_gene := sapply(seq_len(.N), function(i) {
    annot <- eqtl_annotation[i]
    gene  <- ensg[i]
    
    # Check if either annotation or gene are missing or empty
    if (is.na(annot) || is.na(gene) || gene == "") return("No")
    
    # Extract all ENSG IDs from the annotation string
    ids <- unlist(str_extract_all(annot, "ENSG\\d+"))
    
    if (length(ids) == 0) {
      return("No")
    } else if (gene %in% ids) {
      return("Yes")
    } else {
      return("No")
    }
  })]

  fnf_sqtl_eqtlLD_99[, same_gene := sapply(seq_len(.N), function(i) {
    annot <- eqtl_annotation[i]
    gene  <- ensg[i]
    
    # Check if either annotation or gene are missing or empty
    if (is.na(annot) || is.na(gene) || gene == "") return("No")
    
    # Extract all ENSG IDs from the annotation string
    ids <- unlist(str_extract_all(annot, "ENSG\\d+"))
    
    if (length(ids) == 0) {
      return("No")
    } else if (gene %in% ids) {
      return("Yes")
    } else {
      return("No")
    }
  })]


  # Create a summary for the exact match dataset (ld_threshold = 1)
  pbs_exact_summary <- pbs_sqtl_eqtlLD_99 %>%
    mutate(exact_category = case_when(
      eqtl_match_count == 1 & same_gene == "Yes" ~ "Variant_Exact_Match_Gene_Match",
      eqtl_match_count == 1 & same_gene == "No" ~ "Variant_Exact_Match_Gene_Mismatch"
      # If needed, include other conditions (e.g., if eqtl_match_count > 1, you might want to label them separately)
    )) %>%
    filter(!is.na(exact_category)) %>%       # keep only cases with defined category
    group_by(exact_category) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(group = "PBS")

  fnf_exact_summary <- fnf_sqtl_eqtlLD_99 %>%
    mutate(exact_category = case_when(
      eqtl_match_count == 1 & same_gene == "Yes" ~ "Variant_Exact_Match_Gene_Match",
      eqtl_match_count == 1 & same_gene == "No" ~ "Variant_Exact_Match_Gene_Mismatch"
    )) %>%
    filter(!is.na(exact_category)) %>%
    group_by(exact_category) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(group = "FN-f")


  #---------------------------------------------------------------
  # For the multiple independent sQTLs and multiple junctions
  #---------------------------------------------------------------
  # Function to compute minimum LD for a set of variants
  compute_ld_info <- function(variants, ld_dir) {
    # Return a list with the numeric minimum LD and a string with all variant pair details.
    if (length(variants) < 2) return(list(min_ld = NA_real_, ld_details = NA_character_))
    
    ld_values <- c()
    ld_pair_info <- c()
    
    # Loop through each unique pair of variants
    for(i in seq_len(length(variants)-1)) {
      for(j in (i+1):length(variants)) {
        var1 <- variants[i]
        var2 <- variants[j]
        chr_sub <- strsplit(var1, ":")[[1]][1]
        ld_file <- file.path(ld_dir, chr_sub, paste0(var1, ".ld"))
        
        if(file.exists(ld_file)) {
          ld_data <- fread(ld_file)
          # Look for the second variant in the LD file (assume columns "SNP_B" and "R2")
          r_entry <- ld_data[SNP_B == var2, R2]
          if(length(r_entry) > 0) { 
            r_val <- as.numeric(r_entry[1])
            ld_values <- c(ld_values, r_val)
            ld_pair_info <- c(ld_pair_info,
              paste0(var1, "-", var2, " (R2=", round(r_val, 3), ")")
            )
          }
        }
      }
    }
    if (length(ld_values) == 0) {
      return(list(min_ld = NA_real_, ld_details = NA_character_))
    }
    
    return(list(min_ld = min(ld_values, na.rm = TRUE),
                ld_details = paste(ld_pair_info, collapse = "; ")))
  }

  # Compute heterogeneity incorporating LD among the sQTL variants
  genes_with_multiple_sqtl_pbs <- response_pbs_results |> 
    group_by(ensg) |> 
    filter(n() > 1) |> 
    ungroup()


  pbs_allelic_heterogeneity_results <- genes_with_multiple_sqtl_pbs %>%
    group_by(ensg) %>%
    summarize(
      num_variants = n_distinct(var_id),
      num_junctions = n_distinct(phe_id),
      ld_info_out = list(compute_ld_info(unique(var_id), 
                                ld_dir = "/work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld")),
      min_ld = ld_info_out[[1]]$min_ld,
      ld_details = ld_info_out[[1]]$ld_details,
      heterogeneity_type = case_when(
        num_variants == 1 | (num_variants > 1 & !is.na(min_ld) & min_ld >= 0.7) ~ "No heterogeneity",
        num_variants > 1 & num_junctions == 1 ~ "Same junction heterogeneity",
        num_variants > 1 & num_junctions > 1 & (!is.na(min_ld) && min_ld < 0.7) ~ "Multiple junction heterogeneity",
        TRUE ~ "Unclassified"
      ),
      .groups = "drop"
    )


  genes_with_multiple_sqtl_fnf <- response_fnf_results |> 
    group_by(ensg) |> 
    filter(n() > 1) |> 
    ungroup()

  fnf_allelic_heterogeneity_results <- genes_with_multiple_sqtl_fnf %>%
    group_by(ensg) %>%
    summarize(
      num_variants = n_distinct(var_id),
      num_junctions = n_distinct(phe_id),
      ld_info_out = list(compute_ld_info(unique(var_id), 
                                ld_dir = "/work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld")),
      min_ld = ld_info_out[[1]]$min_ld,
      ld_details = ld_info_out[[1]]$ld_details,
      heterogeneity_type = case_when(
        num_variants == 1 | (num_variants > 1 & !is.na(min_ld) & min_ld >= 0.7) ~ "No heterogeneity",
        num_variants > 1 & num_junctions == 1 ~ "Same junction heterogeneity",
        num_variants > 1 & num_junctions > 1 & (!is.na(min_ld) && min_ld < 0.7) ~ "Multiple junction heterogeneity",
        TRUE ~ "Unclassified"
      ),
      .groups = "drop"
    )

    save(pbs_allelic_heterogeneity_results, 
        fnf_allelic_heterogeneity_results, 
        file = "/work/users/s/e/seyoun/CQTL_sQTL/output/revision/data/sqtl_allelic_heterogeneity.RData")

  load("/work/users/s/e/seyoun/CQTL_sQTL/output/revision/data/sqtl_allelic_heterogeneity.RData")
  #CAST (R2=0.023)
  #chr5:96768778:G:A 
  #chr5:96785448:A:C

  #MICA R2=0.085
  #chr6:31437351:T:A
  #chr6:31392564:T:C

  #---------------------------------------------------------------
  # Conditional analysis 
  #---------------------------------------------------------------
  # Load the data for covariates, sQTL and genotype data and splicing 
  # snp_make the subset first from the data

  # Extract all the ld_details strings (remove NA values first)
  ld_details_pbs <- na.omit(pbs_sQTL_eqtl_by_rank_two_category$eqtl_annotation)
  ld_details_fnf <- na.omit(fnf_sQTL_eqtl_by_rank_two_category$eqtl_annotation)

# For PBS
pbs_variant_ids <- unique(unlist(lapply(ld_details_pbs, function(x) {
  # Extract all matches of the pattern: the first four colon-separated fields.
  # This pattern stops at the fourth colon.
  str_extract_all(x, "(chr[0-9XY]+:[0-9]+:[ACGT]+:[ACGT]+)")[[1]]
})))

# For FNF
fnf_variant_ids <- unique(unlist(lapply(ld_details_fnf, function(x) {
  str_extract_all(x, "(chr[0-9XY]+:[0-9]+:[ACGT]+:[ACGT]+)")[[1]]
})))

  # Flatten the list into a single vector and take unique values
  #eqtl_variant_ids <- unique(unlist(c(pbs_variants_list, fnf_variants_list)))

  eqtl_sqtl_variants_overlapsIds <- unique(c(response_pbs_results$var_id, 
  response_fnf_results$var_id, pbs_variant_ids, fnf_variant_ids))


  write.table(eqtl_sqtl_variants_overlapsIds, 
  file = "/work/users/s/e/seyoun/CQTL_sQTL/output/revision/data/e_sqtl_variant_ids.txt", 
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  #---------------------------------------------------------------
  #Load the geno file
  #---------------------------------------------------------------

  geno_dir <- "/work/users/s/e/seyoun/CQTL_sQTL/output/geno/"

  pbs_geno_matrix <- fread(paste0(geno_dir,"pbs_geno/07.subset_sigSNps_eqtl/recodeA_pbs.traw"))
  fnf_geno_matrix <- fread(paste0(geno_dir,"fnf_geno/07.subset_sigSNps_eqtl/recodeA_fnf.traw"))

  # Define your key columns
  key_cols <- c("CHR", "SNP", "(C)M", "POS", "COUNTED", "ALT")

  # Define your key columns
  key_cols <- c("CHR", "SNP", "(C)M", "POS", "COUNTED", "ALT")

  # Merge the PBS and FNF genotype matrices by the key columns.
  combined_geno_matrix <- full_join(pbs_geno_matrix,   fnf_geno_matrix,
                                    by = key_cols)
#----------------------------------------------------------------
# Load the covariate data
#---------------------------------------------------------------
# Load the covariate data

pc5_cov_all <- fread("gtex_cluster/qtltools_prep/covariates_PC5") %>%
  t() %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  `colnames<-`(.[1, ]) %>%
  .[-1, ] %>%
  cbind(sampleID = rownames(.), .) %>%
  mutate(across(-sampleID, as.numeric))

pc4_cov_all <- fread("gtex_cluster/qtltools_prep/covariates_PC4") %>%
  t() %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  `colnames<-`(.[1, ]) %>%
  .[-1, ] %>%
  cbind(sampleID = rownames(.), .) %>%
  mutate(across(-sampleID, as.numeric))

#meta #exclude samples for  ['AM7352', 'AM7244']
load("combined_meta_data.RData")
selected_columns <- c("Donor","ID","Condition","Sex", "Age","FragmentBatch","RIN","RNAextractionKitBatch","RNAshippedDate") #select column needed it based
meta_data <- combined_data %>% dplyr::select(all_of(selected_columns))
#Ancestry 
ancestry_df <- fread("/proj/phanstiel_lab/Data/processed/CQTL/geno/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_COA6_COA7/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_COA6_COA7_predictedAncestry.csv")
ancestry_df$Donor <- sub("^.*_(AM[0-9]+)_.*$", "\\1", ancestry_df$Donor)

ancestry_OA_df <- fread("/proj/phanstiel_lab/Data/processed/CQTL/geno/COA8_OA/ancestry/CQTL_COA8_predictedAncestry.csv")
ancestry_OA_df_filtered <- ancestry_OA_df %>%
  dplyr::filter(grepl("^OA\\d+", Donor))
ancestry_OA_df_filtered$Donor <- gsub("_r2","",ancestry_OA_df_filtered$Donor)

ancestry_cqtl <-rbind(ancestry_df,ancestry_OA_df_filtered)
meta_cqtl <- merge(meta_data,ancestry_cqtl,by="Donor",all.x=TRUE)
#write.table(meta_cqtl, file = "output/clu_fnf/meta_cqtl",sep='\t',quote=F,row.names=F,col.names=T)
meta_cqtl <- fread("clu_fnf/meta_cqtl")
meta_catl_simple <- meta_cqtl %>%
  dplyr::filter(Condition %in% c("CTL","FNF")) %>%
  dplyr::select(-Sex,-Age,-FragmentBatch,-RIN,-RNAextractionKitBatch,-RNAshippedDate,-Predicted_Ancestry) %>%
  dplyr::filter(!str_detect(Donor, 'AM7352|AM7244'))

# Normalized PSI value
ctl_fnf_ratio <- fread("clu_fnf/ratio_fnf.txt") |> as.data.frame()
rownames(ctl_fnf_ratio) <- ctl_fnf_ratio$Junction
ratios_fnf <- ctl_fnf_ratio[,-1]

pbs_QTL_cond_final <- readRDS("01.qtltools_re/pbs_cond_sig_beta_maf_added.rds")
fnf_QTL_cond_final <- readRDS("01.qtltools_re/fnf_cond_sig_beta_maf_added.rds")

intron_test <- c(pbs_QTL_cond_final$phe_id, fnf_QTL_cond_final$phe_id) |> unique()
psi_matrix <- ratios_fnf[rownames(ratios_fnf) %in% intron_test,] |> as.data.frame()
intronID <- rownames(psi_matrix)
psi_ratio <- apply(psi_matrix, 2, as.numeric)
rownames(psi_ratio) <- intronID
psi_ratio.df <- t(psi_ratio) |> as.data.frame()


#----------------------------------------------------------------
#Function for the Conditional analysis
#----------------------------------------------------------------

analyze_conditional_qtl <- function(qtl_results, norm_counts, geno_qtl, meta_df, covariates_df, cond) {
  # Initialize output columns
  qtl_results$eQTL_test <- NA
  qtl_results$eQTL_r2 <- NA
  qtl_results$anova_conditional_pval <- NA
  qtl_results$f_stat <- NA
  qtl_results$cond_summary_pvalue <- NA
  qtl_results$cond_summary_beta <- NA
  qtl_results$cond_summary_se <- NA

  for (i in 1:nrow(qtl_results)) {
    cat("Processing", i, "\n")
    FeatureID <- qtl_results$phe_id[i]
    sqtl_variant <- qtl_results$var_id[i]

    if (!is.na(qtl_results$eqtl_annotation[i])) {
      annos <- unlist(strsplit(qtl_results$eqtl_annotation[i], ";"))
      matches <- annos[grepl(qtl_results$ensg[i], annos)]
      chosen_annotation <- if (length(matches) > 0) matches[1] else annos[1]
      anno_parts <- unlist(strsplit(chosen_annotation, ":"))
      if (length(anno_parts) >= 7) {
        eqtl_variant <- paste(anno_parts[1:4], collapse = ":")
        eqtl_r2 <- anno_parts[7]
      } else {
        cat("Warning: Invalid annotation format for iteration", i, "\n")
        next
      }
    } else {
      cat("Warning: Missing annotation for iteration", i, "\n")
      next
    }

    qtl_results$eQTL_test[i] <- eqtl_variant
    qtl_results$eQTL_r2[i] <- eqtl_r2

    # Extract normalized expression
    norm_counts_data <- norm_counts[FeatureID, , drop = FALSE]
    norm_counts_df <- data.frame(
      sampleID = colnames(norm_counts_data),
      psi = as.double(norm_counts_data[1, ])
    )

    # Get sQTL genotype
    geno_sqtl_data <- geno_qtl %>% filter(SNP == sqtl_variant)
    if (nrow(geno_sqtl_data) == 0) {
      cat("Warning: sQTL variant", sqtl_variant, "not found for iteration", i, "\n")
      next
    }
    geno_sqtl_data <- geno_sqtl_data %>% column_to_rownames("SNP")
    geno_sqtl_df <- data.frame(
      sampleID = colnames(geno_sqtl_data),
      sqtl_genotype = as.numeric(geno_sqtl_data[1, ])
    )

    # Get eQTL genotype
    geno_eqtl_data <- geno_qtl %>% filter(SNP == eqtl_variant)
    if (nrow(geno_eqtl_data) == 0) {
      cat("Warning: eQTL variant", eqtl_variant, "not found for iteration", i, "\n")
      next
    }
    geno_eqtl_data <- geno_eqtl_data %>% column_to_rownames("SNP")
    geno_eqtl_df <- data.frame(
      sampleID = colnames(geno_eqtl_data),
      eqtl_genotype = as.numeric(geno_eqtl_data[1, ])
    )

    # Merge data
    meta_work <- meta_df %>% dplyr::rename(sampleID = ID) %>% dplyr::select(sampleID, Donor, Condition)
    df_merged <- geno_sqtl_df %>%
      left_join(geno_eqtl_df, by = "sampleID") %>%
      left_join(norm_counts_df, by = "sampleID") %>%
      left_join(meta_work, by = "sampleID") %>%
      left_join(covariates_df, by = "sampleID") %>%
      filter(Condition == cond)

    if (nrow(df_merged) == 0) {
      cat("Warning: No samples for condition", cond, "in iteration", i, "\n")
      next
    }

    df_merged <- df_merged %>%
      mutate(
        Donor = as.factor(Donor),
        sqtl_genotype = as.numeric(sqtl_genotype),
        eqtl_genotype = as.numeric(eqtl_genotype),
        sampleID = as.character(sampleID)
      )

    # Build formulas
    covar_columns <- colnames(df_merged)[!(colnames(df_merged) %in%
      c("sampleID", "Donor", "Condition", "psi", "sqtl_genotype", "eqtl_genotype"))]
    covar_formula <- if (length(covar_columns) > 0) {
      paste(covar_columns, collapse = " + ")
    } else {
      ""
    }

    reduced_formula <- as.formula(if (covar_formula == "") {
      "psi ~ sqtl_genotype"
    } else {
      paste("psi ~ sqtl_genotype +", covar_formula)
    })

    conditional_formula <- as.formula(if (covar_formula == "") {
      "psi ~ sqtl_genotype + eqtl_genotype"
    } else {
      paste("psi ~ sqtl_genotype + eqtl_genotype +", covar_formula)
    })

    # Fit models
    reduced_model <- lm(reduced_formula, data = df_merged)
    full_model <- lm(conditional_formula, data = df_merged)
    model_comp <- anova(reduced_model, full_model)

    qtl_results$anova_conditional_pval[i] <- model_comp$`Pr(>F)`[2]
    qtl_results$f_stat[i] <- model_comp$F[2]

    # Extract eQTL coefficient info
    coef_table <- summary(full_model)$coefficients
    if ("eqtl_genotype" %in% rownames(coef_table)) {
      qtl_results$cond_summary_pvalue[i] <- coef_table["eqtl_genotype", "Pr(>|t|)"]
      qtl_results$cond_summary_beta[i] <- coef_table["eqtl_genotype", "Estimate"]
      qtl_results$cond_summary_se[i] <- coef_table["eqtl_genotype", "Std. Error"]
    } else {
      cat("Warning: eqtl_genotype not found in coefficients for iteration", i, "\n")
    }
  }

  return(qtl_results)
}




pbs_eqtl_sqtl_overlaps <- pbs_sQTL_eqtl_by_rank_two_category |> 
filter(overlap_category == "Variant_Overlap_with_Gene")

fnf_eqtl_sqtl_overlaps <- fnf_sQTL_eqtl_by_rank_two_category |> 
filter(overlap_category == "Variant_Overlap_with_Gene")

Conditional_pbs_re <- analyze_conditional_qtl(qtl_results=pbs_eqtl_sqtl_overlaps,
                                             norm_counts=psi_ratio,
                                             geno_qtl=combined_geno_matrix,
                                             meta_df=meta_catl_simple,
                                             covariates_df=pc5_cov_all,
                                             cond="CTL")


Conditional_fnf_re <- analyze_conditional_qtl(qtl_results=fnf_eqtl_sqtl_overlaps,
                                             norm_counts=psi_ratio,
                                             geno_qtl=combined_geno_matrix,
                                             meta_df=meta_catl_simple,
                                             covariates_df=pc4_cov_all,
                                             cond="FNF")



save(Conditional_pbs_re, Conditional_fnf_re,
     file = "/work/users/s/e/seyoun/CQTL_sQTL/output/revision/data/conditional_analysis.RData")


#----------------------------------------------------------------
#FDR Correction for the conditional analysis
#----------------------------------------------------------------

Conditional_pbs_re$anova_cond_fdr <- p.adjust(Conditional_pbs_re$anova_conditional_pval, method = "fdr")
Conditional_fnf_re$anova_cond_fdr <- p.adjust(Conditional_fnf_re$anova_conditional_pval, method = "fdr")

Conditional_pbs_re |> filter(anova_cond_fdr < 0.05) 
Conditional_fnf_re |> filter(anova_cond_fdr < 0.05) 


#----------------------------------------------------------------