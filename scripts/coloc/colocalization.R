# Finding the Collocalization
## Author: Seyoun Byun
## Date: 06.27.2024
## Edited:
#-------------------------------------------------------------------------------
setwd("/work/users/s/e/seyoun/CQTL_sQTL/")
library(tidyverse)
library(data.table)
library(coloc)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
source("scripts/utils/utils.R")

#-------------------------------------------------------------------------------
#Input data - sQTL
#-------------------------------------------------------------------------------
#QTL
#sig_response_pbs_nolmer <- fread("output/01.qtltools_re/response_qtl/reseponsesQTL_PBS_significant.csv")
#sig_response_fnf_nolmer <- fread("output/01.qtltools_re/response_qtl/reseponsesQTL_FNF_significant.csv")
pbs_sig_qtl_cond_annot <- readRDS("output/01.qtltools_re/pbs_cond_sig_beta_maf_added.rds")
fnf_sig_qtl_cond_annot <- readRDS("output/01.qtltools_re/fnf_cond_sig_beta_maf_added.rds")
pbs_sig_qtl_cond_annot[, lead_sqtl := process_var_id(pbs_sig_qtl_cond_annot$var_id)]
fnf_sig_qtl_cond_annot[, lead_sqtl := process_var_id(fnf_sig_qtl_cond_annot$var_id)]


#pbs_sig_qtl_nodup <- pbs_sig_qtl_cond_annot %>% 
#  group_by(lead_sqtl) %>%
#  slice_min(order_by= bwd_pval, n=1)%>%
#  ungroup()


#fnf_sig_qtl_nodup <- fnf_sig_qtl_cond_annot %>% 
#  group_by(lead_sqtl) %>%
#  slice_min(order_by= bwd_pval, n=1)%>%
#  ungroup()

sqtl_pbs_list <- list()
sqtl_fnf_list <- list()

for (chr in 1:22) {

qtl_pbs_chr <- pbs_sig_qtl_nodup %>%
  dplyr::filter(var_chr == paste0("chr", chr))

qtl_fnf_chr <- fnf_sig_qtl_nodup %>%
  dplyr::filter(var_chr == paste0("chr", chr))

# Store the QTL data in the list
sqtl_pbs_list[[paste0("chr", chr)]] <- qtl_pbs_chr
sqtl_fnf_list[[paste0("chr", chr)]] <- qtl_fnf_chr
# Optional: Print a message to track progress
cat("Processed chromosome", chr, "\n")

}



#-------------------------------------------------------------------------------
#Input data - GWAS
#-------------------------------------------------------------------------------
# Read in corresponding subtype leads and their LD buddies

OAsubtypes <- c("AllOA", "FingerOA", "HandOA", "HipOA", "KneeHipOA", "KneeOA",
                "SpineOA", "THR", "ThumbOA", "TJR", "TKR")

subtype_leads <- fread(paste0("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg38_LD/",
                                 OAsubtypes[1],
                                 "/leads/ALL_",OAsubtypes[1],"_leads_ld_final.csv"))

subtype_leads_oaall_eur <- fread(paste0("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg38_LD/",
                                 OAsubtypes[1],
                                 "/leads/EUR_",OAsubtypes[1],"_leads_ld_final.csv")) |> 
  dplyr::select("CHR:hg38POS","ldbuddy_CHR:hg38POS",  "ldbuddy_R2", "p")

# Get signal region of + or - for the GWAS
#min_region <- gwas_variant$hg38pos - 250000
#max_region <- gwas_variant$hg38pos + 250000

gwas_summary_stats_list <- list()

for (chr in 1:22) {
  # Construct the file path
  file_path <- paste0("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg38_LD/",
                      OAsubtypes[1],
                      "/summary_stats/",
                      OAsubtypes[1],"_chr", chr, ".csv")
  
  # Read the file
  gwas_summary_stat <- fread(file_path, data.table = FALSE)
  # Remove rows where hg38POS is NA
  gwas_summary_stat <- gwas_summary_stat[!is.na(gwas_summary_stat$hg38pos), ]
  
  
  # Store the result in the list, using the chromosome number as the name
  gwas_summary_stats_list[[paste0("chr", chr)]] <- gwas_summary_stat
  
  # Optional: Print a message to track progress
  cat("Processed chromosome", chr, "\n")
}

save(gwas_summary_stats_list,file="output/coloc/gwas_summary_stat_list.rda")

#Remove duplicate by max sample size
library(dplyr)

for (chr in 1:22) {
  chr_name <- paste0("chr", chr)
  
  if (!is.null(gwas_summary_stats_list[[chr_name]])) {
    # Identify duplicate rows based on `CHR:hg38POS`
    duplicates <- gwas_summary_stats_list[[chr_name]] %>%
      filter(duplicated(`CHR:hg38POS`) | duplicated(`CHR:hg38POS`, fromLast = TRUE))
    
    # Handle duplicates by keeping the row with max SampleSize
    deduplicated <- duplicates %>%
      group_by(`CHR:hg38POS`) %>%
      slice_max(order_by = SampleSize, n = 1, with_ties = FALSE) %>%
      ungroup()
    
    # Non-duplicate rows
    non_duplicates <- gwas_summary_stats_list[[chr_name]] %>%
      filter(!`CHR:hg38POS` %in% duplicates$`CHR:hg38POS`)
    
    # Combine deduplicated rows with non-duplicate rows
    gwas_summary_stats_list[[chr_name]] <- bind_rows(non_duplicates, deduplicated)
    
    # Print progress
    cat("Processed", chr_name, "\n")
    cat("Rows after handling duplicates:", nrow(gwas_summary_stats_list[[chr_name]]), "\n")
    cat("Unique positions:", length(unique(gwas_summary_stats_list[[chr_name]]$`CHR:hg38POS`)), "\n\n")
  } else {
    cat("No data found for", chr_name, "\n\n")
  }
}


for (chr in 1:22) {
  chr_name <- paste0("chr", chr)
  
  if (!is.null(gwas_summary_stats_list[[chr_name]])) {
    # Convert EAF to MAF and filter out rows where MAF is 0
    gwas_summary_stats_list[[chr_name]] <- gwas_summary_stats_list[[chr_name]] %>%
      mutate(MAF = ifelse(EAF < 0.5, EAF, 1 - EAF)) %>%
      filter(MAF != 0)
    
    # Print progress
    cat("Processed", chr_name, "\n")
    cat("Rows after MAF conversion and filtering:", nrow(gwas_summary_stats_list[[chr_name]]), "\n\n")
  } else {
    cat("No data found for", chr_name, "\n\n")
  }
}




case_control_sizes <- read_csv("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg38_LD/Case_Control_sampleSizes.csv")
#Fraction sizes
#fraction_cases <- case_control_sizes$Max_Cases[1]/(case_control_sizes$Max_Cases[1] + case_control_sizes$Max_Controls[1])
# Usage
fraction_cases <- calculate_case_fraction(case_control_sizes)

#-------------------------------------------------------------------------------
#Run coloc
#-------------------------------------------------------------------------------
n <- 101
coloc_result_chr1 <- coloc.abf(dataset1 = list(pvalues = sqtl_fnf_list$chr1$bwd_pval,
                                          N = n,
                                          MAF = sqtl_fnf_list$chr1$MAF,
                                          type = "quant",
                                          beta = sqtl_fnf_list$chr1$bwd_slope,
                                          snp = sqtl_fnf_list$chr1$lead_sqtl),
                          dataset2 = list(beta = gwas_summary_stats_list$chr1$BETA,
                                          s = fraction_cases$fraction_cases,
                                          N = fraction_cases$total_samples,
                                          type = "cc",
                                          snp = gwas_summary_stats_list$chr1$`CHR:hg38POS`,
                                          MAF = gwas_summary_stats_list$chr1$MAF,
                                          pvalues = gwas_summary_stats_list$chr1$p))



sensitivity(coloc_result_chr1,rule="H4 > 0.2") 


# Fnf

coloc_results <- list()

for (chr in 1:22) {
  chr_name <- paste0("chr", chr)
  
  if (!is.null(sqtl_fnf_list[[chr_name]]) && !is.null(gwas_summary_stats_list[[chr_name]])) {
    coloc_result <- coloc.abf(
      dataset1 = list(
        pvalues = sqtl_fnf_list[[chr_name]]$bwd_pval,
        N = n,
        MAF = sqtl_fnf_list[[chr_name]]$MAF,
        type = "quant",
        beta = sqtl_fnf_list[[chr_name]]$bwd_slope,
        snp = sqtl_fnf_list[[chr_name]]$lead_sqtl
      ),
      dataset2 = list(
        beta = gwas_summary_stats_list[[chr_name]]$BETA,
        s = fraction_cases$fraction_cases,
        N = fraction_cases$total_samples,
        type = "cc",
        snp = gwas_summary_stats_list[[chr_name]]$`CHR:hg38POS`,
        MAF = gwas_summary_stats_list[[chr_name]]$MAF,
        pvalues = gwas_summary_stats_list[[chr_name]]$p
      )
    )
    
    # Store the result in the list
    coloc_results[[chr_name]] <- coloc_result
    
    # Print progress
    cat("Processed", chr_name, "\n")
  } else {
    cat("No data found for", chr_name, "\n")
  }
}


save(coloc_results, file = "output/coloc/coloc_results_fnf.rda")



# PBS

coloc_results <- list()

for (chr in 1:22) {
  chr_name <- paste0("chr", chr)
  
  if (!is.null(sqtl_fnf_list[[chr_name]]) && !is.null(gwas_summary_stats_list[[chr_name]])) {
    coloc_result <- coloc.abf(
      dataset1 = list(
        pvalues = sqtl_fnf_list[[chr_name]]$bwd_pval,
        N = n,
        MAF = sqtl_fnf_list[[chr_name]]$MAF,
        type = "quant",
        beta = sqtl_fnf_list[[chr_name]]$bwd_slope,
        snp = sqtl_fnf_list[[chr_name]]$lead_sqtl
      ),
      dataset2 = list(
        beta = gwas_summary_stats_list[[chr_name]]$BETA,
        s = fraction_cases$fraction_cases,
        N = fraction_cases$total_samples,
        type = "cc",
        snp = gwas_summary_stats_list[[chr_name]]$`CHR:hg38POS`,
        MAF = gwas_summary_stats_list[[chr_name]]$MAF,
        pvalues = gwas_summary_stats_list[[chr_name]]$p
      )
    )
    
    # Store the result in the list
    coloc_results[[chr_name]] <- coloc_result
    
    # Print progress
    cat("Processed", chr_name, "\n")
  } else {
    cat("No data found for", chr_name, "\n")
  }
}





#senstitivity
print_sensitivity <- function(coloc_results) {
  for (chr in 1:22) {
    chr_name <- paste0("chr", chr)
    
    if (!is.null(coloc_results[[chr_name]])) {
      result <- coloc_results[[chr_name]] 
      sensitivity_H3 <- sensitivity(result, rule = "H3 > 0.5")
      sensitivity_H4 <- sensitivity(result, rule = "H4 > 0.5")
      
      cat("Chromosome:", chr_name, "\n")
      cat("  Sensitivity H3 > 0.5:", sensitivity_H3, "\n")
      cat("  Sensitivity H4 > 0.5:", sensitivity_H4, "\n\n")
    } else {
      cat("No data found for", chr_name, "\n\n")
    }
  }
}

print_sensitivity(coloc_results)
