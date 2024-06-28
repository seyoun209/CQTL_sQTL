#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Exactly three arguments must be supplied: base_ld_dir, output_dir, and ld_threshold", call.=FALSE)
}

base_ld_dir <- args[1]
output_dir <- args[2]
ld_threshold <- args[3]

setwd("/work/users/s/e/seyoun/CQTL_sQTL/")
library(tidyverse)
library(data.table)
library(coloc)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(dplyr)
source("scripts/utils/utils.R")

# Function to extract data from .ld files
extract_ld_data <- function(ld_file) {
  ld_data <- read.table(ld_file, header = TRUE)
  lead_qtl <- sub(".ld", "", basename(ld_file))
  ld_data <- ld_data %>%
    dplyr::select(c(SNP_A, SNP_B, R2)) %>%
    dplyr::rename("leadsnp" = SNP_A, "ldbuddy" = SNP_B, "ldbuddy_R2" = R2)
  return(ld_data)
}

# Function to process LD data for a given directory
process_ld_data <- function(base_ld_dir, output_dir, ld_threshold) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Process and save combined LD data for each chromosome
  for (chr in 1:22) {
    chr_ld_dir <- file.path(base_ld_dir, paste0("chr", chr))
    ld_files <- list.files(chr_ld_dir, pattern = "\\.ld$", full.names = TRUE)
    chr_ld_data <- lapply(ld_files, extract_ld_data)
    combined_ld_data <- bind_rows(chr_ld_data)
    
    output_file <- file.path(output_dir, paste0("combined_ld_data_chr", chr, ".csv"))
    write.csv(combined_ld_data, output_file, quote=FALSE, col.names=TRUE, row.names=FALSE)
    
    cat("Processed and saved:", output_file, "\n")
  }
  
  # Process and combine LD data for all chromosomes
  combined_ld_data <- list()
  for (chr in 1:22) {
    combined_ld_qtl <- fread(file.path(output_dir, paste0("combined_ld_data_chr", chr, ".csv")))
    combined_ld_qtl$leadsnp <- sapply(combined_ld_qtl$leadsnp, process_var_id)
    combined_ld_qtl$ldbuddy <- sapply(combined_ld_qtl$ldbuddy, process_var_id)
    combined_ld_data[[paste0("chr", chr)]] <- combined_ld_qtl
    
    cat("Processed:", "chr", chr, "\n")
  }
  
  save(combined_ld_data, file=file.path(output_dir, paste0("ld", ld_threshold, "_combined_leadsnps.rda")))
  cat("Saved combined LD data to:", file.path(output_dir, paste0("ld", ld_threshold, "_combined_leadsnps.rda")), "\n")
}



  