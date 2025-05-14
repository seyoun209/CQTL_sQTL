## Prep for the bedtools and covariates for the QTLtools - Modified for PBS vs OA
## Author: Seyoun Byun (Modified version)
#-------------------------------------------------------------------------------
setwd("/work/users/s/e/seyoun/CQTL_sQTL")
#library(leafcutter)
library(data.table)
library(dplyr)
library(yaml)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
print(args)
qqnorm_path <- args[1]
bed_dir <- args[2]
pc_dir <- args[3]
pbs_oa_prefix <- ifelse(length(args) >= 4, args[4], "pbs_oa")

if (!dir.exists(bed_dir)) {
  dir.create(bed_dir, recursive = TRUE)
}

# Create a clean tabix script
cat("#!/bin/bash\n", file = paste0(bed_dir,"/","tabix.sh"))

# Pattern for PBS vs OA
pattern <- paste0(pbs_oa_prefix, "_perind.counts.filtered.gz.qqnorm*")
file_list <- list.files(qqnorm_path, pattern = pattern, full.names = TRUE)
filtered_files <- grep("\\.tbi$", file_list, value = TRUE, invert = TRUE)
print(filtered_files)

chromosome_names <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                      "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")

# First loop: check files, similar to original script
for (chr_name in chromosome_names) {
  filtered_fi <- grep(paste0(chr_name, "\\b"), filtered_files, value = TRUE)
  filtered_fi <- gsub("//", "/", filtered_fi)
  if (length(filtered_fi) > 0) {
    bed_file <- fread(filtered_fi[1]) |> as.data.frame()  # Ensure at least one file is there
    # Just validation, no processing yet
  } else {
    cat("No files found for", chr_name, "\n")
  }
}

# Create a list to track which chromosomes we've already processed
processed_chrs <- c()

# Second loop: actually process each chromosome
for (chr_name in chromosome_names) {
  # Skip if we've already processed this chromosome
  if (chr_name %in% processed_chrs) {
    next
  }
  
  filtered_fi <- grep(paste0(chr_name,"\\b"), filtered_files, value = TRUE)
  if (length(filtered_fi) == 0) {
    cat("No files found for", chr_name, "\n")
    next
  }
  
  filtered_fi <- gsub("//", "/", filtered_fi[1])  # Take only the first matching file
  bed_file <- fread(filtered_fi) |> as.data.frame()
  
  bed_file$strand <- do.call(rbind,strsplit(do.call(rbind,strsplit(bed_file$ID, ":"))[,4],"_"))[,3]
  bed_file$clusterID <- do.call(rbind,strsplit(bed_file$ID, ":"))[,4]
  bed_file_reordered <- bed_file[, c("#Chr", "start", "end", "ID", "clusterID", "strand", setdiff(colnames(bed_file), c("#Chr", "start", "end", "ID", "clusterID", "strand")))]
  
  chr_bed <- bed_file_reordered %>%
    filter(`#Chr` == chr_name)
  
  # Define BED file name for PBS vs OA
  bed_file_name <- paste0(pbs_oa_prefix, "_qqnorm_", chr_name, ".bed")
  bed_file_path <- file.path(bed_dir, bed_file_name)
  write.table(chr_bed, bed_file_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote=F)
  
  gzip_command <- paste("bgzip", bed_file_path)
  cat(gzip_command, "\n", file = paste0(bed_dir,"/","tabix.sh"), append = TRUE)
  tabix_command <- paste0("tabix -p bed ", bed_file_path, ".gz")
  cat(tabix_command, "\n", file = paste0(bed_dir,"/","tabix.sh"), append = TRUE)
  
  # Mark this chromosome as processed
  processed_chrs <- c(processed_chrs, chr_name)
}
#-------------------------------------------------------------------------------
# adding covariates
#-------------------------------------------------------------------------------

# samplesheet
aligned_samplesheet_path <- "aligned_samplesheet.txt"
donor_samples_path <- "donor_samples.txt"
rna_extraction_path <- "rna_extraction.txt"
conditions_to_include <- c("CTL", "OA")  # Modified for PBS vs OA

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

# Process ID column
combined_data <- combined_data %>%
  mutate(ID = paste(Donor, Condition, Tech_Rep, Sex, sep = "_"))

# Omit samples and filter rows based on configuration and conditions
combined_data <- combined_data[!combined_data$Donor %in% config$samples_to_omit, ] %>%
  filter(Condition %in% conditions_to_include)

# Select columns
selected_columns <- c("ID")
final_meta_data <- combined_data %>% dplyr::select(all_of(selected_columns))

# PCA for PBS vs OA
pca_raw <- fread(pc_dir) |>
  t() |>
  as.data.frame()
pca_raw <- pca_raw[-1,]
colnames(pca_raw) <- c(sprintf("PC%s",seq(1:20)))
pca_raw_nm <- cbind(rownames(pca_raw),pca_raw)
colnames(pca_raw_nm)[1] <-c("ID")
pc20_pbs_oa <- pca_raw_nm[,1:21]

# Merge with metadata
splicingPCA_df <- merge(final_meta_data, pc20_pbs_oa, by="ID", all=FALSE)

# Change to factor
for (col in colnames(splicingPCA_df)) {
  # Convert only columns not starting with "PC" to factor
  if (!grepl("^PC", col)) {
    splicingPCA_df[[col]] <-  as.factor(splicingPCA_df[[col]])
  }
}

# Rename columns
#colnames(splicingPCA_df)[which(colnames(splicingPCA_df) == 'RNAshippedDate')] <- c("SeqeuncingBatch")
splicingPCA_numeric <- sapply(splicingPCA_df[,2:21],as.numeric) |> as.data.frame() 
splicingPCA_num <- cbind(splicingPCA_df[,1],splicingPCA_numeric)

# Read in genotype PCA data
pca_pbs <- fread("./output/geno/pbs_geno/05.pca/CQTL.eigenvec")
pca_oa <- fread("./output/geno/oa_geno/05.pca/CQTL.eigenvec")  # Updated to OA

pca_pbs <- pca_pbs[,-1]
# set names
names(pca_pbs)[1] <- "ID"
names(pca_pbs)[2:ncol(pca_pbs)] <- paste0("GenoPC", 1:(ncol(pca_pbs)-1))

pca_oa <- pca_oa[,-1]
# set names
names(pca_oa)[1] <- "ID"
names(pca_oa)[2:ncol(pca_oa)] <- paste0("GenoPC", 1:(ncol(pca_oa)-1))

# Combine PBS and OA genotype PCs
geno_pcs <- rbind(pca_pbs[,1:17], pca_oa)
geno_pcs_df <- geno_pcs[,1:5] |> as.data.frame() # Change the number of the Geno PCs to 4

# Create covariate files for each PC
for (i in 1:20) {
  print(i)
  col_n <- i
  cov_merged <- merge(splicingPCA_num[,1:(i+1)], geno_pcs_df, by='ID', all.x=TRUE)
  cov_transposed <- t(cov_merged[, -1])
  colnames(cov_transposed) <- cov_merged$ID
  cov_fi <- cbind(rownames(cov_transposed), cov_transposed)
  colnames(cov_fi)[1] <- 'id' 
  write.table(cov_fi, file=paste0(bed_dir,"/","covariates_PC",col_n), sep='\t', quote=F, row.names=F, col.names=T)
}