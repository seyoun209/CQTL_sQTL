## Prep for the bedtools and covariates for the QTLtools
## Author: Seyoun Byun
#-------------------------------------------------------------------------------
setwd("/work/users/s/e/seyoun/CQTL_sQTL")
library(leafcutter)
library(data.table)
library(dplyr)
library(yaml)
library(readr)

# args <- c("/work/users/s/e/seyoun/CQTL_sQTL/output/gtex_cluster_wasp/ctl_fnf.leafcutter.bed.gz" ,
#           "/work/users/s/e/seyoun/CQTL_sQTL/output/gtex_cluster_wasp/qtltools_prep",
#  "/work/users/s/e/seyoun/CQTL_sQTL/output/gtex_cluster_wasp/ctl_fnf.leafcutter.PCs.txt",
#  "wasp"
#  )
# 
# args <- c("/work/users/s/e/seyoun/CQTL_sQTL/output/gtex_cluster/",
#           "ctl_fnf_perind.counts.filtered.gz.qqnorm*",
#           "/work/users/s/e/seyoun/CQTL_sQTL/output/gtex_cluster/qtltools_prep",
#           "/work/users/s/e/seyoun/CQTL_sQTL/output/gtex_cluster/ctl_fnf.leafcutter.PCs.txt",
#           "wasp"
# )


args <- commandArgs(trailingOnly = TRUE)
print(args)
qqnorm_path <- args[1]
#bed_file <- fread(args[1]) |> as.data.frame()
pattern <- "ctl_fnf_perind.counts.filtered.gz.qqnorm*"
bed_dir <- args[2]
pc_dir <- args[3]
#pc_dir_wasp <- args[4]
use_wasp_id  <- tail(args, 1) == "wasp"


if (!dir.exists(bed_dir)) {
  dir.create(bed_dir, recursive = TRUE)
}

file_list <- list.files(qqnorm_path, pattern = pattern, full.names = TRUE)
filtered_files <- grep("\\.tbi$", file_list, value = TRUE, invert = TRUE)
print(filtered_files)

chromosome_names <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                      "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")

for (chr_name in chromosome_names) {
  filtered_fi <- grep(paste0(chr_name, "\\b"), filtered_files, value = TRUE)
  filtered_fi <- gsub("//", "/", filtered_fi)
  if (length(filtered_fi) > 0) {
    bed_file <- fread(filtered_fi[1]) |> as.data.frame()  # Ensure at least one file is there
    # Remaining processing...
  } else {
    cat("No files found for", chr_name, "\n")
  }
}

for (chr_name in chromosome_names) {
  filtered_fi <- grep(paste0(chr_name,"\\b"), filtered_files, value = TRUE)
  filtered_fi <- gsub("//", "/", filtered_fi)
  bed_file <- fread(filtered_fi) |> as.data.frame()
  bed_file$strand <- do.call(rbind,strsplit(do.call(rbind,strsplit(bed_file$ID, ":"))[,4],"_"))[,3]
  bed_file$clusterID <- do.call(rbind,strsplit(bed_file$ID, ":"))[,4]
  bed_file_reordered <- bed_file[, c("#Chr", "start", "end", "ID", "clusterID", "strand", setdiff(colnames(bed_file), c("#Chr", "start", "end", "ID", "clusterID", "strand")))]
  
  chr_bed <- bed_file_reordered %>%
    filter(`#Chr` == chr_name)
  
  # Define BED file name
  bed_file_name <- paste0("ctrvsfnf_qqnorm_", chr_name, ".bed")
  bed_file_path <- file.path(bed_dir, bed_file_name)
  write.table(chr_bed, bed_file_path, sep = "\t", row.names = FALSE, col.names = TRUE,quote=F)
  
  gzip_command <- paste("bgzip", bed_file_path)
  cat(gzip_command, "\n", file = paste0(bed_dir,"/","tabix.sh"), append = TRUE)
  tabix_command <- paste0("tabix -p bed ", bed_file_path, ".gz")
  cat(tabix_command, "\n", file = paste0(bed_dir,"/","tabix.sh"), append = TRUE)
}



#-------------------------------------------------------------------------------
# adding covariates
#-------------------------------------------------------------------------------

# samplesheet
aligned_samplesheet_path <- "aligned_samplesheet.txt"
donor_samples_path <- "donor_samples.txt"
rna_extraction_path <- "rna_extraction.txt"
conditions_to_include <- unlist(strsplit("CTL FNF", " ")) # Split the first argument into condition values
#use_wasp_id <- tail(args, 1) == "wasp" # Check if the last argument is "wasp"

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
# Omit samples and filter rows based on configuration and conditions
combined_data <- combined_data[!combined_data$Donor %in% config$samples_to_omit, ] %>%
  filter(Condition %in% conditions_to_include)

# Select and rename columns dynamically, based on the use of wasp ID
selected_columns <- c("ID","FragmentBatch","RNAextractionKitBatch","RNAshippedDate")
selected_columns_wasp <- c("ID_wasp","FragmentBatch","RNAextractionKitBatch","RNAshippedDate")

final_meta_data <- combined_data %>% dplyr::select(all_of(selected_columns))
final_meta_data_wasp <- combined_data %>% dplyr::select(all_of(selected_columns_wasp))


#pca for ctl vs fnf
pca_raw <- fread(pc_dir) |>
  t() |>
  as.data.frame()
pca_raw <- pca_raw[-1,]
colnames(pca_raw) <- c(sprintf("PC%s",seq(1:20)))
pca_raw_nm <- cbind(rownames(pca_raw),pca_raw)
colnames(pca_raw_nm)[1] <-c("ID")
pc20_ctl_fnf <- pca_raw_nm[,1:21]

#Merge with the
if(use_wasp_id) {
  splicingPCA_df <- merge(final_meta_data_wasp, pc20_ctl_fnf, by.x="ID_wasp" ,by.y="ID", all=FALSE) # Exclude last argument if it's "wasp"
} else {
  splicingPCA_df <- merge(final_meta_data,pc20_ctl_fnf, by="ID" , all=FALSE) # Include all remaining arguments otherwise
}


#Change to factor

# Assuming splicingPCA_df is your dataframe
for (col in colnames(splicingPCA_df)) {
  # Convert only columns starting with "PC" to numeric
  if (!grepl("^PC", col)) {
    splicingPCA_df[[col]] <-  as.factor(splicingPCA_df[[col]])
  }
}


#Seperate out the condition
# Separate the dataframe into two based on the Condition column

colnames(splicingPCA_df)[which(colnames(splicingPCA_df) == 'RNAshippedDate')] <- c("SeqeuncingBatch")
splicingPCA_numeric <- sapply(splicingPCA_df[,2:24],as.numeric) |> as.data.frame() 
splicingPCA_num <- cbind(splicingPCA_df[,1],splicingPCA_numeric)
# Convert columns to numeric, assuming first column is ID

if(use_wasp_id) {
  splicingPCA_num$ID_wasp <- gsub("_wasp","",splicingPCA_num$ID_wasp)
  colnames(splicingPCA_num)[1] <- "ID"
} 


# read in data
pca_pbs <- fread("./output/geno/pbs_geno/05.pca/CQTL.eigenvec")
pca_fnf <- fread("./output/geno/fnf_geno/05.pca/CQTL.eigenvec")

pca_pbs <- pca_pbs[,-1]
# set names
names(pca_pbs)[1] <- "ID"
names(pca_pbs)[2:ncol(pca_pbs)] <- paste0("GenoPC", 1:(ncol(pca_pbs)-1))

pca_fnf <- pca_fnf[,-1]
# set names
names(pca_fnf)[1] <- "ID"
names(pca_fnf)[2:ncol(pca_fnf)] <- paste0("GenoPC", 1:(ncol(pca_fnf)-1))


geno_pcs <-  rbind(pca_pbs,pca_fnf)
geno_pcs_df <- geno_pcs[,1:11] |> as.data.frame()
 
for (i in 1:20) {
  print(i)
  col_n <- i
  cov_merged <- merge(splicingPCA_num[,1:(i+4)],geno_pcs_df, by='ID',all.x=TRUE)
  cov_transposed <- t(cov_merged[, -1])
  colnames(cov_transposed) <- cov_merged$ID
  cov_fi <- cbind(rownames(cov_transposed), cov_transposed)
  colnames(cov_fi)[1] <- 'id' 
  write.table(cov_fi,file=paste0(bed_dir,"/","covariates_PC",col_n),sep='\t',quote=F,row.names=F,col.names=T)
}

