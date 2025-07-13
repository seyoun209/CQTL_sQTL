# Make qtl file
setwd("/work/users/s/e/seyoun/CQTL_sQTL/")
library(dplyr)
library(data.table)
library(ggplot2)
library(readr)
library(tidyr)
library(purrr)

#response_pbs_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_pbs_results.rds") 
#response_fnf_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_fnf_results.rds")
response_pbs_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_pbs_resInv.rds")
response_fnf_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_fnf_resInv.rds")


pbs_bed_temp <- response_pbs_results |> dplyr::filter(rank == 0) |> dplyr::select("var_chr","var_from", "var_id", "phe_id", "phe_strd")
pbs_bed <- pbs_bed_temp |> dplyr::mutate(var_start  =  var_from -1) |> dplyr::rename(var_end = "var_from") |>
  dplyr::select(var_chr, var_start, var_end, var_id, phe_id, phe_strd)
write.table(pbs_bed, file="output/Enrichment/rbp/data_prep/significant_pbs_rank0.bed",row.names = F,col.names = F,quote = F,sep = "\t")

# Prepare a bed for phenotype
pbs_pheno_temp <- response_pbs_results |> dplyr::filter(rank == 0) |> 
  dplyr::select("phe_chr","phe_from", "phe_to", "phe_id", "var_id","phe_strd")

pbs_pheno_bed <- pbs_pheno_temp |> dplyr::mutate(phe_start  =  phe_from -1) |>
  dplyr::select("phe_chr","phe_start", "phe_to", "phe_id", "var_id","phe_strd")
write.table(pbs_pheno_bed, file="output/Enrichment/rbp/data_prep/significant_pbs_quantified_rank0.bed",row.names = F,col.names = F,quote = F,sep = "\t")


#FNF ####

fnf_bed_temp <- response_fnf_results |> dplyr::filter(rank == 0) |> 
  dplyr::select("var_chr","var_from", "var_id", "phe_id", "phe_strd")
fnf_bed <- fnf_bed_temp |> dplyr::mutate(var_start  =  var_from -1) |> dplyr::rename(var_end = "var_from") |>
  dplyr::select(var_chr, var_start, var_end, var_id, phe_id, phe_strd)
write.table(fnf_bed, file="output/Enrichment/rbp/data_prep/significant_fnf_rank0.bed",row.names = F,col.names = F,quote = F,sep = "\t")

# Prepare a bed for phenotype
fnf_pheno_temp <- response_fnf_results |> dplyr::filter(rank == 0)|> 
  dplyr::select("phe_chr","phe_from", "phe_to", "phe_id", "var_id","phe_strd")

fnf_pheno_bed <- fnf_pheno_temp |> dplyr::mutate(phe_start  =  phe_from -1) |>
  dplyr::select("phe_chr","phe_start", "phe_to", "phe_id", "var_id","phe_strd")
write.table(fnf_pheno_bed, file="output/Enrichment/rbp/data_prep/significant_fnf_quantified_rank0.bed",row.names = F,col.names = F,quote = F,sep = "\t")


rbp_bed <- fread("tools/rbp_db.txt.gz")

 
write.table(rbp_bed, file="output/Enrichment/rbp/data_prep/rbp.bed",row.names = F,col.names = F,quote = F,sep = "\t")



#-------------------------------------------------------------------------------
#RBP 

#RBP sorted by 
rbp_db <- fread("output/Enrichment/rbp/data_prep/rbp.bed")
colnames(rbp_db) <- c("chr_binding","start","end","peakID","strand","RBPname","experiment_method",
                      "tissue","accession","confidence_score")


# Analyze experiment methods
experiment_summary <- rbp_db %>%
  group_by(experiment_method) %>%
  summarise(count = n()) %>%
  arrange(desc(count))


ggplot(experiment_summary, aes(x = reorder(experiment_method, count), y = count)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Distribution of Experiment Methods",
       x = "Experiment Method",
       y = "Count") +
  theme_minimal()


# Analyze tissue types
tissue_summary <- rbp_db %>%
  group_by(tissue) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

# Visualize top 20 tissue types
ggplot(tissue_summary %>% top_n(20, count), aes(x = reorder(tissue, count), y = count)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Distribution of Top 20 Tissue Types",
       x = "Tissue Type",
       y = "Count") +
  theme_minimal()

# Analyze confidence scores
summary(rbp_db$confidence_score)

ggplot(rbp_db, aes(x = confidence_score)) +
  geom_histogram(bins = 50) +
  labs(title = "Distribution of Confidence Scores",
       x = "Confidence Score",
       y = "Count") +
  theme_minimal()

# Analyze relationship between experiment method and confidence score
ggplot(rbp_db, aes(x = fct_reorder(experiment_method, confidence_score), y = confidence_score)) +
  geom_boxplot() +
  coord_flip() +
  labs(title = "Confidence Scores by Experiment Method",
       x = "Experiment Method",
       y = "Confidence Score") +
  theme_minimal()

# Print top RBPs by frequency
top_rbps <- rbp_db %>%
  group_by(RBPname) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  top_n(20)

#Get the Only eclip data sets


# I will divided First Brain and calculate and make bed file based on the RBP 

# Step 1: Remove PIP-seq experiments
rbp_db_filtered <- rbp_db %>%
  filter(experiment_method != "PIP-seq")


# Step 2: Function to create bed files for a given tissue
create_bed_files <- function(data, tissue_name) {
  data %>%
    dplyr::filter(tissue == tissue_name) %>%
    group_by(RBPname) %>%
    do({
      bed_data <- select(., "chr_binding","start","end","peakID","strand","RBPname","experiment_method",
                         "tissue","accession","confidence_score")
      file_name <- paste0("output/Enrichment/rbp/data_prep/", tissue_name, "_", .$RBPname[1], ".bed")
      write_tsv(bed_data, file_name, col_names = FALSE, quote = "none")
      tibble()
    })
}

# Step 3: Create bed files for each tissue
tissues_to_process <- c("HEK293T", "HEK293", "HeLa")
tissues_to_process <- c("Brain")


for (tissue in tissues_to_process) {
  create_bed_files(rbp_db_filtered, tissue)
  cat("Completed processing for", tissue, "\n")
}

cat("All bed files have been created.\n")

# RBP ALL DATA

# Function to create bed files for each RBP
create_rbp_bed_files <- function(data) {
  data %>%
    group_by(RBPname) %>%
    do({
      bed_data <- select(., chr_binding, start, end, peakID, strand, RBPname, experiment_method,
                         tissue, accession, confidence_score)
      file_name <- paste0("output/Enrichment/rbp_all/data_prep/", unique(.$RBPname), ".bed")
      write_tsv(bed_data, file_name, col_names = FALSE, quote = "none")
      tibble()
    })
  cat("RBP bed files created in output/Enrichment/rbp_all/data_prep/\n")
}

# Create RBP bed files
create_rbp_bed_files(rbp_db_filtered)


#-------------------------------------------------------------------------------
# Using only encode eCLIP data 

#-------------------------------------------------------------------------------
# chromHMM

chromhmm_temp <- fread("tools/chromhmm/E049_15_coreMarks_hg38lift_dense.bed.gz")
colnames(chromhmm_temp) <- c("chr","start","end","name","score","strand","ThickStart","ThickEnd","ItemRgb")

##15 states in ChromHMM data
states = c("1_TssA","2_TssAFlnk","3_TxFlnk","4_Tx","5_TxWk","6_EnhG","7_Enh","8_ZNF/Rpts","9_Het","10_TssBiv","11_BivFlnk","12_EnhBiv","13_ReprPC","14_ReprPCWk","15_Quies");

# Create output directory
output_dir <- "output/Enrichment/chromhmm/data_prep"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# Create output directory
output_dir <- "output/Enrichment/chromhmm/data_prep"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Process each state and save as bed file
for (state in states) {
  state_data <- chromhmm_temp %>%
    filter(name == state) %>%
    select(chr, start, end)
  
  # Replace "/" with "_" in the filename
  safe_state_name <- gsub("/", "_", state)
  
  output_file <- file.path(output_dir, paste0(safe_state_name, ".bed"))
  fwrite(state_data, output_file, sep = "\t", col.names = FALSE)
  
  cat("Saved", state, "to", output_file, "\n")
}




