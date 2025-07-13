## 01. Merge two VCF file
setwd("/work/users/s/e/seyoun/CQTL_sQTL/output/geno")
#Reindex vcf file
tabix -p vcf /work/users/s/e/seyoun/CQTL_sQTL/output/geno/pbs_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz
tabix -p vcf /work/users/s/e/seyoun/CQTL_sQTL/output/geno/fnf_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz
# mkdir merged (This will only output sites that are identical between the two files, effectively skipping any conflicting sites.)
bcftools merge --merge none /work/users/s/e/seyoun/CQTL_sQTL/output/geno/pbs_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz /work/users/s/e/seyoun/CQTL_sQTL/output/geno/fnf_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz -o /work/users/s/e/seyoun/CQTL_sQTL/output/geno/merged/merged_pbs_fnf_output.vcf.gz

# view the final vcf file 
bcftools stats  merged_pbs_fnf_output.vcf.gz > merged_pbs_fnf_output_stats.txt

## 02. Make the new bed format for the adding the condition into the covariate format and run through pc1-pc20
## I need to modify /work/users/s/e/seyoun/CQTL_sQTL/scripts/sQTL_rscripts/bedfile_seperatebyCHR.R
## Saved the new covariates for the covariate_adding_condition_to_cov.R

## 03. Re-run the qtltools pc1-pc20 norm, perm. 
## sh qtlrun.sh --> to run the sbatch qtl_re_run_new_cov_condition.sbatch $i
## 03.-1then merge to one sh /work/users/s/e/seyoun/CQTL_sQTL/output/02.qtltools_condition_cov/qtltools_merge.sh

## 04. Count and compared

setwd("/work/users/s/e/seyoun/CQTL_sQTL")
library(tidyr)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(ensembldb)
library(AnnotationHub)
library(AnnotationFilter)
library(org.Hs.eg.db)
library(ggtext)

# Define the possible arguments files, significant p-values etc
adjusted_beta_p <- 0.05
qtl_dir <- "output/02.qtltools_condition_cov"

#-------------------------------------------------------------------------------
process_pbs_files <- function(filepath_pattern, pattern_suffix) {
  filepaths <- list.files(filepath_pattern, pattern = pattern_suffix, full.names = TRUE)   # List files based on the pattern
  
  # Sort filepaths numerically based on the PC number
  sorted_filepaths <- filepaths[order(as.numeric(gsub(".*pc(\\d+)_.*", "\\1", filepaths)))]
  
  # Read files and process them
  list_processed <- lapply(sorted_filepaths, function(filepath) {
    df <- read.table(filepath, header = FALSE, stringsAsFactors = FALSE)
    
    # Filter out rows with "chrX" or "chrY" from the specified column
    df_filtered <- df[!(df[, 2] %in% c("chrX", "chrY")), ]
    
    # Add adjusted p-values using FDR
    df_filtered$p_adjusted <- p.adjust(df_filtered[, 20], method = "fdr")
    
    # Filter based on significance threshold
    df_sig <- df_filtered[df_filtered$p_adjusted <= adjusted_beta_p, ]
    
    return(df_sig)
  })
  
  return(list_processed)
}


perm_processed <- process_pbs_files(qtl_dir, "_allchr.perm")


#dotplot to see the counts of snps only for the --------------------------------
# Count unique values in column 1 of each list element
unique_counts <- lapply(perm_processed, function(x) {
  length(unique(x$V1))-1
})

# Create a data frame with the list elements and names
df <- data.frame(pc = paste0("pc", 1:20), permuation_counts = unlist(unique_counts))
df$pc <- factor(df$pc, levels = paste0("pc", 1:20))

# Create the dot plot-----------------------------------------------------------
counts_dot_plot  <- ggplot(df, aes(x = pc, y = permuation_counts)) +
  geom_point(size = 2, color = "#2057A7") +
  labs(y = "Permutation counts") +
  scale_y_continuous(limits = c(6000, 10000)) +
  theme_classic() +
  geom_text(data = df[which.max(df$permuation_counts),], 
            aes(label = paste0(pc, ": ", permuation_counts)),
            vjust = -1, hjust = 0.5, size = 3.5, color = "#2057A7") +
  theme(text = element_text(family = "Helvetica"),
        panel.background = element_rect(fill = 'transparent', color = "transparent"),
        plot.background = element_rect(fill = 'transparent', color = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.background = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 9),
        axis.line = element_line(color = "black", linewidth = 0.25),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.ticks.y = element_line(color = "black", linewidth = 0.25),
        axis.text = element_text(color = "black", size = 8),
        axis.ticks.x = element_line(color = "black", linewidth = 0.25),
        axis.text.x = element_text(color = "black", size = 8, angle = 45, hjust = 1),
        plot.margin = margin(t = 5, r = 10, b = 5, l = 5))

##05. After the count with the max. This is for the computational analysis but also overall counts and the thresholds
# Rscript runFDR_cis.R /work/users/s/e/seyoun/CQTL_sQTL/output/02.qtltools_condition_cov/perm/pc1_allchr.perm 0.05 /work/users/s/e/seyoun/CQTL_sQTL/output/02.qtltools_condition_cov/significant_threshold/perm_0.05_pc1 

##06. Comparison of the current model. what are the new and what are the overlaps. compare with the response QTL?

perm_pc1 <- read.table("output/02.qtltools_condition_cov/perm/pc1_allchr.perm", header = FALSE, stringsAsFactors = FALSE)
header <- read.table("scripts/sQTL_rscripts/header.txt")

# Identify the position of the 'slope' column
slope_pos <- which(header == "slope")

# Create the new header by inserting 'slope_se' after 'slope' and adding 'qval' and 'threshold' at the end
new_header <- c(
  as.character(header[1:slope_pos]),
  "slope_se",
  as.character(header[(slope_pos+1):ncol(header)]),
  "qval",
  "threshold"
)

sigQTL_cond <- read.table("output/02.qtltools_condition_cov/significant_threshold/perm_0.05_pc1.significant.txt") |> as.data.frame()
colnames(sigQTL_cond) <- new_header

response_pbs_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_pbs_resInv.rds")
response_fnf_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_fnf_resInv.rds")

#Step1 Take out the conditional analysis to have the fair comparison

primary_pbs_sqtl <- response_pbs_results |> dplyr::filter(rank == "0") |> 
  dplyr::filter(!is.na(genomicLoc))
primary_fnf_sqtl <- response_fnf_results |> dplyr::filter(rank == "0") |> 
  dplyr::filter(!is.na(genomicLoc))

#Step2 Venn diagram to see how many sqtl-sIntron_junction overlaps in the response QTL

# Venndiagram for the response QTL 

# Subset significant reponse sQTL only
pbs_resQTL <- primary_pbs_sqtl %>%
  dplyr::filter(interaction_pval  < 0.05) |>
  arrange(interaction_pval)

fnf_resQTL <- primary_fnf_sqtl %>%
  dplyr::filter(interaction_pval  < 0.05) |>
  arrange(interaction_pval)

# Subset high confidence re-sQTls ####

pbs_highConf_resQtL <- primary_pbs_sqtl %>%
  dplyr::filter(interaction_pval  < 0.05) %>%
  dplyr::filter(abs(delta_beta ) > 0.2) |>
  dplyr::filter(minor_alle_count >= 5) |>
  arrange(interaction_pval)

fnf_highConf_resQtL <- primary_fnf_sqtl %>%
  dplyr::filter(interaction_pval  < 0.05) %>%
  dplyr::filter(abs(delta_beta ) > 0.2) |>
  dplyr::filter(minor_alle_count >= 5) |>
  arrange(interaction_pval)

pbs_specific_resQTL <- pbs_highConf_resQtL |>
  dplyr::filter(!clusterID %in% fnf_highConf_resQtL$clusterID )

fnf_specific_resQTL <- fnf_highConf_resQtL |>
  dplyr::filter(!clusterID %in% pbs_highConf_resQtL$clusterID )


# Load required libraries
library(ggplot2)
library(ggforce)
library(ggvenn)

# Function to extract unique QTL pairs
extract_qtl_pairs <- function(data) {
  unique(paste(data$phe_id, data$var_id, sep = "_"))
}

# Extract unique QTL pairs for each dataset
pbs_qtls <- extract_qtl_pairs(pbs_resQTL)
fnf_qtls <- extract_qtl_pairs(fnf_resQTL)
cond_qtls <- extract_qtl_pairs(sigQTL_cond)

# Create a list of the three sets
qtl_sets <- list(PBS_resQTL = pbs_qtls, FNF_resQTL = fnf_qtls, sigQTL_cond = cond_qtls)


ggvenn(
  qtl_sets, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4
)


# What about the overall sqtl finding?
# Extract unique QTL pairs for each dataset
pbs_qtls <- extract_qtl_pairs(primary_pbs_sqtl)
fnf_qtls <- extract_qtl_pairs(primary_fnf_sqtl)
cond_qtls <- extract_qtl_pairs(sigQTL_cond)

# Create a list of the three sets
sQTL_primary_sets <- list(PBS_sQTL = primary_pbs_sqtl$phe_id, FNF_sQTL = primary_fnf_sqtl$phe_id, sigQTL_conditionCov = sigQTL_cond$phe_id)

ggvenn(
  sQTL_primary_sets, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.2, set_name_size = 4
)

# Curious question. Why the 553 for the PBS and the FNF wasn't overlap with the sQTL for the condition?
get_exclusive <- function(set, other_sets) {
  setdiff(set, unlist(other_sets))
}

# Get exclusive elements for each set
PBS_exclusive <- get_exclusive(sQTL_primary_sets$PBS_sQTL, 
                               list(sQTL_primary_sets$FNF_sQTL, sQTL_primary_sets$sigQTL_conditionCov))
FNF_exclusive <- get_exclusive(sQTL_primary_sets$FNF_sQTL, 
                               list(sQTL_primary_sets$PBS_sQTL, sQTL_primary_sets$sigQTL_conditionCov))
sigQTL_exclusive <- get_exclusive(sQTL_primary_sets$sigQTL_conditionCov, 
                                  list(sQTL_primary_sets$PBS_sQTL, sQTL_primary_sets$FNF_sQTL))

pbs_exclusive_df <- primary_pbs_sqtl |> dplyr::filter(phe_id %in% PBS_exclusive) 
fnf_exclusive_df <- primary_fnf_sqtl |> dplyr::filter(phe_id %in% FNF_exclusive) 

#Make a plot to see it? for some cases. 
## plotted it using the utils.R createplot and also the enrichement directory meta_information.R for the main_figure4.R and create_intron_usage_plot() function. 
fnf_exclusive_df |> dplyr::filter(interaction_pval > 0.05) |> arrange(FNF_p)
PBS_exclusive_df |> dplyr::filter(interaction_pval > 0.05) |> arrange(PBS_p)
#-------------------------------------------------------------------------------
# Use the invNorm_expression of Nicole's


# I need to calculate PC for the 
pbs_invNorm_geneAnno <- fread("/work/users/s/e/seyoun/dbGap/nicole/Expression/PBS_CPMadjTMM_invNorm_geneAnnotations.txt")
pbs_invNorm <- fread("/work/users/s/e/seyoun/dbGap/nicole/Expression/PBS_CPMadjTMM_invNorm_expression.tsv")

fnf_invNorm_geneAnno <- fread("/work/users/s/e/seyoun/dbGap/nicole/Expression/FNF_CPMadjTMM_invNorm_geneAnnotations.txt")
fnf_invNorm <- fread("/work/users/s/e/seyoun/dbGap/nicole/Expression/FNF_CPMadjTMM_invNorm_expression.tsv")
