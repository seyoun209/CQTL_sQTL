setwd("/work/users/s/e/seyoun/CQTL_sQTL")
library(data.table)
# Garfield
# make bed file for the Garfield input
# Bash script
## subset bedtools intersect for the RBP binding sites with the sQTl postion
#step1: make bed format for the response-sQTL

sig_response_pbs_nolmer <- fread("output/01.qtltools_re/response_qtl/reseponsesQTL_PBS_significant.csv")
sig_response_fnf_nolmer <- fread("output/01.qtltools_re/response_qtl/reseponsesQTL_FNF_significant.csv")
rbp_bed <- fread("/work/users/s/e/seyoun/CQTL_sQTL/tools/human.txt.gz")
colnames(rbp_bed) <- c("chr","start","end","peak_id","strand","RBP_name","experiment_method","sample_tissue","accession","confidence_score")
# For the sQTL enrichement, It is better to have high-resolutation and high-confidence. 

methods_table <- tibble(
  Method = c("4SU-iCLIP", "eCLIP", "HITS-CLIP", "iCLAP", "PAR-CLIP", "PIP-seq"),
  `Best Use Case` = c(
    "High-resolution, transient interactions",
    "Comprehensive, high-confidence mapping",
    "Exploratory, broad patterns",
    "Specific, high-affinity interactions",
    "High-resolution, high-confidence mapping",
    "Dynamic interactions across conditions"
  ),
  `Suggested Computational Tools` = c(
    "PureCLIP, Piranha",
    "Piranha, CLIPper",
    "CIMS, CLIPper",
    "PureCLIP, Piranha",
    "PARalyzer",
    "MiClip, Piranha"
  )
)

# Print the table
print(methods_table)

"eCli"


dir.create("output/garfield/sqtl/rsQTL_pbs",recursive = TRUE, showWarnings = FALSE)
dir.create("output/garfield/sqtl/rsQTL_fnf",recursive = TRUE, showWarnings = FALSE)
dir.create("output/garfield/rbp",recursive = TRUE, showWarnings = FALSE)

#Filter and save BED files for each chromosome
for (i in 1:22) {
  bed_data <- sig_response_pbs_nolmer %>%
    dplyr::filter(var_chr == paste0("chr", i)) %>%
    dplyr::select(var_chr, var_from, var_to)
  
  write.table(bed_data, file = paste0("output/garfield/sqtl/rsQTL_pbs/rsqtl_chr", i, ".bed"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

for (i in 1:22) {
  bed_data <- sig_response_fnf_nolmer %>%
    dplyr::filter(var_chr == paste0("chr", i)) %>%
    dplyr::select(var_chr, var_from, var_to)
  
  write.table(bed_data, file = paste0("output/garfield/sqtl/rsQTL_fnf/rsqtl_chr", i, ".bed"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

for (i in 1:22) {
  bed_data <- rbp_bed %>%
    dplyr::filter(var_chr == paste0("chr", i)) %>%
    dplyr::select(var_chr, var_from, var_to)
  
  write.table(bed_data, file = paste0("output/garfield/sqtl/rsQTL_fnf/rsqtl_chr", i, ".bed"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}


## annotation
# make bed file for the 


## maftssd

## pval

## tags

