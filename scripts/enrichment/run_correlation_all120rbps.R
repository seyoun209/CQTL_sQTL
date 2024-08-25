library(GenomicRanges)
library(data.table)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(broom)
library(stringr)
library(tidyverse)
library(ggtext)

setwd("/work/users/s/e/seyoun/CQTL_sQTL")

response_pbs_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_pbs_resInv.rds")
response_fnf_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_fnf_resInv.rds")

load("external_data/encode_rbp/rbp_prep_rbpOnly/pbs_rbp_elip.Rdata") #pbs_rbp_eclip 
load("external_data/encode_rbp/rbp_prep_rbpOnly/fnf_rbp_elip.Rdata") #fnf_rbp_eclip

pbs_rbp_eclip_subset <- pbs_rbp_eclip |> group_by(rbp) %>%
  dplyr::filter(observed > 5) |>
  ungroup()  |>
  arrange(desc(odd_med))

fnf_rbp_eclip_subset <- fnf_rbp_eclip |> group_by(rbp) %>%
  dplyr::filter(observed > 5) |>
  ungroup()  |>
  arrange(desc(odd_med))

combined_data_rbp_eclip_sig <- rbind( pbs_rbp_eclip_subset,fnf_rbp_eclip_subset )
all_sig_enriched_rbps  <- unique(combined_data_rbp_eclip_sig$rbp)

# First, subset the all sig enriched rbps bed file and make grange to find the overlaps with the sQTLs. 

#qtl
pbs_sQtl_bed <- fread("output/Enrichment/rbp/data_prep/significant_pbs_rank0.bed")
colnames(pbs_sQtl_bed) <- c("var_chr", "var_start", "var_end", "var_id", "phe_id", "phe_strd")
fnf_sQtl_bed <- fread("output/Enrichment/rbp/data_prep/significant_fnf_rank0.bed")
colnames(fnf_sQtl_bed) <-  c("var_chr", "var_start", "var_end", "var_id", "phe_id", "phe_strd")

# Convert sQTL data to GRanges objects
pbs_gr <- with(pbs_sQtl_bed, GRanges(var_chr, IRanges(var_start, var_end)))
fnf_gr <- with(fnf_sQtl_bed, GRanges(var_chr, IRanges(var_start, var_end)))


load("output/combined_meta_data.RData")
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
meta_cqtl <- fread("output/clu_fnf/meta_cqtl")
meta_catl_simple <- meta_cqtl %>%
  dplyr::filter(Condition %in% c("CTL","FNF")) %>%
  dplyr::select(-Sex,-Age,-FragmentBatch,-RIN,-RNAextractionKitBatch,-RNAshippedDate,-Predicted_Ancestry) %>%
  dplyr::filter(!str_detect(Donor, 'AM7352|AM7244'))


#This is ratios and meta data
ctl_fnf_ratio <- fread("output/clu_fnf/ratio_fnf.txt") |> as.data.frame()
rownames(ctl_fnf_ratio) <- ctl_fnf_ratio$Junction
ratios_fnf <- ctl_fnf_ratio[,-1]
intron_test <- c(response_pbs_results$phe_id, response_fnf_results$phe_id) |> unique()
psi_matrix <- ratios_fnf[rownames(ratios_fnf) %in% intron_test,] |> as.data.frame()
intronID <- rownames(psi_matrix)
psi_ratio <- apply(psi_matrix, 2, as.numeric)
rownames(psi_ratio) <- intronID
psi_ratio.df <- t(psi_ratio) |> as.data.frame()

#Function
add_gene_symbols <- function(vsd_genes_exp, org_db = org.Hs.eg.db) {
  
  # Add gene symbols to DE results
  annots <- AnnotationDbi::select(org_db, 
                                  vsd_genes_exp$gene_id,
                                  columns = c("SYMBOL"), 
                                  keytype = "ENSEMBL")
  annots <- dplyr::rename(annots, gene_id = ENSEMBL)
  annots_unique <- annots %>% distinct(gene_id, .keep_all = TRUE)
  de_results_with_symbol <- left_join(vsd_genes_exp, annots_unique, by = "gene_id")
  
  
  
  return(de_results_with_symbol)
}


vsd_gene <- fread("output/quant/normalized_vst_gene.txt") |>
  dplyr::mutate(gene_id = gsub("\\..*$", "" , ENSG))

vsd_geneExp <- add_gene_symbols(vsd_gene)

#rbp_bed

rbp_sig.list <- list()
for(i in all_sig_enriched_rbps){
  print(i)
  rbp_bed <- fread(paste0("external_data/encode_rbp/rbp_prep_rbpOnly/",i,".bed"))
  colnames(rbp_bed) <- c("chr","start","end","strand","rbp_nm","cell_line","rep_n","-log2FC","-log10pval")
  rbp_sig.list[[i]] <- rbp_bed
}


# Function to find overlaps and return a data frame
find_overlaps <- function(sqtl_gr, rbp_bed) {
  rbp_gr <- with(rbp_bed, GRanges(chr, IRanges(start, end)))
  overlaps <- findOverlaps(sqtl_gr, rbp_gr)
  data.frame(
    sqtl_index = queryHits(overlaps),
    rbp_index = subjectHits(overlaps)
  )
}


pbs_overlaps <- lapply(rbp_sig.list, function(rbp) find_overlaps(pbs_gr, rbp))
fnf_overlaps <- lapply(rbp_sig.list, function(rbp) find_overlaps(fnf_gr, rbp))

#chr16: 4624826 4690972 
#4683932	4688796

#wwp2
#chr16 69925485 69925980
#mgrn1 <- GRanges(seqnames = "chr16", IRanges(4661789,4661790))

#wwp2 <- GRanges(seqnames = "chr16", IRanges(69921901,69921902))

#find_overlaps_mgrn1 <- find_overlaps(mgrn1, rbp_all.df)
#find_overlaps_wwp2 <- find_overlaps(wwp2, rbp_all.df)

#Now, 
create_overlap_matrix <- function(overlaps, sqtl_bed) {
  rbp_names <- names(overlaps)
  matrix_data <- sapply(rbp_names, function(rbp) {
    overlap_data <- overlaps[[rbp]]
    sqtl_bed$phe_id[overlap_data$sqtl_index]
  })
  
  
  for (i in seq_along(rbp_names)) {
    psi_matrix <- ctl_fnf_ratio |> dplyr::filter(Junction %in% matrix_data[[i]])
  }
  
  psi_matrix
}

pbs_matrix <- create_overlap_matrix(pbs_overlaps, pbs_sQtl_bed)
fnf_matrix <- create_overlap_matrix(fnf_overlaps, fnf_sQtl_bed)

#-------------------------------------------------------------------------------
#Trying to get all the significant RBP from 120 RBPs. 

#Getting combined to one that which rbp overlaps with
# create_overlap_matrix <- function(qtl_rbp_overlaps, sqtl_bed,
#                                   ctl_fnf_ratio, vsd_geneExp, response_results, 
#                                   target_rbps = c("EFTUD2", "SRSF9", "PRPF8", "SF3B4", "DROSHA"), 
#                                   significant_rbps = c("UTP3", "EFTUD2", "SRSF9", "PRPF8", "SF3B4", "NSUN2", "DROSHA"),
#                                   r2_threshold = 0.2, p_value_threshold = 0.05) {
#   
#   
#   # Combine PBS and FNF overlaps
#   all_overlaps <- c(qtl_rbp_overlaps)
#   all_sqtl_bed <- rbind(sqtl_bed)
#   
#   # Get variants overlapping with target RBPs
#   var_id_counts <- all_overlaps %>%
#     keep(names(.) %in% target_rbps) %>%
#     imap_dfr(~ tibble(
#       rbp = .y,
#       var_id = all_sqtl_bed$var_id[.x$sqtl_index],
#       phe_id = all_sqtl_bed$phe_id[.x$sqtl_index]
#     )) %>%
#     group_by(var_id, phe_id) %>%
#     summarise(
#       rbp_count = n_distinct(rbp),
#       overlapping_rbps = list(unique(rbp)),
#       .groups = "drop"
#     ) %>%
#     arrange(desc(rbp_count))
#   
#   # Get PSI data
#   psi_data <- ctl_fnf_ratio %>%
#     pivot_longer(cols = -Junction, names_to = "sample", values_to = "PSI") %>%
#     filter(Junction %in% var_id_counts$phe_id)
#   
#   # Get RBP expression data
#   rbp_exp_data <- vsd_geneExp %>%
#     filter(SYMBOL %in% target_rbps) %>%
#     pivot_longer(cols = -c(ENSG, gene_id, SYMBOL), names_to = "sample", values_to = "expression")
#   
#   # Combine PSI and RBP expression data
#   combined_data <- psi_data %>%
#     left_join(rbp_exp_data, by = "sample")
#   
#   # Calculate R-squared using linear regression for each variant-RBP pair
#   correlations <- var_id_counts %>%
#     rowwise() %>%
#     mutate(
#       correlations = list(
#         map_dfr(overlapping_rbps, function(rbp) {
#           data <- combined_data %>%
#             filter(Junction == phe_id, SYMBOL == rbp)
#           if(nrow(data) > 0) {
#             model <- lm(PSI ~ expression, data = data)
#             tibble(
#               SYMBOL = rbp,
#               r_squared = summary(model)$r.squared,
#               p_value = summary(model)$coefficients[2,4]
#             )
#           } else {
#             tibble(SYMBOL = rbp, r_squared = NA, p_value = NA)
#           }
#         })
#       )
#     ) %>%
#     unnest(correlations) %>%
#     filter(r_squared >= r2_threshold, p_value <= p_value_threshold, !is.na(r_squared)) %>% 
#     dplyr::rename(rbp_name = SYMBOL)
#   
#   # Join with response results and format overlapping_rbps
#   result <- correlations %>%
#     left_join(response_results, by = c("var_id", "phe_id")) %>%
#     group_by(var_id, phe_id) %>%
#     summarise(
#       rbp_count = n(),
#       overlapping_rbps = paste(rbp_name, collapse = ", "),
#       significant_rbps = paste(intersect(rbp_name, significant_rbps), collapse = ", "),
#       .groups = "drop"
#     ) %>%
#     left_join(response_results, by = c("var_id", "phe_id")) %>%
#     arrange(desc(rbp_count))
#   
#   return(result)
# }
# 
# 
# 
# 
# 
# 
# #This will be all of them 
# fnf_results_120rbp <- create_overlap_matrix(fnf_overlaps, fnf_sQtl_bed,  ctl_fnf_ratio, vsd_geneExp,  response_fnf_results, 
#                                             target_rbps = all_sig_enriched_rbps,
#                                             significant_rbps = c("UTP3", "EFTUD2", "SRSF9", "PRPF8", "SF3B4", "NSUN2", "DROSHA"),
#                                             r2_threshold = 0.2, p_value_threshold = 0.05)
# 
# save(fnf_results_120rbp, file="output/Enrichment/fnf_results_120rbp_overlapping.rda")
# 
# pbs_results_120rbp <- create_overlap_matrix(pbs_overlaps, pbs_sQtl_bed, 
#                                             ctl_fnf_ratio, vsd_geneExp, response_pbs_results, 
#                                             target_rbps = all_sig_enriched_rbps, 
#                                             significant_rbps = c("EFTUD2", "PRPF8", "SF3B4","DDX55","TROVE2","SATU2","U2AF2","TRA2A","CSTF2","SRSF7","SND1","AATF"),
#                                             r2_threshold = 0.2, p_value_threshold = 0.05)
# save(pbs_results_120rbp, file="output/Enrichment/pbs_results_120rbp_overlapping.rda")





#-------------------------------------------------------------------------------
#Total rbps 
create_overlap_matrix <- function(qtl_rbp_overlaps, sqtl_bed,
                                  ctl_fnf_ratio, vsd_geneExp, response_results, 
                                  target_rbps = c("EFTUD2", "SRSF9", "PRPF8", "SF3B4", "DROSHA"), 
                                  significant_rbps = c("UTP3", "EFTUD2", "SRSF9", "PRPF8", "SF3B4", "NSUN2", "DROSHA"),
                                  r2_threshold = 0.2, p_value_threshold = 0.05) {
  
  # Combine PBS and FNF overlaps
  all_overlaps <- c(qtl_rbp_overlaps)
  all_sqtl_bed <- rbind(sqtl_bed)
  
  # Get variants overlapping with target RBPs
  var_id_counts <- all_overlaps %>%
    keep(names(.) %in% target_rbps) %>%
    imap_dfr(~ tibble(
      rbp = .y,
      var_id = all_sqtl_bed$var_id[.x$sqtl_index],
      phe_id = all_sqtl_bed$phe_id[.x$sqtl_index]
    )) %>%
    group_by(var_id, phe_id) %>%
    summarise(
      rbp_count = n_distinct(rbp),
      overlapping_rbps = list(unique(rbp)),
      .groups = "drop"
    ) %>%
    arrange(desc(rbp_count))
  
  # Get PSI data
  psi_data <- ctl_fnf_ratio %>%
    pivot_longer(cols = -Junction, names_to = "sample", values_to = "PSI") %>%
    filter(Junction %in% var_id_counts$phe_id)
  
  # Get RBP expression data
  rbp_exp_data <- vsd_geneExp %>%
    filter(SYMBOL %in% target_rbps) %>%
    pivot_longer(cols = -c(ENSG, gene_id, SYMBOL), names_to = "sample", values_to = "expression")
  
  # Combine PSI and RBP expression data
  combined_data <- psi_data %>%
    left_join(rbp_exp_data, by = "sample")
  
  correlations <- var_id_counts %>%
    rowwise() %>%
    mutate(
      correlations = list(
        map_dfr(overlapping_rbps, function(rbp) {
          data <- combined_data %>%
            filter(Junction == phe_id, SYMBOL == rbp)
          if(nrow(data) > 0) {
            pearson_result <- cor.test(data$PSI, data$expression, method = "pearson")
            tibble(
              SYMBOL = rbp,
              pearson_r = pearson_result$estimate,
              r_squared = pearson_result$estimate^2,
              p_value = pearson_result$p.value
            )
          } else {
            tibble(SYMBOL = rbp, pearson_r = NA, r_squared = NA, p_value = NA)
          }
        })
      )
    ) %>%
    unnest(correlations) %>%
    dplyr::rename(rbp_name = SYMBOL)
  
  # Create a function to format RBP info
  format_rbp_info <- function(rbp_name, pearson_r, r_squared, p_value) {
    sprintf("%s (r=%.3f, R²=%.3f, p=%.3e)", rbp_name, pearson_r, r_squared, p_value)
  }
  
  # Join with response results and format RBP information
  result <- correlations %>%
    group_by(var_id, phe_id) %>%
    summarise(
      rbp_count = n(),
      all_overlapping_rbps = paste(format_rbp_info(rbp_name, pearson_r, r_squared, p_value), collapse = "; "),
      significant_overlapping_rbps = paste(format_rbp_info(
        rbp_name[r_squared >= r2_threshold & p_value <= p_value_threshold], 
        pearson_r[r_squared >= r2_threshold & p_value <= p_value_threshold],
        r_squared[r_squared >= r2_threshold & p_value <= p_value_threshold], 
        p_value[r_squared >= r2_threshold & p_value <= p_value_threshold]
      ), collapse = "; "),
      enriched_rbps = paste(intersect(rbp_name, significant_rbps), collapse = ", "),
      .groups = "drop"
    ) %>%
    left_join(response_results, by = c("var_id", "phe_id")) %>%
    arrange(desc(rbp_count))
  
  return(result)
}

all_sig_enriched_rbps <- all_sig_enriched_rbps |> data.frame() |>   filter(all_sig_enriched_rbps != "LIN28B") %>%  # Remove LIN28B
  mutate(all_sig_enriched_rbps = ifelse(all_sig_enriched_rbps == "TROVE2", "RO60", all_sig_enriched_rbps)) |> pull(all_sig_enriched_rbps)

fnf_results_allCorrelation_info <- create_overlap_matrix(qtl_rbp_overlaps=fnf_overlaps, sqtl_bed= fnf_sQtl_bed, ctl_fnf_ratio= ctl_fnf_ratio, 
                                                         vsd_geneExp=vsd_geneExp,  response_results= response_fnf_results, 
                                            target_rbps = all_sig_enriched_rbps,
                                            significant_rbps = c("UTP3", "EFTUD2", "SRSF9", "PRPF8", "SF3B4", "NSUN2", "DROSHA"),
                                            r2_threshold = 0.15, p_value_threshold = 0.05)


save(fnf_results_allCorrelation_info, file="output/Enrichment/fnf_results_allCorrelation_info.rda")

print("Done_fnf")

pbs_results_allCorrelation_info <- create_overlap_matrix(pbs_overlaps, pbs_sQtl_bed,
                                            ctl_fnf_ratio, vsd_geneExp, response_pbs_results,
                                            target_rbps = all_sig_enriched_rbps,
                                            significant_rbps = c("EFTUD2", "PRPF8", "SF3B4","DDX55","TROVE2","SATU2","U2AF2","TRA2A","CSTF2","SRSF7","SND1","AATF"),
                                            r2_threshold = 0.15, p_value_threshold = 0.05)
save(pbs_results_allCorrelation_info, file="output/Enrichment/pbs_results_allCorrelation_info.rda")


print("Done_pbs")