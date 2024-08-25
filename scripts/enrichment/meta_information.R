library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(broom)
library(stringr)
library(tidyverse)
library(ggtext)
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

# This is the way to inverse the PSI value



inverseNormGene <- function(geneRow){
  normValues <- qnorm((rank(as.numeric(geneRow),
                            na.last = "keep") - 0.5)/sum(!is.na(as.numeric(geneRow))))
  return(normValues)
}

# Function to process data for a specific condition
process_condition <- function(data, condition) {
  # Select columns for the condition
  condition_cols <- grep(condition, colnames(data), value = TRUE)
  
  # Apply inverse normal transformation
  norm_data_subset <- data[,condition_cols] 
  inv_norm_data <- as.data.frame(t(apply(norm_data_subset, 1, inverseNormGene)))
  
  colnames(inv_norm_data) <- condition_cols
  
  # Combine with other information
  #result <- data %>%
  #  dplyr::select(seqnames, start, end, width, strand, gene_id, geneChr, geneStart, geneEnd, geneLength, geneStrand, SYMBOL) %>%
  #  bind_cols(inv_norm_data) %>%
  #  rename(chr = seqnames) %>%
  #  arrange(chr, start)
  
  return(inv_norm_data)
}

# Process for CTL condition
CTL_invNorm <- process_condition(psi_ratio, "CTL")

# Process for FNF condition
FNF_invNorm <- process_condition(psi_ratio, "FNF")
CTL_invNorm$id <- rownames(CTL_invNorm)
FNF_invNorm$id <- rownames(FNF_invNorm)

all_psi_inv <- left_join(CTL_invNorm,FNF_invNorm,by="id")
rownames(all_psi_inv) <- all_psi_inv$id
all_psi_inv <- all_psi_inv |> dplyr::select(-'id')


pbsGeno_raw <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/pbs_geno/06.subset_sigSNps/recodeA_pbs.raw")
colnames(pbsGeno_raw) <- sub("_.*", "", colnames(pbsGeno_raw))

fnfGeno_raw <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/fnf_geno/06.subset_sigSNps/recodeA_fnf.raw")
colnames(fnfGeno_raw) <- sub("_.*", "", colnames(fnfGeno_raw))
all_geno_matrix <- bind_rows(pbsGeno_raw,fnfGeno_raw)
all_geno_transpose_df <- all_geno_matrix %>%
  column_to_rownames(var = "IID")  %>%
  dplyr::select(-FID, -MAT,-PAT, -SEX, -PHENOTYPE) %>% 
  as.data.frame() %>%  # Ensure it's still a data frame
  t() 



# This is plotting for the boxplot for the sqtl and RBP only

create_intron_usage_plot <- function(gene_name, subset_phe_id, psi_ratio, all_geno_transpose_df, meta_catl_simple, response_results) {
  # Get the most significant intron junction
  mostSig_phe_id <- response_results |> 
    dplyr::filter(phe_id %in% subset_phe_id) |> 
    slice_min(FNF_p) |> 
    pull(phe_id)
  
  # Get minor allele and variant info
  minorAllele <- response_results |> 
    dplyr::filter(phe_id %in% subset_phe_id) |> 
    slice_min(FNF_p) |> 
    pull(minor_allele)
  
  variantID <- response_results |> 
    dplyr::filter(phe_id %in% subset_phe_id) |> 
    slice_min(FNF_p) |> 
    pull(var_id)
  
  rsID <- response_results |> 
    dplyr::filter(phe_id %in% subset_phe_id) |> 
    slice_min(FNF_p) |> 
    pull(rsID)
  
  # Process allele information
  alleles <- unlist(strsplit(variantID, ":"))
  ref_allele <- alleles[3]
  alt_allele <- alleles[4]
  protective_allele <- ifelse(minorAllele == ref_allele, alt_allele, ref_allele)
  
  # Get PSI data
  psi_matrix <- psi_ratio[rownames(psi_ratio) %in% mostSig_phe_id, ]
  psi_data <- data.frame(sampleID = names(psi_matrix), psi = as.double(psi_matrix))
  
  geno_data <- all_geno_transpose_df[rownames(all_geno_transpose_df) %in% variantID, ] %>%
    as.data.frame() %>%
    rownames_to_column(var = "sampleID")
  colnames(geno_data)[2] <- "genotype"
  
  meta_data_fi <- geno_data %>%
    left_join(meta_catl_simple,by = c("sampleID"="ID")) %>%
    left_join(psi_data,by = c("sampleID"="sampleID"))
  
  meta_combined_all <-  meta_data_fi %>%
    mutate(
      genotype = factor(genotype, levels = c("0", "1", "2")),
      Condition = ifelse(Condition == "CTL", "PBS", ifelse(Condition == "FNF", "FN-f", NA)),
      Condition = factor(Condition, levels = c("PBS", "FN-f")),
      genotype = case_when(
        genotype == "2" ~ paste(minorAllele, minorAllele, sep = "/"),
        genotype == "1" ~ paste(minorAllele, protective_allele, sep = "/"),
        genotype == "0" ~ paste(protective_allele, protective_allele, sep = "/")
      )
    )  %>%
    mutate(
      genotype = factor(genotype, levels = c(
        paste(protective_allele, protective_allele, sep = "/"),
        paste(minorAllele, protective_allele, sep = "/"),
        paste(minorAllele, minorAllele, sep = "/")
      ))
    )
  
  
  maxrange <- range(meta_combined_all$psi, na.rm = TRUE)
  yaxis_range_min <- min(0,round(maxrange[1],1))
  yaxis_range_max <- round(maxrange[2], 1)
  
  if(yaxis_range_min == 0 & yaxis_range_max == 0){
    yaxis_range_min <-  maxrange[1]
    yaxis_range_max <-  maxrange[2]
  }
  
  ggplot(meta_combined_all, aes(x = genotype, y = psi)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5)+
    geom_boxplot(outlier.shape = NA,
                 linewidth = 0.25, alpha = 0.7, fill = "grey80") +
    geom_jitter(color = "grey30",position = position_jitter(width = 0.2)  ,size = 0.25) +
    #geom_smooth(method = "lm", se = FALSE, color = "#FFB81C",aes(group = 1)) +
    labs(x = "Genotype", y = "Normalized PSI") +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(strip.placement = "outside",
          axis.line.y =  element_line(linewidth = 0.25),
          axis.line.x =  element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(color = "black", linewidth = 0.25),
          axis.ticks.length.y = unit(-0.1, "cm"),
          axis.title.x = element_blank(),
          axis.title.y = element_markdown(size = 6, family= "Helvetica",
                                          margin = margin(r = -0.5)),
          text = element_text(family = "Helvetica"),
          axis.text.y = element_text(color = "black", size = 6),
          axis.text.x = element_text(color = "black", size = 6, margin = margin(b = 0)),
          strip.background = element_blank(),
          strip.text = element_blank(),  # This removes the facet titles
          panel.background = element_rect(fill = "transparent", color = "transparent"),
          plot.background = element_rect(fill = "transparent", color = "transparent"),
          panel.grid = element_blank(),
          legend.position = "none",  # This removes the legend
          panel.spacing.x = unit(0.5, "cm")) 
}



snhg29_Genoboxplot <- create_intron_usage_plot("SNHG29", snhg29_phe_id, all_psi_inv, all_geno_transpose_df,meta_catl_simple,response_fnf_results)
create_intron_usage_plot("CALU", MICA_results$phe_id, all_psi_inv,all_geno_transpose_df,meta_catl_simple,response_fnf_results)


save(snhg29_Genoboxplot,file = "output/results_plots/rbp/rbp_snhg29_Genoboxplot.rda")
save(calu_Genoboxplot, file = "output/results_plots/rbp/rbp_calu_Genoboxplot.rda")
mica_Genoboxplot <- create_intron_usage_plot("MICA", MICA_results$phe_id,all_psi_inv ,all_geno_transpose_df,meta_catl_simple,response_pbs_results)
save(mica_Genoboxplot, file = "output/results_plots/rbp/mica_Genoboxplot.rda")

PCYT1A_Genoboxplot <- create_intron_usage_plot("PCYT1A", pcyt1a_results$phe_id[2],all_psi_inv ,all_geno_transpose_df,meta_catl_simple,response_fnf_results)

save(PCYT1A_Genoboxplot, file = "output/results_plots/rbp/PCYT1A_Genoboxplot.rda")

#Annotating genes Symbol for the normalzied ones


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

#Visualize RBP vs Intron junction

visualize_rbp_psi_correlation <- function(gene_name, subset_phe_id, psi_ratio, response_results, rbp_exp, vsd_genes_exp,p_threshold = 0.05,r2_threshold = 0.2) {

  
  # Get the most significant intron junction
  mostSig_phe_id <- response_results |> 
    dplyr::filter(phe_id %in% subset_phe_id) |> 
    slice_min(FNF_p) |> 
    pull(phe_id)
  
  # Get PSI data
  psi_matrix <- psi_ratio[rownames(psi_ratio) %in% mostSig_phe_id, ]
  psi_data <- data.frame(sampleID = names(psi_matrix), psi = as.double(psi_matrix))
  
  # Process RBP expression data
  rbp_exp_data <- vsd_geneExp |> 
    filter(SYMBOL %in% rbp_exp$SYMBOL) |> 
    pivot_longer(cols = -c(ENSG, gene_id,SYMBOL), names_to = "sampleID", values_to = "expression")
  
  # Combine PSI and expression data
  combined_data <- psi_data |> 
    left_join(rbp_exp_data, by = "sampleID")
  
  combined_data <- combined_data %>%
    mutate(Group = case_when(
      grepl("CTL", sampleID) ~ "PBS",
      grepl("FNF", sampleID) ~ "Fn-f",
      TRUE ~ "Unknown"
    ))
  
  # Calculate R² values and p-values for each group
  r_squared <- combined_data %>%
    group_by(SYMBOL, Group) %>%
    do({
      model <- lm(psi ~ expression, data = .)
      summary <- glance(model)
      data.frame(
        r.squared = summary$r.squared,
        p.value = summary$p.value
      )
    }) %>%
    mutate(label = sprintf("%s: R² = %.2f", Group, r.squared)) %>%
    arrange(desc(Group)) %>%  # This ensures PBS comes first
    group_by(SYMBOL) %>%
    mutate(y_position = 1 - (0.05 * (row_number() - 1)))
  
  # Filter for RBPs with at least one significant group and R² > threshold
  significant_rbps <- r_squared %>%
    group_by(SYMBOL) %>%
    summarise(
      any_significant = any(p.value < p_threshold),
      any_high_r2 = any(r.squared > r2_threshold)
    ) %>%
    filter(any_significant & any_high_r2) %>%
    pull(SYMBOL)
  
  # Filter the data for significant RBPs
  combined_data_filtered <- combined_data %>% filter(SYMBOL %in% significant_rbps)
  r_squared_filtered <- r_squared %>% filter(SYMBOL %in% significant_rbps)
  
  # Get y-axis limits for each facet
  y_limits <- combined_data_filtered %>%
    group_by(SYMBOL) %>%
    summarise(ymin = min(psi), ymax = max(psi))
  
  # Join y_limits to r_squared
  r_squared_filtered <- r_squared_filtered %>%
    left_join(y_limits, by = "SYMBOL") %>%
    mutate(y_label = ymin + (ymax - ymin) * y_position)
  
  ggplot(combined_data_filtered, aes(x = expression, y = psi, color = Group)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, aes(group = Group)) +
    facet_wrap(~ SYMBOL, scales = "free") +
    scale_color_manual(values = c("PBS" = "#1e87a5", "Fn-f" = "#FFB81C")) +
    geom_text(data = r_squared_filtered,
              aes(x = -Inf, y = y_label, label = label, color = Group),
              hjust = -0.1, vjust = 1, size = 3)  +  # Adjusted hjust to move text slightly left
    theme_minimal() +
    labs(x = "RBP Expression", y = "Intron Usage (PSI)") +
    theme(strip.placement = "outside",
          axis.line.y =  element_line(linewidth = 0.25),
          axis.line.x =  element_line(linewidth = 0.25),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(color = "black", linewidth = 0.25),
          axis.ticks.length.y = unit(-0.1, "cm"),
          axis.title.x = element_markdown(size = 6, family= "Helvetica",
                                          margin = margin(r = -0.5)),
          axis.title.y = element_markdown(size = 6, family= "Helvetica",
                                          margin = margin(r = -0.5)),
          text = element_text(family = "Helvetica"),
          axis.text.y = element_text(color = "black", size = 6),
          axis.text.x = element_text(color = "black", size = 6, margin = margin(b = 0)),
          strip.background = element_blank(),
          panel.background = element_rect(fill = "transparent", color = "transparent"),
          plot.background = element_rect(fill = "transparent", color = "transparent"),
          panel.grid = element_blank(),
          legend.position = "none",
          panel.spacing.x = unit(1, "cm"))
}

visualize_rbp_psi_correlation("serpina5", serpina5_phe_id, psi_ratio, response_fnf_results, serpina5_rbp_exp, vsd_geneExp,
                              p_threshold = 1, r2_threshold = 0)
visualize_rbp_psi_correlation("LRRFIP2", schm1_results$phe_id, psi_ratio, response_fnf_results, schm1_results$rbp_exp, vsd_geneExp,p_threshold = 0.05, r2_threshold = 0.2)

visualize_rbp_psi_correlation("CALU", calu_results$phe_id, psi_ratio, response_fnf_results, calu_results$rbp_exp, vsd_geneExp,p_threshold = 0.05, r2_threshold = 0.2)


visualize_rbp_psi_correlation("PCYT1A", MICA_results$phe_id, psi_ratio, response_pbs_results, MICA_results$rbp_exp, vsd_geneExp,
                              p_threshold = 0.05, r2_threshold = 0)

visualize_rbp_psi_correlation <- function(gene_name, subset_phe_id, psi_ratio_data, response_results, rbp_exp, vsd_genes_exp, 
                                          correlation_method = "spearman", threshold = 0.4, p_threshold = 0.01, selected_rbps = NULL) {
  # Get the most significant intron junction
  mostSig_phe_id <- response_results |> 
    dplyr::filter(phe_id %in% subset_phe_id) |> 
    slice_min(FNF_p) |> 
    pull(phe_id)
  
  # Get PSI data
  psi_matrix <- psi_ratio_data[rownames(psi_ratio_data) %in% mostSig_phe_id, ]
  psi_data <- data.frame(sampleID = names(psi_matrix), psi = as.double(psi_matrix))
  
  # Process RBP expression data
  rbp_exp_data <- vsd_genes_exp |> 
    filter(gene_id %in% rbp_exp$gene_id) |> 
    pivot_longer(cols = -c(ENSG, gene_id, SYMBOL), names_to = "sampleID", values_to = "expression")
  
  # Combine all data
  combined_data <- psi_data |> 
    left_join(rbp_exp_data, by = "sampleID")
  
  # Calculate correlation values and p-values
  model_stats <- combined_data %>%
    group_by(SYMBOL) %>%
    summarise(
      correlation = cor(psi, expression, method = correlation_method, use = "pairwise.complete.obs"),
      p.value = tryCatch({
        cor.test(psi, expression, method = correlation_method, exact = FALSE, continuity = TRUE)$p.value
      }, warning = function(w) NA_real_, error = function(e) NA_real_)
    ) %>%
    mutate(
      label = sprintf("%s = %.2f%s", 
                      ifelse(correlation_method == "spearman", "Rho", "R²"),
                      correlation, 
                      ifelse(!is.na(p.value) & p.value < p_threshold, "*", "")),
      y_position = 1
    )
  
  # Filter for selected RBPs if provided
  if (!is.null(selected_rbps)) {
    combined_data <- combined_data %>% filter(SYMBOL %in% selected_rbps)
    model_stats <- model_stats %>% filter(SYMBOL %in% selected_rbps)
  }
  
  # Filter for RBPs with correlation > threshold
  significant_rbps <- model_stats %>%
    filter(abs(correlation) >= threshold) %>%
    pull(SYMBOL)
  
  # Filter the data for significant RBPs
  combined_data_filtered <- combined_data %>% filter(SYMBOL %in% significant_rbps)
  model_stats_filtered <- model_stats %>% filter(SYMBOL %in% significant_rbps)
  
  # Create the plot
  ggplot(combined_data_filtered, aes(x = expression, y = psi)) +
    geom_point(alpha = 0.6, size = 1, color = "grey30") +
    geom_smooth(method = "lm", se = FALSE, color = "#FFB81C") +
    facet_wrap(~ SYMBOL, scales = "free", ncol = 3) +
    geom_text(data = model_stats_filtered,
              aes(x = -Inf, y = Inf, label = label),
              hjust = -0.1, vjust = 1, size = 3)  +
    labs(x = "RBP Expression", y = "Intron junction (PSI)") +
    theme_minimal() +
    theme(#strip.text = element_blank(),
      axis.line.y =  element_line(linewidth = 0.25),
      axis.line.x =  element_line(linewidth = 0.25),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(color = "black", linewidth = 0.25),
      axis.ticks.length.y = unit(-0.1, "cm"),
      axis.title.x = element_markdown(size = 6, family= "Helvetica", margin = margin(r = -0.5)),
      axis.title.y = element_markdown(size = 6, family= "Helvetica", margin = margin(r = -0.5)),
      text = element_text(family = "Helvetica"),
      axis.text.y = element_text(color = "black", size = 6),
      axis.text.x = element_text(color = "black", size = 6, margin = margin(b = 0)),
      strip.background = element_blank(),
      panel.background = element_rect(fill = "transparent", color = "transparent"),
      plot.background = element_rect(fill = "transparent", color = "transparent"),
      panel.grid = element_blank(),
      legend.position = "none",
      panel.spacing.x = unit(1, "cm")
    )
}


#correlation_method = "spearman" or correlation_method =  "pearson"
snhg29_rbpCorrPlot <- visualize_rbp_psi_correlation("SNHG29", snhg29_results$phe_id,psi_ratio, response_fnf_results, snhg29_results$rbp_exp, vsd_geneExp,correlation_method = "pearson",
                              threshold = 0.5, p_threshold = 0.05,selected_rbps=c('AATF'))
save(snhg29_rbpCorrPlot, file="output/results_plots/rbp/snhg29_rbpCorrPlot.rda")

mica_tra2a_rbpCorr_plot <- visualize_rbp_psi_correlation("MICA", MICA_results$phe_id, psi_ratio, response_pbs_results, MICA_results$rbp_exp, vsd_geneExp,correlation_method = "spearman",
                              threshold = 0.2, p_threshold = 0.05,selected_rbps=c("TRA2A"))
save(mica_tra2a_rbpCorr_plot, file="output/results_plots/rbp/mica_tra2a_rbpCorr_plot.rda")
mica_aatf_rbpCorr_plot <- visualize_rbp_psi_correlation("MICA", MICA_results$phe_id, psi_ratio, response_pbs_results, MICA_results$rbp_exp, vsd_geneExp,correlation_method = "spearman",
                              threshold = 0.2, p_threshold = 0.05,selected_rbps=c("AATF"))
save(mica_aatf_rbpCorr_plot, file="output/results_plots/rbp/mica_aatf_rbpCorr_plot.rda")
mica_cstf2_rbpCorr_plot <- visualize_rbp_psi_correlation("MICA", MICA_results$phe_id, psi_ratio, response_pbs_results, MICA_results$rbp_exp, vsd_geneExp,correlation_method = "spearman",
                              threshold = 0.2, p_threshold = 0.05,selected_rbps=c("CSTF2"))
save(mica_cstf2_rbpCorr_plot, file="output/results_plots/rbp/mica_cstf2_rbpCorr_plot.rda")
calu_rbpCorrPlot <- visualize_rbp_psi_correlation("CALU", calu_results$phe_id, psi_ratio, response_fnf_results, calu_results$rbp_exp, vsd_geneExp,correlation_method = "spearman",
                              threshold = 0.3, p_threshold = 0.05)

save(calu_rbpCorrPlot, file="output/results_plots/rbp/calu_rbaCorrPlot.rda")
PCYT1A_rbpCorrPlot <- visualize_rbp_psi_correlation("PCYT1A", pcyt1a_results$phe_id[2], psi_ratio, response_fnf_results, pcyt1a_results$rbp_exp, vsd_geneExp,correlation_method = "pearson",
                                                  threshold = 0.2, p_threshold = 0.05,selected_rbps =  c("EFTUD2"))

save(PCYT1A_rbpCorrPlot, file="output/results_plots/rbp/PCYT1A_rbpCorrPlot.rda")

visualize_rbp_psi_correlation("PCYT1A", MICA_results$phe_id, psi_ratio, response_fnf_results, MICA_results$rbp_exp, vsd_geneExp,correlation_method = "pearson",
                              threshold = 0, p_threshold = 0.05)


visualize_rbp_psi_correlation <- function(gene_name, subset_phe_id, psi_ratio, response_results, rbp_exp, vsd_genes_exp, all_geno_transpose_df, p_threshold = 0.05, r2_threshold = 0.2) {
  
  # Get the most significant intron junction
  mostSig_phe_id <- response_results |> 
    dplyr::filter(phe_id %in% subset_phe_id) |> 
    slice_min(FNF_p) |> 
    pull(phe_id)
  
  # Get variant info
  variantID <- response_results |> 
    dplyr::filter(phe_id %in% subset_phe_id) |> 
    slice_min(FNF_p) |> 
    pull(var_id)
  
  minorAllele <- response_results |> 
    dplyr::filter(phe_id %in% subset_phe_id) |> 
    slice_min(FNF_p) |> 
    pull(minor_allele)
  
  # Process allele information
  alleles <- unlist(strsplit(variantID, ":"))
  ref_allele <- alleles[3]
  alt_allele <- alleles[4]
  protective_allele <- ifelse(minorAllele == ref_allele, alt_allele, ref_allele)
  
  # Get PSI data
  psi_matrix <- psi_ratio[rownames(psi_ratio) %in% mostSig_phe_id, ]
  psi_data <- data.frame(sampleID = names(psi_matrix), psi = as.double(psi_matrix))
  
  # Get genotype data
  geno_data <- all_geno_transpose_df[rownames(all_geno_transpose_df) %in% variantID, ] %>%
    as.data.frame() %>%
    rownames_to_column(var = "sampleID")
  colnames(geno_data)[2] <- "genotype"
  
  # Process RBP expression data
  rbp_exp_data <- vsd_geneExp |> 
    filter(gene_id %in% rbp_exp$gene_id) |> 
    pivot_longer(cols = -c(ENSG, gene_id, SYMBOL), names_to = "sampleID", values_to = "expression")
  
  # Combine all data
  combined_data <- psi_data |> 
    left_join(rbp_exp_data, by = "sampleID") |>
    left_join(geno_data, by = "sampleID")
  
  combined_data <- combined_data %>%
    mutate(
      Group = case_when(
        grepl("CTL", sampleID) ~ "PBS",
        grepl("FNF", sampleID) ~ "Fn-f",
        TRUE ~ "Unknown"
      ),
      genotype = factor(genotype, levels = c("0", "1", "2")),
      genotype = case_when(
        genotype == "2" ~ paste(minorAllele, minorAllele, sep = "/"),
        genotype == "1" ~ paste(minorAllele, protective_allele, sep = "/"),
        genotype == "0" ~ paste(protective_allele, protective_allele, sep = "/")
      )
    ) %>%
    mutate(
      genotype = factor(genotype, levels = c(
        paste(protective_allele, protective_allele, sep = "/"),
        paste(minorAllele, protective_allele, sep = "/"),
        paste(minorAllele, minorAllele, sep = "/")
      ))
    )
  
  # Calculate R² values
  # r_squared <- combined_data %>%
  #   group_by(SYMBOL, Group) %>%
  #   do(glance(lm(psi ~ expression, data = .))) %>%
  #   mutate(label = sprintf("%s R² = %.3f", Group, r.squared))
  # Calculate R² values and p-values for each genotype
  r_squared <- combined_data %>%
    group_by(SYMBOL, genotype) %>%
    do({
      model <- lm(psi ~ expression, data = .)
      summary <- glance(model)
      data.frame(
        r.squared = summary$r.squared,
        p.value = summary$p.value
      )
    }) %>%
    mutate(label = sprintf("%s: R² = %.2f", genotype, r.squared)) %>%
    arrange(genotype) %>%
    group_by(SYMBOL) %>%
    mutate(y_position = 1 - (0.05 * (row_number() - 1)))
  
  # Filter for RBPs with at least one significant genotype and R² > threshold
  significant_rbps <- r_squared %>%
    group_by(SYMBOL) %>%
    summarise(
      any_significant = any(p.value < p_threshold),
      any_high_r2 = any(r.squared > r2_threshold),
      minor_allele_r2 = r.squared[genotype == paste(minorAllele, minorAllele, sep = "/")]
    ) %>%
    filter(any_significant & any_high_r2 & minor_allele_r2 > 0.2) %>%
    pull(SYMBOL)
  
  # Filter the data for significant RBPs
  combined_data_filtered <- combined_data %>% filter(SYMBOL %in% significant_rbps)
  r_squared_filtered <- r_squared %>% filter(SYMBOL %in% significant_rbps)
  
  # Get y-axis limits for each facet
  y_limits <- combined_data_filtered %>%
    group_by(SYMBOL) %>%
    summarise(ymin = min(psi), ymax = max(psi))
  
  # Join y_limits to r_squared
  r_squared_filtered <- r_squared_filtered %>%
    left_join(y_limits, by = "SYMBOL") %>%
    mutate(y_label = ymin + (ymax - ymin) * y_position)
  
  ggplot(combined_data_filtered, aes(x = expression, y = psi, color = genotype)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = FALSE, aes(group = genotype)) +
    facet_wrap(~ SYMBOL, scales = "free", ncol = 3) +
    scale_color_manual(values = c("#1e87a5", "#FFB81C", "#E41A1C")) +
    geom_text(data = r_squared_filtered,
              aes(x = -Inf, y = y_label, label = label, color = genotype),
              hjust = -0.1, vjust = 1, size = 3)  +
    labs(x = "RBP Expression", y = "Intron Usage (PSI)") +
    theme_minimal() +
    theme(strip.placement = "outside",
          axis.line.y =  element_line(linewidth = 0.25),
          axis.line.x =  element_line(linewidth = 0.25),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(color = "black", linewidth = 0.25),
          axis.ticks.length.y = unit(-0.1, "cm"),
          axis.title.x = element_markdown(size = 6, family= "Helvetica",
                                          margin = margin(r = -0.5)),
          axis.title.y = element_markdown(size = 6, family= "Helvetica",
                                          margin = margin(r = -0.5)),
          text = element_text(family = "Helvetica"),
          axis.text.y = element_text(color = "black", size = 6),
          axis.text.x = element_text(color = "black", size = 6, margin = margin(b = 0)),
          strip.background = element_blank(),
          panel.background = element_rect(fill = "transparent", color = "transparent"),
          plot.background = element_rect(fill = "transparent", color = "transparent"),
          panel.grid = element_blank(),
          legend.position = "none",
          panel.spacing.x = unit(1, "cm"))
}


# visualize_rbp_psi_correlation("MGRN1", mgrn1_phe_id, psi_ratio, response_fnf_results, mgrn1_rbp_exp, vsd_geneExp,all_geno_transpose_df)
# visualize_rbp_psi_correlation("OS9", os9_phe_id, psi_ratio, response_fnf_results, OS9_rbp_exp[1:6,], vsd_geneExp,all_geno_transpose_df)
#visualize_rbp_psi_correlation("SNHG29", snhg29_phe_id, psi_ratio, response_fnf_results, snhn29_rbp_exp[7:10,], vsd_geneExp,all_geno_transpose_df,
#                              p_threshold = 0.5, r2_threshold = 0.1)
#visualize_rbp_psi_correlation("CAST", os9_phe_id, all_psi_inv, response_fnf_results, cast_rbp_exp, vsd_geneExp, all_geno_transpose_df,
#                              p_threshold = 0.01, r2_threshold = 0)
#visualize_rbp_psi_correlation("CALU", calu_results$phe_id, psi_ratio, response_fnf_results, calu_results$rbp_exp, vsd_geneExp,all_geno_transpose_df,
#                              p_threshold = 0.05, r2_threshold = 0.2)


  