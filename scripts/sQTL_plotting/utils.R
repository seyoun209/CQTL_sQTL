# Utils for the sqtl 
## Author: Seyoun Byun
## Date: 04.22.2024
## Edited:
#-------------------------------------------------------------------------------
library(rtracklayer)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(broom)
library(stringr)
library(tidyverse)
setwd("/work/users/s/e/seyoun/CQTL_sQTL")

venn_font <- function(p, font){
  
  grep_grob <- function(gt, lab){
    which(sapply(gt, function(x) grepl(lab, x$name)))
  }
  
  p2 <- ggplot_gtable(ggplot_build(p))
  mygrobs <- p2$grobs
  # Break down grobs into panel, venn object, and text pieces
  panel_grob <- mygrobs[[grep_grob(mygrobs, "panel")]]
  venn_grob <- panel_grob$children[[grep_grob(panel_grob$children, "venn")]]
  text_grobs <- venn_grob$children[grep_grob(venn_grob$children, "text")]
  # Make both new font family
  text_grobs <- do.call(grid::gList, 
                        lapply(text_grobs, 
                               function(x) {x$gp$fontfamily <- font; 
                               x}))
  # Make titles bold
  text_grobs[[1]]$gp$fontface <- "bold"
  
  # Add grobs back
  venn_grob$children[grep_grob(venn_grob$children, "text")] <- text_grobs
  panel_grob$children[[grep_grob(panel_grob$children, "venn")]] <- venn_grob
  mygrobs[[grep_grob(mygrobs, "panel")]] <- panel_grob
  p2$grobs <- mygrobs
  grid::grid.newpage()
  grid::grid.draw(p2)
  
}


hg38_intron <- fread("/work/users/s/e/seyoun/CQTL_sQTL/tools/leafcutter/gencode_hg38_all_introns.bed.gz")
colnames(hg38_intron) <- c("chr","start","end","symbol","ensg","strand","enst","id","region","info")
#hg38_intron_filtered <- hg38_intron %>%
#  filter(!startsWith(symbol, "ENSG"))

hg38_intron_sub <- hg38_intron %>%
  mutate(genomicLoc = paste(chr, start, end, sep = ":")) %>%
  dplyr::select(genomicLoc, ensg)

hg38_intron_sub_select_first <- hg38_intron_sub %>%
  group_by(genomicLoc) %>%
  arrange(ensg) %>%  # Ensure that the data is sorted by ensg within each genomicLoc
  dplyr::slice(1) %>%       # Select the first ensg for each genomicLoc
  ungroup()


leafcutter_pheno <- fread("output/gtex_cluster/ctl_fnf.leafcutter.phenotype_groups.txt",header=F)
colnames(leafcutter_pheno) <- c("id_ensg","ensg")
leafcutter_pheno_mutated <- leafcutter_pheno %>%
  mutate(
    chr=gsub("\\..*$", "", sapply(strsplit(as.character(id_ensg), ":"), `[`, 1)),
    start=gsub("\\..*$", "", sapply(strsplit(as.character(id_ensg), ":"), `[`, 2)),
    end=gsub("\\..*$", "", sapply(strsplit(as.character(id_ensg), ":"), `[`, 3)),
    clusterID=gsub("\\..*$", "", sapply(strsplit(as.character(id_ensg), ":"), `[`, 4)),
    phe_id=sapply(strsplit(as.character(id_ensg), ":"), function(x) paste(x[1:4], collapse = ":")),
    genomicLoc=sapply(strsplit(as.character(id_ensg), ":"), function(x) paste(x[1:3], collapse = ":"))
  )
leafcutter_pheno_subset <- leafcutter_pheno_mutated %>%
  dplyr::select(ensg, clusterID) %>%
  group_by(clusterID) %>%
  arrange(ensg) %>%  # Ensure that the data is sorted by ensg within each genomicLoc
  dplyr::slice(1) %>%       # Select the first ensg for each genomicLoc
  ungroup()




gtf_path <- "/users/s/e/seyoun/Ref/genome/gencode.v45.annotation.gtf"
gtf_data <- import(gtf_path)
gene_info <- gtf_data[gtf_data$type == "gene"]

genes_df <- data.frame(
  gene_id = mcols(gtf_data)$gene_id,
  gene_name = mcols(gtf_data)$gene_name,
  stringsAsFactors = FALSE
)

genes_df_unique <- genes_df %>%
  dplyr::rename(ensg = gene_id) %>%
  dplyr::distinct(ensg, gene_name, .keep_all = TRUE)  


txdb <- loadDb("/work/users/s/e/seyoun/crispr/02.test_seq/gencode.v45.annotation.TxDb")
txdb_genes <- genes(txdb)


annotate_genes <- function(cond_pbs, hg38_intron_sub_select_first, leafcutter_pheno_subset, txdb_genes) {
  # Processing initial data
  conditional_pbs <- cond_pbs %>%
    mutate(
      clusterID = gsub("\\..*$", "", sapply(strsplit(as.character(phe_id), ":"), `[`, 4)),
      genomicLoc = sapply(strsplit(as.character(phe_id), ":"), function(x) paste(x[1:3], collapse = ":"))
    )
  
  # First annotation trial
  conditional_sGene_pbs_anno_v1 <- conditional_pbs %>%
    left_join(hg38_intron_sub_select_first, by = "genomicLoc") 
  conditional_sGene_pbs_NO_anno_v1 <- conditional_sGene_pbs_anno_v1[is.na(conditional_sGene_pbs_anno_v1$ensg),]
  
  # Second annotation trial
  conditional_sGene_pbs_anno_v2 <- conditional_sGene_pbs_NO_anno_v1 %>%
    left_join(leafcutter_pheno_subset, by = "clusterID") %>%
    mutate(ensg = coalesce(ensg.x, ensg.y)) %>%
    dplyr::select(-ensg.x, -ensg.y)
  conditional_sGene_pbs_NO_anno_v2 <- conditional_sGene_pbs_anno_v2[is.na(conditional_sGene_pbs_anno_v2$ensg),]
  
  # Third annotation trial with genomic ranges
  conditional_sGene_pbs_gr <- GRanges(
    seqnames = conditional_sGene_pbs_NO_anno_v2$phe_chr,
    ranges = IRanges(start = conditional_sGene_pbs_NO_anno_v2$phe_from-1, end = conditional_sGene_pbs_NO_anno_v2$phe_to)
  )
  overlaps <- mergeByOverlaps(conditional_sGene_pbs_gr, txdb_genes, type="within") 
  overlaps_df <- data.frame(
    gene_id = overlaps$gene_id,
    genomicLoc = paste0(seqnames(overlaps$conditional_sGene_pbs_gr), ":", start(overlaps$conditional_sGene_pbs_gr), ":", end(overlaps$conditional_sGene_pbs_gr))
  )
  overlaps_df_selected_first <- overlaps_df %>% 
    group_by(genomicLoc) %>%
    arrange(gene_id) %>%
    dplyr::slice(1) %>%
    ungroup()
  
  conditional_sGene_pbs_anno_v3 <- conditional_sGene_pbs_NO_anno_v2 %>%
    left_join(overlaps_df_selected_first, by = "genomicLoc") %>%
    mutate(ensg = coalesce(gene_id, ensg)) %>%
    dplyr::select(-gene_id) 
  
  # Merging all
  conditional_pbs_final <- rbind(
    conditional_sGene_pbs_anno_v1[!is.na(conditional_sGene_pbs_anno_v1$ensg),],
    conditional_sGene_pbs_anno_v2[!is.na(conditional_sGene_pbs_anno_v2$ensg),], 
    conditional_sGene_pbs_anno_v3
  )
  
  return(conditional_pbs_final)
}


#making a pie chart for the response QTL
create_pie_chart <- function(data, colors, group) {
  # Compute percentages and cumulative positions
  data$Percentage <- round((data$Freq / sum(data$Freq)), 3) * 100
  data$ymax <- cumsum(data$Freq / sum(data$Freq))
  data$ymin <- c(0, head(data$ymax, n = -1))
  data$labelPosition <- (data$ymax + data$ymin) / 2
  
  # Create the pie chart
  p <- ggplot(data, aes(x = "", y = Percentage, fill = group)) +
    geom_bar(stat = "identity") +
    theme_void() +
    scale_fill_manual(values = colors) +
    theme(legend.position = "none", plot.background = element_rect(fill = 'transparent', color = NA))
  
  if (group == "PBS") {
    p <- p + coord_polar("y", direction = -1, start = 0)
  } else if (group == "FNF") {
    p <- p + coord_polar("y", direction = 1, start = 0)
  }
  
  return(p)
}



#Dot from for the line plot 
get_y_for_median <- function(data, median_x, binwidth) {
  bins <- hist(data$abs_distance, breaks = seq(0, max(data$abs_distance) + binwidth, by = binwidth), plot = FALSE)
  bin_index <- which(bins$mids >= median_x)[1]
  return(bins$counts[bin_index])
}



# Change the ranks if it's exist two times compared to sQTL data
adjust_ranks_for_zero <- function(data) {
  
  # Identify the id groups that have a rank 0
  ids_with_rank_0 <- data %>%
    group_by(phe_id) %>%
    dplyr::filter(any(rank == 0)) %>%
    ungroup() %>%
    distinct(phe_id) %>%
    dplyr::pull(phe_id)
  
  # Adjust ranks only for these id groups
  adjusted_data <- data %>%
    group_by(phe_id) %>%
    dplyr::mutate(
      new_rank = if_else(phe_id %in% ids_with_rank_0,
                         rank + 1,
                        rank)
    ) %>%
    ungroup()
  
  # Replace the original rank column
  adjusted_data <- adjusted_data %>% 
    dplyr::select(-rank) %>% 
    dplyr::rename("rank" = "new_rank")
  
  return(adjusted_data)
}

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




#make boxplot for the selected gene 

create_boxplot <- function(gene_name, highConf_resQtL) {

  # Filter data for the specified gene
  test_boxplotInfo <- highConf_resQtL |> dplyr::filter(SYMBOL %in% gene_name)
  if (nrow(test_boxplotInfo) != 1) {
    test_boxplotInfo <- test_boxplotInfo |> slice_min(interaction_pval)
  }
  minorAllele <- unique(test_boxplotInfo$minor_allele)
  
  # Process PSI data
  psi_matrix <- all_psi_inv[rownames(all_psi_inv) %in% test_boxplotInfo$phe_id, ]
  psi_data <- data.frame(sampleID = names(psi_matrix), psi = as.double(psi_matrix))
  
  # Process variant data
  variantID <- unique(test_boxplotInfo$var_id)
  rsID <- unique(test_boxplotInfo$rsID)
  alleles <- unlist(strsplit(variantID, ":"))
  ref_allele <- alleles[3]
  alt_allele <- alleles[4]
  protective_allele <- ifelse(minorAllele == ref_allele, alt_allele, ref_allele)
  
  # Process genotype data
  geno_data <- all_geno_transpose_df[rownames(all_geno_transpose_df) %in% variantID, ] %>%
    as.data.frame() %>%
    rownames_to_column(var = "sampleID")
  colnames(geno_data)[2] <- "genotype"
  
  # Combine metadata
  meta_data_fi <- geno_data %>%
    left_join(meta_catl_simple, by = c("sampleID" = "ID")) %>%
    left_join(psi_data, by = c("sampleID" = "sampleID"))
  
  # Prepare final dataset
  meta_combined_all <- meta_data_fi %>%
    mutate(
      genotype = factor(genotype, levels = c("0", "1", "2")),
      Condition = ifelse(Condition == "CTL", "PBS", ifelse(Condition == "FNF", "FN-f", NA)),
      Condition = factor(Condition, levels = c("PBS", "FN-f")),
      genotype = case_when(
        genotype == "2" ~ paste(minorAllele, minorAllele, sep = "/"),
        genotype == "1" ~ paste(protective_allele, minorAllele, sep = "/"),
        genotype == "0" ~ paste(protective_allele, protective_allele, sep = "/")
      )
    ) %>%
    mutate(
      genotype = factor(genotype, levels = c(
        paste(protective_allele, protective_allele, sep = "/"),
        paste(protective_allele, minorAllele, sep = "/"),
        paste(minorAllele, minorAllele, sep = "/")
      ))
    )
  
  # Calculate y-axis range
  maxrange <- range(meta_combined_all$psi, na.rm = TRUE)
  yaxis_range_min <- min(0, round(maxrange[1], 1))
  yaxis_range_max <- round(maxrange[2], 1)
  
  if (yaxis_range_min == 0 & yaxis_range_max == 0) {
    yaxis_range_min <- maxrange[1]
    yaxis_range_max <- maxrange[2]
  }
  
  # Prepare text data
  text_data <- data.frame(
    Condition = c("PBS", "FN-f"),
    beta = c(test_boxplotInfo$PBS_beta[1], test_boxplotInfo$FNF_beta[1]),
    pvalue = c(test_boxplotInfo$PBS_p[1], test_boxplotInfo$FNF_p[1])
  ) |> mutate(Condition = factor(Condition, levels = c("PBS", "FN-f")))
  
  # Create plot
  p <- ggplot(meta_combined_all, aes(x = genotype, y = psi, fill = Condition)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.25, alpha = 0.7) +
    scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +  
    geom_point(color = "grey40", position = position_jitterdodge(), size = 0.25) +
    labs(x = "Genotype", y = "Intron usage (PSI)") +
    scale_fill_manual(values = c('#0067B9', '#FCCE52')) +
    scale_color_manual(values = c('#0067B9', '#FCCE52')) +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(strip.placement = "outside",
          axis.line.y = element_line(linewidth = 0.25),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(color = "black", linewidth = 0.25),
          axis.ticks.length.y = unit(-0.1, "cm"),
          axis.title.x = element_blank(),
          axis.title.y = element_markdown(size = 6, family = "Helvetica",
                                          margin = margin(r = -0.5)),
          text = element_text(family = "Helvetica"),
          axis.text.y = element_text(color = "black", size = 6),
          axis.text.x = element_text(color = "black", size = 6, margin = margin(b = 0)),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.background = element_rect(fill = "transparent", color = "transparent"),
          plot.background = element_rect(fill = "transparent", color = "transparent"),
          panel.grid = element_blank(),
          legend.position = "none",
          panel.spacing.x = unit(0.5, "cm")) +
    facet_wrap(~ Condition, ncol = 2) +
    geom_text(data = text_data,
              aes(x = 1.5, y = Inf ,
                  label = sprintf("%s\n Beta: %.3f\n p-value: %.2e",
                                  Condition, beta, pvalue),
                  color = Condition),
              hjust = 0.5, vjust = 1, size = 2, family = "Helvetica",lineheight = 0.75,
              check_overlap = TRUE)
  
  return(p)
}

pbs_norm_qtl_pc5_list <- readRDS("output/nominal_1mb/nominal_pbs/pbs_norm_qtl_pc5_list.rds")
fnf_norm_qtl_pc4_list <- readRDS("output/nominal_1mb/nominal_fnf/fnf_norm_qtl_pc4_list.rds")

load("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/freq/maf_id_add.rds") #maf_id_add  is the name
maf_subset_rsID <- maf_id_add |> dplyr::select(-c("minor_allele","MAF"))

#Locus zoom plot gene level



plot_sQTL_manhattan <- function(gene_name, primary_dataset, x_start, y_start, width, height, zoom_range = 200000) {
  # Validate input
  if (!primary_dataset %in% c("PBS", "FNF")) {
    stop("Primary dataset must be either 'PBS' or 'FNF'")
  }
  
  # Select the appropriate dataset
  highConf_resQtL <- if(primary_dataset == "PBS") response_pbs_results else response_fnf_results
  
  # Filter data for the specified gene
  test_boxplotInfo <- highConf_resQtL |> dplyr::filter(SYMBOL %in% gene_name) 
  if (nrow(test_boxplotInfo) != 1) {
    test_boxplotInfo <- test_boxplotInfo |> slice_min(!!sym(paste0(primary_dataset, "_p")))
  }
  
  # Extract necessary information
  chrom <- test_boxplotInfo$phe_chr
  variantID <- unique(test_boxplotInfo$var_id)
  rsID <- test_boxplotInfo$rsID
  intronID <- test_boxplotInfo$phe_id
  
  # Load LD data
  ld_file <- paste0("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld/", chrom, "/", variantID, ".ld")
  ld_calc <- fread(ld_file) |> dplyr::select(c("SNP_B","R2")) |>
    dplyr::rename("var_id" = "SNP_B")
  
  # Load and process QTL data for both PBS and FNF
  pbs_norm_qtl <- pbs_norm_qtl_pc5_list[[chrom]]
  fnf_norm_qtl <- fnf_norm_qtl_pc4_list[[chrom]]
  
  pbs_qtl_region <- pbs_norm_qtl |> dplyr::filter(phe_id %in% intronID)
  fnf_qtl_region <- fnf_norm_qtl |> dplyr::filter(phe_id %in% intronID)
  
  pbs_qtl_pval_rsid <- left_join(pbs_qtl_region, maf_subset_rsID, by="var_id")
  fnf_qtl_pval_rsid <- left_join(fnf_qtl_region, maf_subset_rsID, by="var_id")
  
  leftjoin_pbs <- inner_join(pbs_qtl_pval_rsid, ld_calc, by = c("var_id" = "var_id"))
  leftjoin_fnf <- inner_join(fnf_qtl_pval_rsid, ld_calc, by = c("var_id" = "var_id"))
  
  prepare_locus_plot <- function(data) {
    data |> 
      mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1))) |>
      mutate(LDgrp = addNA(LDgrp)) |>
      dplyr::rename(chrom = "var_chr",
                    pos = "var_from",
                    p = "nom_pval",
                    LD = "R2",
                    snp = "rsID") |>
      dplyr::select("chrom", "pos", "p", "snp", "LD", "LDgrp", "phe_id", "var_id") |> 
      mutate(LDgrp = factor(LDgrp, levels = c(NA,"(0.8,1]","(0.6,0.8]", "(0.4,0.6]","(0.2,0.4]","(0,0.2]"), ordered = TRUE)) |>
      filter(!is.na(LDgrp)) |>
      arrange(desc(LDgrp)) |>  # Sort by LDgrp so that higher LD values are at the end
      data.frame()
  }
  
  # Prepare plotting data for both PBS and FNF
  pbs_locus_plot <- prepare_locus_plot(leftjoin_pbs)
  fnf_locus_plot <- prepare_locus_plot(leftjoin_fnf)
  
  # Set plot region with adjustable zoom
  minregion <- test_boxplotInfo$phe_from - zoom_range
  maxregion <- test_boxplotInfo$phe_to + zoom_range
  
  # Set up plot parameters
  region_pg <- pgParams(assembly = "hg38", chrom = chrom,
                        chromstart = minregion,
                        chromend = maxregion,
                        x = x_start, y = y_start, width = width, height = height)
  
  # Calculate y-axis limits
  pbs_ylim <- ceiling(max(log10(pbs_locus_plot$p)*-1)) + 1
  fnf_ylim <- ceiling(max(log10(fnf_locus_plot$p)*-1)) + 1
  ylim_pg <- max(pbs_ylim, fnf_ylim)
  
  # Function to create Manhattan plot
  create_manhattan_plot <- function(data, y_offset, label) {
    plot_height <- (height - 0.2) / 2  # Reduce height to make room for space between plots
    locus_plot <- plotManhattan(data = data,
                                params = region_pg,
                                range = c(0, ylim_pg),
                                fill = colorby("LDgrp",
                                               palette = colorRampPalette(c("#DD3931",  
                                                                            "#EEA741", 
                                                                            "#499A53",  
                                                                            "#98CDED",  
                                                                            "#262C74" 
                                               ))),
                                y = y_start + y_offset, height = plot_height,
                                snpHighlights = data.frame(snp = rsID,
                                                           pch = c(24),
                                                           cex = c(0.5),
                                                           col = c("black")))
    
    # Add y-axis
    annoYaxis(plot = locus_plot, at = seq(0, ylim_pg, 2),
              axisLine = TRUE, fontsize = 6)
    
    # Add label
    plotText(label = paste(label, "sQTL"), 
             x = x_start + 0.1, y = y_start + y_offset, 
             just = c("left", "top"),
             fontfamily = "Helvetica", fontsize = 8)
    
    return(locus_plot)
  }
  
  # Create both PBS and FNF plots with space between them
  pbs_plot <- create_manhattan_plot(pbs_locus_plot, 0, "PBS")
  fnf_plot <- create_manhattan_plot(fnf_locus_plot, (height - 0.2) / 2 + 0.2, "FN-f")
  
  # Add y-axis label
  plotText(label = "-log10(p-value)", 
           x = x_start - 0.3, y = y_start + height/2, 
           rot = 90, fontsize = 6, just = "center",
           default.units = "inches", fontfamily = "Helvetica")
  
  # Add legend (adjusted position)
  legend_x <- x_start + width - 0.5
  legend_y <- y_start + height - 2  # Moved up to accommodate space between plots
  plotLegend(legend = c("0.8 - 1.0", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "0.0 - 0.2"),
             fill = c("#DD3931", "#EEA741", "#499A53", "#98CDED", "#262C74"),
             x = legend_x, y = legend_y, width = 0.1, height = 0.35, border = FALSE,
             fontsize = 6)
  
  # Add rsID label
  plotText(label = rsID, 
           x = x_start + width/2, y = y_start - 0.1, 
           just = c("center", "top"),
           fontfamily = "Helvetica", fontsize = 8)
  
  # Add gene track (adjusted position)
  gene_hl <- test_boxplotInfo$SYMBOL
  plotgenes <- plotGenes(params = region_pg, 
                         y = y_start + height + 0.1,
                         height = 0.4,
                         geneHighlights = data.frame("gene" = gene_hl,
                                                     "color" = "#37a7db"), fontsize = 6,geneOrder=gene_hl)
  annoGenomeLabel(plot = plotgenes, params = region_pg, fontsize = 6, 
                  y = y_start + height + 0.55)
}


