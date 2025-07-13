#This is the locus zoom plot for the 
setwd("/work/users/s/e/seyoun/CQTL_sQTL")
library(plotgardener)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(data.table)
library(grid)
library(ggpubr)
library(ggforce)

response_pbs_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_pbs_resInv.rds")
response_fnf_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_fnf_resInv.rds")

pbs_norm_qtl_pc5_list <- readRDS("output/nominal_1mb/nominal_pbs/pbs_norm_qtl_pc5_list.rds")
fnf_norm_qtl_pc4_list <- readRDS("output/nominal_1mb/nominal_fnf/fnf_norm_qtl_pc4_list.rds")

load("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/freq/maf_id_add.rds") #maf_id_add  is the name
maf_subset_rsID <- maf_id_add |> dplyr::select(-c("minor_allele","MAF"))

pbs_signals_not_Prime <- response_pbs_results |> dplyr::filter(rank != 0) |>
  arrange(PBS_p)
fnf_signals_not_Prime <- response_fnf_results |> dplyr::filter(rank != 0) |>
  arrange(FNF_p)
# Example for CFL1 -PBS chr11:65856242:65857093:clu_6119_- or SRP14 chr7:38178420:38201714:clu_33847_+ HLA-C chr6:31269525:31355107:clu_30858_- 



# Example for GBP3 chr5:180301608:180304772:clu_29793_-
fnf_signals_not_Prime |> dplyr::filter(!SYMBOL %in% pbs_signals_not_Prime$SYMBOL ,minor_alle_count >5, abs(dist_phe_var) > 5000)

plot_condtional_locuszoom <- function(phe_id_test, primary_dataset, x_start, y_start, width, height, zoom_range = 200000) {
  # Validate input
  if (!primary_dataset %in% c("PBS", "FNF")) {
    stop("Primary dataset must be either 'PBS' or 'FNF'")
  }
  
  # Select the appropriate dataset
  highConf_resQtL <- if(primary_dataset == "PBS") response_pbs_results else response_fnf_results
  
  # Filter data for the specified gene(s)
  test_boxplotInfo <- highConf_resQtL |> dplyr::filter(phe_id %in% phe_id_test) 
  
  # Extract necessary information
  chrom <- test_boxplotInfo$phe_chr
  variantID <- test_boxplotInfo$var_id
  rsID <- test_boxplotInfo$rsID
  intronID <- test_boxplotInfo$phe_id
  ranks <- test_boxplotInfo$rank
  
  locus_plot_datasets <- NULL
  signal_plot_datasets <- list()
  
  # Prepare rank 0 data from highConf_resQtL
  rank0_data <- test_boxplotInfo |> dplyr::filter(rank == 0)
  
  # Load LD data for rank 0
  ld_file <- paste0("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld/", chrom[1], "/", variantID[1], ".ld")
  ld_calc <- fread(ld_file) |> 
    dplyr::select(c("SNP_B", "R2")) |>
    dplyr::rename("var_id" = "SNP_B")
  
  # Prepare plotting data for rank 0
  norm_qtl_dataset <- if(primary_dataset == "PBS") pbs_norm_qtl_pc5_list[[chrom[1]]] else fnf_norm_qtl_pc4_list[[chrom[1]]]
  qtl_region <- norm_qtl_dataset |> dplyr::filter(phe_id == intronID[1])
  qtl_pval_rsid <- left_join(qtl_region, maf_subset_rsID, by = "var_id")
  leftjoin_qtl_data <- inner_join(qtl_pval_rsid, ld_calc, by = "var_id")
  
  prepare_locus_plot <- function(data) {
    data |> 
      mutate(
        LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)),
        LDgrp = addNA(LDgrp),
        LDgrp = factor(LDgrp, levels = c(NA, "(0.8,1]", "(0.6,0.8]", "(0.4,0.6]", "(0.2,0.4]", "(0,0.2]"), ordered = TRUE)
      ) |>
      dplyr::rename(
        chrom = "var_chr",
        pos = "var_from",
        p = "nom_pval",
        LD = "R2",
        snp = "rsID"
      ) |>
      dplyr::select("chrom", "pos", "p", "snp", "LD", "LDgrp", "phe_id", "var_id") |> 
      filter(!is.na(LDgrp)) |>
      arrange(desc(LDgrp)) |>
      as.data.frame()
  }
  
  locus_plot_datasets <- prepare_locus_plot(leftjoin_qtl_data)
  
  # Prepare data for conditional analysis (rank 1 and above)
  cond_dir <- if (primary_dataset == "PBS") "/work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/cond_norm_pbs/" else "/work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/cond_norm_fnf/"
  
  for (i in seq_along(ranks)) {  
    signal_file <- paste0(cond_dir, intronID[i], "_rank", ranks[i])
    signal_data <- fread(signal_file)
    colnames(signal_data) <- c("phe_id", "phe_chr", "phe_from", "phe_to", "phe_strd",
                               "n_var_in_cis", "dist_phe_var", "var_id", "var_chr", "var_from",
                               "var_to", "nom_pval", "r_squared", "slope", "slope_se", "best_hit")
    best_hit <- signal_data |> dplyr::filter(best_hit == 1) |> dplyr::pull(var_id)
    
    ld_file <- paste0("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld/", chrom[i], "/", variantID[i], ".ld")
    ld_calc <- fread(ld_file) |> 
      dplyr::select(c("SNP_B", "R2")) |>
      dplyr::rename("var_id" = "SNP_B")
    
    leftjoin_signal_data <- signal_data |> 
      dplyr::select("phe_id", "var_id", "var_chr", "var_from", "nom_pval", "slope") |> 
      dplyr::filter(phe_id == intronID[i]) |>
      left_join(maf_subset_rsID, by = "var_id") |>
      inner_join(ld_calc, by = "var_id")
    
    signal_plot_datasets[[i]] <- prepare_locus_plot(leftjoin_signal_data)
  }
  
  # Set plot region with adjustable zoom
  minregion <- min(test_boxplotInfo$phe_from) - zoom_range
  maxregion <- max(test_boxplotInfo$phe_to) + zoom_range
  
  # Set up plot parameters
  region_pg <- pgParams(
    assembly = "hg38",
    chrom = chrom[1],
    chromstart = minregion,
    chromend = maxregion,
    x = x_start,
    y = y_start,
    width = width,
    height = height
  )
  
  # Calculate y-axis limits
  ylim_pg <- ceiling(max(-log10(locus_plot_datasets$p))+1)
  
  total_plot_height <- height - 1.2  # Reserve more space for legend and gene track
  
  # Function to create Manhattan plot
  # Function to create Manhattan plot
  create_manhattan_plot <- function(data, y_offset, label, rsid=NULL) {
    plot_height <- (total_plot_height - 0.1 * (length(signal_plot_datasets))) / (1 + length(signal_plot_datasets))
    
    # Add label to the right of the plot
    plotText(
      label = paste(label,"-", rsid), 
      x = x_start ,
      y = y_start + y_offset, 
      just = c("left", "bottom"),
      fontfamily = "Helvetica",
      fontsize = 8
    )
    
    locus_plot <- plotManhattan(
      data = data,
      params = region_pg,
      range = c(0, ylim_pg),
      fill = colorby("LDgrp", palette = colorRampPalette(c("#DD3931", "#EEA741", "#499A53", "#98CDED", "#262C74"))),
      y = y_start + y_offset,
      height = plot_height,
      snpHighlights = data.frame(
        snp = rsID,
        pch = rep(23, length(rsID)),
        cex = rep(0.5, length(rsID)),
        col = rep("#DD3931", length(rsID))
      )
    )
    
    # Add y-axis
    annoYaxis(plot = locus_plot, at = seq(0, ylim_pg, 2), axisLine = TRUE, fontsize = 6)
    
    return(locus_plot)
  }
  
  # Calculate plot heights and offsets
  plot_height <- (total_plot_height - 0.3 * length(signal_plot_datasets)) / (1 + length(signal_plot_datasets))
  y_offsets <- seq(from = 0.7, by = plot_height + 0.3, length.out = 1 + length(signal_plot_datasets))
  
  # Create plots for conditional analysis results (in reverse order)
  signal_labels <- c("primary sQTL signal", "secondary sQTL signal", "tertiary sQTL signal", "quaternary sQTL signal")
  for (i in rev(seq_along(signal_plot_datasets))) {
    label <- ifelse(i <= length(signal_labels), 
                    signal_labels[i], 
                    paste("signal", i))
    create_manhattan_plot(signal_plot_datasets[[i]], y_offsets[i+1], label, rsID[i])
  }
  
  # Create plot for rank 0 (locus plot dataset)
  create_manhattan_plot(locus_plot_datasets, y_offsets[1], "sQTL signal before conditional analysis")
  
  legend_x <- x_start + width-0.5 
  legend_y <- y_start+0.5   # Moved up to accommodate space between plots
           
  # Add gene track (adjusted position)
  gene_hl <- test_boxplotInfo$SYMBOL
  plotgenes <- plotGenes(params = region_pg, 
                         y = y_start + height-0.25,
                         height = 0.4,
                         geneHighlights = data.frame("gene" = gene_hl,
                                                     "color" = "#37a7db"), fontsize = 6, geneOrder=gene_hl)
  annoGenomeLabel(plot = plotgenes, params = region_pg, fontsize = 6, 
                  y = y_start + height+0.2 )
  
  
  plotText(
    label = paste("Conditional sQTL mapping-" ,primary_dataset), 
    x = x_start+(width/2) ,
    y = y_start+0.35 , 
    just = "center",
    fontfamily = "Helvetica",
    fontsize = 10
  )
  
  # Add legend (adjusted position)

  plotLegend( legend = c("0.8 - 1.0", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "0.0 - 0.2"),
              fill = c("#DD3931", "#EEA741", "#499A53", "#98CDED", "#262C74"),
              x = legend_x, y = legend_y, width = 0.1, height = 0.35, border = FALSE,
              fontsize = 6)
}

# Example usage
pdf(file = "output/results_plots/Supplementary_figures/supp_fig5_example_conditional.pdf",   # The directory you want to save the file in
    width = 8.5, # The width of the plot in inches
    height = 4.5)
pageCreate(width = 8.5, height = 4.5, default.units = "inches", showGuides = FALSE)

plotText("a", x = 0.1, y = 0.1, just = c("center", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

# PBS example for great!
plot_condtional_locuszoom(
  phe_id_test = "chr5:96768431:96776404:clu_29391_-",
  primary_dataset = "PBS",
  x_start = 0.7,
  y_start = 0.25,
  width = 3.2,
  height = 3.5,  # Increased height to accommodate all elements
  zoom_range = 100000
)


plotText("b", x = 4.4, y = 0.1, just = c("center", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

plot_condtional_locuszoom(
  phe_id_test = "chr6:31410797:31411209:clu_31910_+",
  primary_dataset = "FNF",
  x_start = 4.9,
  y_start = 0.25,
  width = 3.2,
  height = 3.5,  # Increased height to accommodate all elements
  zoom_range = 100000
)

dev.off()



