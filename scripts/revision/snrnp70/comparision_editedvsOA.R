#rMATs comparison between edited and the other set OA 
setwd("/work/users/s/e/seyoun/CQTL_sQTL/")
library(maser)
library(rtracklayer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)
library(colorspace)

#-------------------------------------------------------------------------------
#OA
# Comparison wd vs hetKD
path <- "output/rmats_oa_101/"

oa_pbs <- maser(path, c("oa", "pbs"), ftype = "JCEC")
geneEvents(oa_pbs, "SNRNP70")
head(summary(oa_pbs, type = "SE")[, 1:8])
oa_filt <- filterByCoverage(oa_pbs, avg_reads = 5)
sig_oa <- topEvents(oa_filt, fdr = 0.05, deltaPSI = 0.1)
#summary(sig_oa,'SE') |> dim()

#geneEvents(sig_kd, "SNRNP70")
rmats_oa_genes <- c(summary(sig_oa,'SE')[,"geneSymbol"],
                    summary(sig_oa,'RI')[,"geneSymbol"],
                    summary(sig_oa,'A3SS')[,"geneSymbol"],
                    summary(sig_oa,'A5SS')[,"geneSymbol"],
                    summary(sig_oa,'MXE')[,"geneSymbol"]) |> unique()
save(oa_pbs, sig_oa,oa_filt, rmats_oa_genes,
     file = "output/rmats_oa_101/sig_oa_pbs.RData")

#-------------------------------------------------------------------------------
#fnf
fnf_path <- "output/rmats_fnf_101_no_paired/"
fnf_pbs <- maser(fnf_path, c("fnf", "pbs"), ftype = "JCEC")
geneEvents(fnf_pbs, "SNRNP70")
head(summary(fnf_pbs, type = "SE")[, 1:8])
fnf_filt <- filterByCoverage(fnf_pbs, avg_reads = 5)
sig_fnf <- topEvents(fnf_filt, fdr = 0.05, deltaPSI = 0.1)
geneEvents(sig_fnf, "SNRNP70")
rmats_fnf_genes <- c(summary(sig_fnf,'SE')[,"geneSymbol"],
                    summary(sig_fnf,'RI')[,"geneSymbol"],
                    summary(sig_fnf,'A3SS')[,"geneSymbol"],
                    summary(sig_fnf,'A5SS')[,"geneSymbol"],
                    summary(sig_fnf,'MXE')[,"geneSymbol"]) |> unique()
save(fnf_pbs, sig_fnf, fnf_filt,rmats_fnf_genes,
     file = "output/rmats_fnf_101_no_paired/sig_fnf_pbs.RData")

#-------------------------------------------------------------------------------
#Making the psi boxplot between OA and edited

compare_and_plot <- function(test_obj, compared_obj) {
  valid_types <- c("SE", "RI", "A3SS", "A5SS", "MXE")
  plot_data_list <- NULL
  
  for (i in valid_types) {
    cat("Processing type:", i, "\n")
    
    # Get summary data
    summary_data <- summary(test_obj, type = i)
    compared_summary_data <- summary(compared_obj, type = i)
    
    cat("Rows in summary_data:", nrow(summary_data), "\n")
    cat("Rows in compared_summary_data:", nrow(compared_summary_data), "\n")
    
    # Define matching columns based on event type
    match_cols <- switch(i,
                         "SE" = c("exon_target", "exon_upstream", "exon_downstream"),
                         "RI" = c("exon_ir", "exon_upstream", "exon_downstream"),
                         "A3SS" = c("exon_long", "exon_short", "exon_flanking"),
                         "A5SS" = c("exon_long", "exon_short", "exon_flanking"),
                         "MXE" = c("exon_1", "exon_2", "exon_upstream", "exon_downstream"))
    
    cat("Matching columns:", paste(match_cols, collapse = ", "), "\n")
    
    # Ensure columns exist and are comparable
    summary_data <- summary_data %>% mutate(across(all_of(match_cols), as.character))
    compared_summary_data <- compared_summary_data %>% mutate(across(all_of(match_cols), as.character))
    
    # Merge with inner_join
    matched_data <- inner_join(
      summary_data,
      compared_summary_data,
      by = c("Chr", "Strand", match_cols),
      suffix = c("_sig", "_oa")
    )
    
    cat("Rows in matched_data:", nrow(matched_data), "\n")
    
    if (nrow(matched_data) > 0) {
      # Categorize IncLevelDifference from test_obj
      matched_data <- matched_data %>%
        mutate(IncLevelCategory = ifelse(IncLevelDifference_sig > 0, "Higher", "Lower")) %>%
        mutate(EventType = i) # Add event type for tracking
      
      plot_data_list[[i]] <- matched_data
    } else {
      cat("No matches found for type:", i, "\n")
    }
  }
  
  # Combine all data
  plot_data <- bind_rows(plot_data_list)
  return(plot_data)
}


kd_comparedtoFNF_plot_data <- compare_and_plot(sig_kd, fnf_filt)
kd_comparedtoOA_plot_data <- compare_and_plot(sig_kd, oa_filt)

#fnf_comparedtoOA_plot_data <-plot_data
#oa_comparedtoFNF_plot_data <- plot_data
#fnf_comparedtoKD_plot_data <- plot_data
#OA_comparedtoKD_plot_data <- plot_data

# To calculate the Wilcox test for the p-value


generate_plot_and_stats <- function(plot_data, y_label, colors, y_limits = c(-0.3, 0.3), highlight_gene = "SNRNP70") {
  # Filter data for higher and lower categories
  deltaPSI_high <- plot_data |> filter(IncLevelCategory == "Higher") |> dplyr::select(IncLevelDifference_oa)
  deltaPSI_low <- plot_data |> filter(IncLevelCategory == "Lower") |> dplyr::select(IncLevelDifference_oa)
  
  # Perform Wilcoxon tests
  up_test <- wilcox.test(x = as.numeric(deltaPSI_high$IncLevelDifference_oa), mu = 0, alternative = "greater")
  down_test <- wilcox.test(x = as.numeric(deltaPSI_low$IncLevelDifference_oa), mu = 0, alternative = "less")
  
  # Format p-values
  up_pvalue <- format.pval(up_test$p.value, digits = 3)
  down_pvalue <- format.pval(down_test$p.value, digits = 3)
  
# Add a highlight column for the specified gene
  plot_data <- plot_data %>%
    mutate(Highlight = ifelse(geneSymbol_sig == highlight_gene, "Highlighted", "Normal"))
  
  # Generate the plot
  plot <- ggplot(plot_data, aes(x = IncLevelCategory, y = IncLevelDifference_oa, fill = IncLevelCategory)) +
    geom_hline(yintercept = 0, lty = 2, color = "grey25", linewidth = 0.25) +
    geom_jitter(aes(color = Highlight), width = 0.2, size = 0.5) +
geom_text(
      data = plot_data %>% filter(Highlight == "Highlighted"),
      aes(label = geneSymbol_sig),
      color = "black",
      size = 2,
      vjust = -0.5
    ) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.5, alpha = 0.4, width = 0.5, color = colors) +
    stat_boxplot(geom = "errorbar", width = 0.5, color = colors) +
    scale_fill_manual(values = colors) +
scale_color_manual(values = c("Highlighted" = "red", "Normal" = "black")) +
    scale_y_continuous(name = y_label, limits = y_limits, breaks = seq(y_limits[1], y_limits[2], 0.2)) +
    coord_cartesian(clip = "off") +
    annotate("text", x = 1, y = y_limits[2] - 0.02, label = paste("p-value:", up_pvalue), size = 2) +
    annotate("text", x = 2, y = y_limits[2] - 0.02, label = paste("p-value:", down_pvalue), size = 2) +
    theme_minimal(base_family = "Helvetica") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_markdown(size = 8, color = c(darken(colors[1], 0.3), darken(colors[2], 0.3))),
      axis.title.y = element_markdown(size = 8),
      axis.text.y = element_markdown(size = 6, color = "black"),
      axis.line.y = element_line(linewidth = 0.35),
      axis.line.x = element_blank(),
      axis.ticks = element_blank(),
      strip.background = element_blank(),
      strip.text = element_markdown(size = 8),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  
  return(list(plot = plot, up_pvalue = up_pvalue, down_pvalue = down_pvalue))
}

# Example usage for kd_comparedtoOA_plot_data
kd_compared_OA_plot <- generate_plot_and_stats(
  plot_data = kd_comparedtoOA_plot_data,
  y_label = "delta PSI in response to OA/PBS",  # Updated y-axis label
  colors = c("#c8f0bf", '#e7bff0'),
  highlight_gene = "SNRNP70"  # Highlight SNRNP70
)

# Example usage for kd_comparedtoFNF_plot_data
kd_compared_fnf_plot <- generate_plot_and_stats(
  plot_data = kd_comparedtoFNF_plot_data,
  y_label = "delta PSI in response to Fn-F/PBS",  # Updated y-axis label
  colors = c("#FFB81C", '#005587'),
  highlight_gene = "SNRNP70"  # Highlight SNRNP70
)

#------------------------------------------------------------------------
# Save the plots to files
save(kd_comparedtoFNF_plot_data, kd_comparedtoOA_plot_data, file = "output/rmats_edited/sig_kd_comparison.RData")
save(kd_compared_OA_plot, file = "output/results_plots/revision_plots/kd_comparedtoOA_boxplot.rda")
save(kd_compared_fnf_plot, file = "output/results_plots/revision_plots/kd_comparedtoFNF_boxplot.rda")

#------------------------------------------------------------------------
library(plotgardener)

pdf(file = "output/results_plots/revision_plots/rmats_splicing_comparison.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 3.5)

pageCreate(width = 6, height =3.5 , default.units = "inches", showGuides = FALSE)
load("output/results_plots/revision_plots/kd_comparedtoOA_boxplot.rda") #kd_compared_OA_plot
load("output/results_plots/revision_plots/kd_comparedtoFNF_boxplot.rda") #kd_compared_fnf_plot

plotGG(kd_compared_fnf_plot$plot , x = 0, y = 0.5,
       width = 2.7, height =2.5)

plotGG(kd_compared_OA_plot$plot , x = 3, y = 0.5,
       width = 2.7, height =2.5)

plotSegments(
    x0 = 0.65, y0 = 3, x1 = 2.3, y1 = 3,
    default.units = "inches",
    lwd = 0.5, lty = 1
)

plotSegments(
    x0 = 3.65, y0 = 3, x1 = 3.65+1.7, y1 = 3,
    default.units = "inches",
    lwd = 0.5, lty = 1
)

plotText("Splice exon usage in edited SNRNP70 alt ex8", x = 0.65, y = 3.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 6)
plotText("Splice exon usage in edited SNRNP70 alt ex8", x = 3.65, y = 3.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 6)
dev.off()
