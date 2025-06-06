setwd("/work/users/s/e/seyoun/CQTL_sQTL/")
library(maser)
library(rtracklayer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)
library(colorspace)
library(ggtext)
library(ggrepel)
library(ggpmisc)
library(httpgd)
library(plotgardener)


load("output/rmats_edited/sig_kd_wd.RData") #kd_wd, sig_kd, kd_filt, rmats_kd_genes
load("output/rmats_oa_101/sig_oa_pbs.RData") #oa_pbs, sig_oa,oa_filt, rmats_oa_genes
#-------------------------------------------------------------------------------

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
      suffix = c("_KD", "_fnf")
    )
    
    cat("Rows in matched_data:", nrow(matched_data), "\n")
    
    if (nrow(matched_data) > 0) {
      # Categorize IncLevelDifference from test_obj
      matched_data <- matched_data %>%
        mutate(IncLevelCategory = ifelse(IncLevelDifference_KD > 0, "Higher", "Lower")) %>%
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
kd_all_fnf_plot_data <- compare_and_plot(kd_filt, fnf_filt)

#-------------------------------------------------------------------------------
#Making function for the scatter plot 

generate_scatter_plot_and_stats <- function(plot_data, y_label, colors, y_limits = c(-0.3, 0.3), highlight_gene = "SNRNP70") {
  # Filter data for higher and lower categories
  deltaPSI_high <- plot_data |> filter(IncLevelCategory == "Higher") |> dplyr::select(IncLevelDifference_fnf)
  deltaPSI_low <- plot_data |> filter(IncLevelCategory == "Lower") |> dplyr::select(IncLevelDifference_fnf)
  
  # Perform Wilcoxon tests
  up_test <- wilcox.test(x = as.numeric(deltaPSI_high$IncLevelDifference_fnf), mu = 0, alternative = "greater")
  down_test <- wilcox.test(x = as.numeric(deltaPSI_low$IncLevelDifference_fnf), mu = 0, alternative = "less")
  
  # Format p-values
  up_pvalue <- format.pval(up_test$p.value, digits = 3)
  down_pvalue <- format.pval(down_test$p.value, digits = 3)
  
# Add a highlight column for the specified gene
  plot_data <- plot_data %>%
    mutate(Highlight = ifelse(geneSymbol_KD == highlight_gene, "Highlighted", "Normal"))
  
generate_scatter_plot_and_stats <- function(plot_data, x_label, y_label, highlight_gene = "SNRNP70") {
  # Perform Pearson correlation test (if desired; stat_poly_eq will fit its own lm)
  cor_test <- cor.test(plot_data$IncLevelDifference_oa, plot_data$IncLevelDifference_sig)
  cor_r <- round(cor_test$estimate, 3)
  cor_p <- format.pval(cor_test$p.value, digits = 3)
  
  # Create a highlight column: mark the specified gene
  plot_data <- plot_data %>%
    mutate(Highlight = ifelse(geneSymbol_KD == highlight_gene, "Highlighted", "Normal"))

    x_label <- "delta PSI Het-KD/WT"
  y_label <- "delta PSI FNF/PBS"
  # Generate the scatter plot
  scatter_plot <- ggplot(plot_data, aes(x = IncLevelDifference_fnf, y = IncLevelDifference_KD)) +
    geom_point(aes(color = Highlight), size = 1) +
    scale_color_manual(values = c("Highlighted" = "red", "Normal" = "black")) +
    labs(x = x_label,
         y = y_label) +
    stat_poly_eq(
      use_label(c("R2", "P")),
      formula = y ~ x,
      parse = TRUE,
      size = 3,
      color = "black"
    ) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.2) +
    theme(
      strip.placement = "outside",
      axis.line.y = element_line(linewidth = 0.35),
      axis.line.x = element_line(linewidth = 0.35),
      axis.ticks.length.y = unit(-0.1, "cm"),
      axis.title.y = element_markdown(size = 8, family = "Helvetica", margin = margin(r = -15)),
      axis.title.x = element_markdown(size = 8, family = "Helvetica", margin = margin(t = +15)),
      text = element_text(family = "Helvetica"),
      axis.text.y = element_text(color = "black", size = 6),
      axis.text.x = element_text(color = "black", size = 6),
      strip.background = element_blank(),
      strip.text.x.bottom = element_markdown(size = 8, margin = margin(t = 1)),
      panel.background = element_rect(fill = "transparent", color = "transparent"),
      plot.background = element_rect(fill = "transparent", color = "transparent"),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  
  return(list(plot = scatter_plot, correlation_test = cor_test))
}


)



  # Generate the plot
  plot <- ggplot(plot_data, aes(x = IncLevelCategory, y = IncLevelDifference_fnf, fill = IncLevelCategory)) +
    geom_hline(yintercept = 0, lty = 2, color = "grey25", linewidth = 0.25) +
    geom_jitter(aes(color = Highlight), width = 0.2, size = 0.5) +
geom_text(
      data = plot_data %>% filter(Highlight == "Highlighted"),
      aes(label = geneSymbol_KD),
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



# Example usage for kd_comparedtoFNF_plot_data
kd_compared_fnf_plot <- generate_scatter_plot_and_stats(
  plot_data = kd_comparedtoFNF_plot_data,
  y_label = "delta PSI in response to Fn-F/PBS",  # Updated y-axis label
  colors = c("#FFB81C", '#005587'),
  highlight_gene = "SNRNP70"  # Highlight SNRNP70
)