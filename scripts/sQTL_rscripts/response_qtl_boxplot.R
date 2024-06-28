library(plotgardener)
setwd("/work/users/s/e/seyoun/CQTL_sQTL/output")
#make a plot for the response QTL for the FNF REML true
response_pbs_results <- readRDS("01.qtltools_re/conditional_pbs/response_pbs_re_lmer.rds")
response_fnf_results <- readRDS("01.qtltools_re/conditional_fnf/response_fnf_re_lmer.rds")

#make a plot for the response QTL for the FNF REML false
response_pbs_results <- readRDS("01.qtltools_re/conditional_pbs/response_pbs_re_no_lmer.rds")
response_fnf_results <- readRDS("01.qtltools_re/conditional_fnf/response_fnf_re_no_lmer.rds")

sig_re_pbs <- response_pbs_results %>%
  dplyr::filter(interaction_pval < 0.05)


sig_re_fnf <- response_fnf_results %>%
  dplyr::filter(interaction_pval < 0.05)

phe_id_fnf_specific <- setdiff(sig_re_fnf$var_id,sig_re_pbs$var_id)
filtered_sig_re_fnf <- sig_re_fnf %>%
  dplyr::filter(phe_id %in% phe_id_fnf_specific)



fnf_with_source <- cbind(sig_re_fnf, source = "fnf")
pbs_with_source <- cbind(sig_re_pbs, source = "pbs")
combined_both_qtl <- rbind(fnf_with_source, pbs_with_source)

combined_both_qtls_sources <- aggregate(source ~ var_id + phe_id, data = combined_both_qtl, 
                          FUN = function(x) c(count = length(x), sources = paste(x, collapse = ", ")))
specific_to_condtions_qtl <- combined_both_qtls_sources[combined_both_qtls_sources$source[,"count"] == 1, c("var_id", "phe_id", "source")]
specific_to_condtions_qtl <- specific_to_condtions_qtl %>%
  mutate(
    count = as.numeric(source[, "count"]),
    sources = source[, "sources"]
  ) %>%
  dplyr::select(-source) 

var_id_fnf_specific <- specific_to_condtions_qtl %>% dplyr::filter(sources  =="fnf") %>% dplyr::select(var_id,phe_id)
fnf_specific_re_subset <- sig_re_fnf %>%
  inner_join(var_id_fnf_specific, by = c("var_id", "phe_id"))
# #make a boxplot 
# pdf(file = "./results_plots/response_qtl/fnf_specific_noREML.pdf",   # The directory you want to save the file in
#     width = 3.5, # The width of the plot in inches
#     height = 4) # The height of the plot in inche

split_indices <- split(1:nrow(fnf_specific_re_subset), ceiling(seq_along(1:nrow(fnf_specific_re_subset)) / 500))

for (chunk_index in seq_along(split_indices)) {
  chunk <- split_indices[[chunk_index]]
  
  pdf(file = paste0("./results_plots/response_qtl/fnf_specific_noREML_part", chunk_index, ".pdf"), width = 3.5, height = 4)
  for (i in chunk) { 
      print(i)
      intronID <- fnf_specific_re_subset$phe_id[i]
      
      psi_matrix <- psi_ratio[rownames(psi_ratio) %in% intronID, ]
      psi_data <- data.frame(sampleID = names(psi_matrix),
                             psi = as.double(psi_matrix))
      
      variantID <- fnf_specific_re_subset$var_id[i]
      geno_data <- all_geno_transpose_df[rownames(all_geno_transpose_df) %in% variantID, ] %>%
        as.data.frame() %>%
        rownames_to_column(var = "sampleID")
      colnames(geno_data)[2] <- "genotype"
      
      meta_data_fi <- geno_data %>%
        left_join(meta_catl_simple,by = c("sampleID"="ID")) %>%
        left_join(psi_data,by = c("sampleID"="sampleID")) %>%
        dplyr::select(c("sampleID","Condition","psi","genotype"))
      
      meta_combined_all <- meta_data_fi %>%
        mutate(
          genotype = factor(genotype, levels = c("0", "1", "2")),
          Condition = ifelse(Condition == "CTL", "PBS", ifelse(Condition == "FNF", "FN-f", NA)),
          Condition = factor(Condition, levels = c("PBS", "FN-f"))
        )
      
      boxplot_fnf_specific <- ggplot(meta_combined_all, aes(x = genotype, y = psi, fill = Condition, color = Condition)) +
        geom_boxplot() +
        geom_point(aes(colour = Condition), position = position_jitterdodge(), alpha = 0.5,size = 1) +
        geom_smooth(aes(group = Condition, color = Condition), method = "lm", se = FALSE, linetype = "solid") +
        labs(x = "Genotype", y = "Ratio of PSI") +
        scale_fill_manual(values = c("PBS" = "#9FCCE4", "FN-f" = "#FAB394")) +
        scale_color_manual(values = c("PBS" = "#3f7d9e", "FN-f" = "#db7e56")) +  # Set color for both "PBS" and "FNF" groups
        theme_minimal() +
        theme(
          aspect.ratio = 1,
          panel.grid = element_blank(),
          axis.line.x = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
          axis.line.y = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
          plot.background = element_rect(fill = 'transparent', color = NA),
          legend.position = "none"
        ) 
      pageCreate(width = 3.5, height =4 , showGuides = FALSE)
      
      
      plotText(label = "FN-f specific",  x = 1.75,
               y = 0.25,
               fontsize = 12, fontfamily = "Helvetica",
               just="center",
               fontcolor = "black",
               fontface = "bold")
      
      plotText(label =paste0("Interaction_pvalue:",round(as.double(fnf_specific_re_subset$interaction_pval[i]),digits=4) ), 
               x = 1.75,
               y=0.5,
               fontsize = 10, fontfamily = "Helvetica",
               just="center",
               fontcolor = "black")
      
      plotText(label =  paste0("PBS-Beta:",round(as.double(fnf_specific_re_subset$PBS_beta[i]),digits=4)), 
               x = 1.75,
               y=0.65,
               fontsize = 10, fontfamily = "Helvetica",
               just="center",
               fontcolor = "black")
      plotText(label =  paste0("FN-f-Beta:",round(as.double(fnf_specific_re_subset$FNF_beta[i]),digits=4)), 
               x = 1.75,
               y=0.8,
               fontsize = 10, fontfamily = "Helvetica",
               just="center",
               fontcolor = "black")
      
      plotGG(plot =  boxplot_fnf_specific, x = 0.25, y = 0.85, height = 3, width = 3)
      
      plotText(label =  variantID, 
               x = 1.75,
               y=3.85,
               fontsize = 10, fontfamily = "Helvetica",
               just="center",
               fontcolor = "black")
      
    }
    dev.off()
}

#shared QTL 
# Combine the two data frames and add a source column
combined_data_sig <- sig_re_fnf %>%
  mutate(source = "fnf") %>%
  bind_rows(sig_re_pbs %>% mutate(source = "pbs"))

# Group by var_id and phe_id, summarize, and filter for counts not equal to 1
filtered_data_shared <- combined_data_sig %>%
  group_by(var_id, phe_id) %>%
  summarize(count = n(), sources = paste(source, collapse = ", "), .groups = 'drop') %>%
  dplyr::filter(count != 1)

# Join the filtered data back with the original combined data to get the relevant rows
combined_data_sig_shared_subset <- combined_data_sig %>%
  semi_join(filtered_data_shared, by = c("var_id", "phe_id"))

#make a boxplot 
split_indices <- split(1:nrow(combined_data_sig_shared_subset), ceiling(seq_along(1:nrow(combined_data_sig_shared_subset)) / 500))

for (chunk_index in seq_along(split_indices)) {
  chunk <- split_indices[[chunk_index]]
  
  pdf(file = paste0("./results_plots/response_qtl/shard_noREML_part", chunk_index, ".pdf"), width = 3.5, height = 4)
  
  
  for (i in chunk) { 
    print(i)
    intronID <- combined_data_sig_shared_subset$phe_id[i]
    
    psi_matrix <- psi_ratio[rownames(psi_ratio) %in% intronID, ]
    psi_data <- data.frame(sampleID = names(psi_matrix),
                           psi = as.double(psi_matrix))
    
    variantID <- combined_data_sig_shared_subset$var_id[i]
    geno_data <- all_geno_transpose_df[rownames(all_geno_transpose_df) %in% variantID, ] %>%
      as.data.frame() %>%
      rownames_to_column(var = "sampleID")
    colnames(geno_data)[2] <- "genotype"
    
    meta_data_fi <- geno_data %>%
      left_join(meta_df,by = c("sampleID"="ID")) %>%
      left_join(psi_data,by = c("sampleID"="sampleID")) %>%
      left_join(covariates_df, by=c("sampleID"= "id")) %>%
      dplyr::select(c("sampleID","Condition","psi","genotype"))
    
    meta_combined_all <- meta_data_fi %>%
      mutate(
        genotype = factor(genotype, levels = c("0", "1", "2")),
        Condition = ifelse(Condition == "CTL", "PBS", ifelse(Condition == "FNF", "FN-f", NA)),
        Condition = factor(Condition, levels = c("PBS", "FN-f"))
      )
    
    shared_sig <- ggplot(meta_combined_all, aes(x = genotype, y = psi, fill = Condition, color = Condition)) +
      geom_boxplot() +
      geom_point(aes(colour = Condition), position = position_jitterdodge(), alpha = 0.5,size = 1) +
      geom_smooth(aes(group = Condition, color = Condition), method = "lm", se = FALSE, linetype = "solid") +
      labs(x = "Genotype", y = "Ratio of PSI") +
      scale_fill_manual(values = c("PBS" = "#9FCCE4", "FN-f" = "#FAB394")) +
      scale_color_manual(values = c("PBS" = "#3f7d9e", "FN-f" = "#db7e56")) +  # Set color for both "PBS" and "FNF" groups
      theme_minimal() +
      theme(
        aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.line.x = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
        axis.line.y = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
        plot.background = element_rect(fill = 'transparent', color = NA),
        legend.position = "none"
      ) 
    pageCreate(width = 3.5, height =4 , showGuides = FALSE)
    
    
    plotText(label = "Shared",  x = 1.75,
             y = 0.25,
             fontsize = 12, fontfamily = "Helvetica",
             just="center",
             fontcolor = "black",
             fontface = "bold")
    
    plotText(label =paste0("Interaction_pvalue:",round(as.double(combined_data_sig_shared_subset$interaction_pval[i]),digits=4) ), 
             x = 1.75,
             y=0.5,
             fontsize = 10, fontfamily = "Helvetica",
             just="center",
             fontcolor = "black")
    
    plotText(label =  paste0("PBS-Beta:",round(as.double(combined_data_sig_shared_subset$PBS_beta[i]),digits=4)), 
             x = 1.75,
             y=0.65,
             fontsize = 10, fontfamily = "Helvetica",
             just="center",
             fontcolor = "black")
    plotText(label =  paste0("FN-f-Beta:",round(as.double(combined_data_sig_shared_subset$FNF_beta[i]),digits=4)), 
             x = 1.75,
             y=0.8,
             fontsize = 10, fontfamily = "Helvetica",
             just="center",
             fontcolor = "black")
    
    plotGG(plot =  shared_sig, x = 0.25, y = 0.85, height = 3, width = 3)
    
    plotText(label =  variantID, 
             x = 1.75,
             y=3.85,
             fontsize = 10, fontfamily = "Helvetica",
             just="center",
             fontcolor = "black")
    
  }
  dev.off()
}

#-------------------------------------------------------------------------------
#pbs-specific

differences <- sig_re_fnf %>%
  mutate(source = "fnf") %>%
  bind_rows(sig_re_pbs %>% mutate(source = "pbs")) %>%
  group_by(var_id, phe_id) %>%
  dplyr::summarize(count = n(), sources = paste(source, collapse = ", "), .groups = 'drop') %>%
  dplyr::filter(count == 1) %>%
  dplyr::select(var_id, phe_id, sources)

var_id_pbs_specific <- differences %>% dplyr::filter(sources  =="pbs") %>% dplyr::select(var_id,phe_id)
pbs_specific_re_subset <- sig_re_pbs %>%
  inner_join(var_id_pbs_specific, by = c("var_id", "phe_id"))



#make a boxplot 
split_indices <- split(1:nrow(pbs_specific_re_subset), ceiling(seq_along(1:nrow(pbs_specific_re_subset)) / 500))

for (chunk_index in seq_along(split_indices)) {
  chunk <- split_indices[[chunk_index]]
  
  pdf(file = paste0("./results_plots/response_qtl/pbs_specific_noREML_part", chunk_index, ".pdf"), width = 3.5, height = 4)
  
  for (i in chunk) { 
  print(i)
  intronID <- pbs_specific_re_subset$phe_id[i]
  
  psi_matrix <- psi_ratio[rownames(psi_ratio) %in% intronID, ]
  psi_data <- data.frame(sampleID = names(psi_matrix),
                         psi = as.double(psi_matrix))
  
  variantID <- pbs_specific_re_subset$var_id[i]
  geno_data <- all_geno_transpose_df[rownames(all_geno_transpose_df) %in% variantID, ] %>%
    as.data.frame() %>%
    rownames_to_column(var = "sampleID")
  colnames(geno_data)[2] <- "genotype"
  
  meta_data_fi <- geno_data %>%
    left_join(meta_catl_simple,by = c("sampleID"="ID")) %>%
    left_join(psi_data,by = c("sampleID"="sampleID")) %>%
    dplyr::select(c("sampleID","Condition","psi","genotype"))
  
  meta_combined_all <- meta_data_fi %>%
    mutate(
      genotype = factor(genotype, levels = c("0", "1", "2")),
      Condition = ifelse(Condition == "CTL", "PBS", ifelse(Condition == "FNF", "FN-f", NA)),
      Condition = factor(Condition, levels = c("PBS", "FN-f"))
    )
  
  boxplot_PBS_specific <- ggplot(meta_combined_all, aes(x = genotype, y = psi, fill = Condition, color = Condition)) +
    geom_boxplot() +
    geom_point(aes(colour = Condition), position = position_jitterdodge(), alpha = 0.5,size = 1) +
    geom_smooth(aes(group = Condition, color = Condition), method = "lm", se = FALSE, linetype = "solid") +
    labs(x = "Genotype", y = "Ratio of PSI") +
    scale_fill_manual(values = c("PBS" = "#9FCCE4", "FN-f" = "#FAB394")) +
    scale_color_manual(values = c("PBS" = "#3f7d9e", "FN-f" = "#db7e56")) +  # Set color for both "PBS" and "FNF" groups
    theme_minimal() +
    theme(
      aspect.ratio = 1,
      panel.grid = element_blank(),
      axis.line.x = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
      axis.line.y = element_line(colour = 'black', size = 0.5, linetype = 'solid'),
      plot.background = element_rect(fill = 'transparent', color = NA),
      legend.position = "none"
    ) 
  pageCreate(width = 3.5, height =4 , showGuides = FALSE)
  
  
  plotText(label = "PBS specific",  x = 1.75,
           y = 0.25,
           fontsize = 12, fontfamily = "Helvetica",
           just="center",
           fontcolor = "black",
           fontface = "bold")
  
  plotText(label =paste0("Interaction_pvalue:",round(as.double(pbs_specific_re_subset$interaction_pval[i]),digits=4) ), 
           x = 1.75,
           y=0.5,
           fontsize = 10, fontfamily = "Helvetica",
           just="center",
           fontcolor = "black")
  
  plotText(label =  paste0("PBS-Beta:",round(as.double(pbs_specific_re_subset$PBS_beta[i]),digits=4)), 
           x = 1.75,
           y=0.65,
           fontsize = 10, fontfamily = "Helvetica",
           just="center",
           fontcolor = "black")
  plotText(label =  paste0("FN-f-Beta:",round(as.double(pbs_specific_re_subset$FNF_beta[i]),digits=4)), 
           x = 1.75,
           y=0.8,
           fontsize = 10, fontfamily = "Helvetica",
           just="center",
           fontcolor = "black")
  
  plotGG(plot =  boxplot_PBS_specific, x = 0.25, y = 0.85, height = 3, width = 3)
  
  plotText(label =  variantID, 
           x = 1.75,
           y=3.85,
           fontsize = 10, fontfamily = "Helvetica",
           just="center",
           fontcolor = "black")
  
}
dev.off()

}













