# This is only for the PBS to draft visualize the sQTL and locus zoom plot
# Response sQTL
## Author: Seyoun Byun

#-------------------------------------------------------------------------------
setwd("/work/users/s/e/seyoun/CQTL_sQTL")
library(tidyverse)
library(tibble)
library(janitor)
library(vcfR)
library(car)
library(lme4)
library(ggtext)
library(data.table)
library(org.Hs.eg.db)
library(qvalue)
library(plotgardener)
library(grid)
library(RColorBrewer)

source("./scripts/sQTL_rscripts/utils.R")

yl_gn_bu <- brewer.pal(n = 9, name = "YlGnBu")
# Signal data 
dir_merged <- "/work/users/s/e/seyoun/CQTL_sQTL/output/signals/merged_norm/"
ctl_signal <- paste0(dir_merged,"CTL_norm.bw")
fnf_signal <- paste0(dir_merged,"FNF_norm.bw")

norm_header <- c("phe_id", "phe_chr", "phe_from", "phe_to", "phe_strd", 
                 "n_var_in_cis", "dist_phe_var", "var_id", "var_chr", "var_from",
                 "var_to", "nom_pval", "r_squared", "slope", "slope_se", "best_hit")

#make a plot for the response QTL
response_pbs_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_pbs_results.rds") 
response_fnf_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_fnf_results.rds")

pbs_norm_qtl_pc5_list <- readRDS("output/nominal_1mb/nominal_pbs/pbs_norm_qtl_pc5_list.rds")
fnf_norm_qtl_pc4_list <- readRDS("output/nominal_1mb/nominal_fnf/fnf_norm_qtl_pc4_list.rds")


load("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/freq/maf_id_add.rds") #maf_id_add  is the name
maf_subset_rsID <- maf_id_add |> dplyr::select(-c("minor_allele","MAF"))

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
#meta #exclude samples for  ['AM7352', 'AM7244']
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






sig_re_pbs <- response_pbs_results %>%
  dplyr::filter(interaction_pval <= 0.05)


sig_re_fnf <- response_fnf_results %>%
  dplyr::filter(interaction_pval <= 0.05)

# This is the high significance 

pbs_highConfi_re <- response_pbs_results %>%
  dplyr::filter(interaction_pval <= 0.05) |>
  dplyr::filter(MAF >= 0.2) |> 
  dplyr::filter(abs(delta_beta) >= 0.2)


fnf_highConfi_re <- response_fnf_results %>%
  dplyr::filter(interaction_pval <= 0.05) |>
  dplyr::filter(MAF >= 0.2) |> 
  dplyr::filter(abs(delta_beta) >= 0.2)



phe_id_fnf_specific <- setdiff(fnf_highConfi_re$phe_id,pbs_highConfi_re$phe_id)
filtered_sig_re_fnf <- fnf_highConfi_re %>%
  dplyr::filter(phe_id %in% phe_id_fnf_specific)


pbs_with_source <- cbind(pbs_highConfi_re, source = "pbs")
fnf_with_source <- cbind(fnf_highConfi_re, source = "fnf")

combined_both_qtl <- rbind(pbs_with_source, fnf_with_source)

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



#-------------------------------------------------------------------------------
#pbs-specific

differences <- fnf_highConfi_re %>%
  mutate(source = "fnf") %>%
  bind_rows(pbs_highConfi_re %>% mutate(source = "pbs")) %>%
  group_by(var_id, phe_id) %>%
  dplyr::summarize(count = n(), sources = paste(source, collapse = ", "), .groups = 'drop') %>%
  dplyr::filter(count == 1) %>%
  dplyr::select(var_id, phe_id, sources)

var_id_pbs_specific <- differences %>% dplyr::filter(sources  =="pbs") %>% dplyr::select(var_id,phe_id)
pbs_specific_re_subset <- pbs_highConfi_re %>%
  inner_join(var_id_pbs_specific, by = c("var_id", "phe_id"))



#make a boxplot 
split_indices <- split(1:nrow(pbs_specific_re_subset), ceiling(seq_along(1:nrow(pbs_specific_re_subset)) / 500))

for (chunk_index in seq_along(split_indices)) {
  chunk <- split_indices[[chunk_index]]
  
  pdf(file = paste0("output/results_plots/response_qtl/highconfidence_pbs_part", chunk_index, ".pdf"), width = 4.5, height = 8.5)
  for (i in chunk) { 
    print(i)
    intronID <- pbs_specific_re_subset$phe_id[i]
    
    psi_matrix <- psi_ratio[rownames(psi_ratio) %in% intronID, ]
    psi_data <- data.frame(sampleID = names(psi_matrix),
                           psi = as.double(psi_matrix))
    
    variantID <- pbs_specific_re_subset$var_id[i]
    rsID <- pbs_specific_re_subset$rsID[i]
    minregion <- pbs_specific_re_subset$var_from[i] -100000
    maxregion  <-  pbs_specific_re_subset$var_from[i] +100000
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
      theme(strip.placement = "outside",
            axis.line = element_line(linewidth = 0.25),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(color = "black", linewidth = 0.25),
            axis.ticks.length.y = unit(-0.1, "cm"),
            axis.title.y = element_markdown(size = 6, family= "Helvetica",
                                            margin = margin(r = -0.5)),
            axis.title.x = element_markdown(size = 6, family= "Helvetica",
                                            margin = margin(r = -0.5)),
            text = element_text(family = "Helvetica"),
            axis.text.x = element_text(color = "black", size = 6),
            axis.text.y = element_text(color = "black", size = 6, margin = margin(b = 0)),
            strip.background = element_blank(),
            strip.text.x.top = element_text(size = 8, margin = margin(b = 5)),
            panel.background = element_rect(fill = "transparent", color = "transparent"),
            plot.background = element_rect(fill = "transparent", color = "transparent"),
            panel.grid = element_blank(),
            legend.position = "none",
            panel.spacing.y = unit(0.5, "cm"))
    
    chrom <- pbs_specific_re_subset$phe_chr[i]
    ld_file <- paste0("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld/",chrom,"/",variantID,".ld")
    ld_calc <- fread(ld_file) |> dplyr::select(c("SNP_B","R2")) |>
      dplyr::rename("var_id" = "SNP_B")
    
    
    pbs_norminal_qtl <- pbs_norm_qtl_pc5_list[[chrom]]
    pbs_qtl_region <- pbs_norminal_qtl |> dplyr::filter(phe_id %in% intronID)
    pbs_qtl_pval_rsid <- left_join(pbs_qtl_region,maf_subset_rsID  ,by="var_id")
    
    fnf_norminal_qtl <- fnf_norm_qtl_pc4_list[[chrom]]
    fnf_qtl_region <- fnf_norminal_qtl |> dplyr::filter(phe_id %in% intronID)
    fnf_qtl_pval_rsid <- left_join(fnf_qtl_region,maf_subset_rsID  ,by="var_id")
    
    
    leftjoin_pbs <- inner_join(pbs_qtl_pval_rsid, ld_calc, by = c("var_id" = "var_id"))
    leftjoin_fnf <- inner_join(fnf_qtl_pval_rsid, ld_calc, by = c("var_id" = "var_id"))
    
    #-------------------------------------------------------------------------
    #Plotting data
    pbs_ld <- leftjoin_pbs |>  mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))
    pbs_ld$LDgrp <- addNA(pbs_ld$LDgrp)
    
    pbs_locus_plot <- pbs_ld |> dplyr::rename(chrom ="var_chr",
                                              pos = "var_from",
                                              p="nom_pval",
                                              LD ="R2",
                                              snp ="rsID") |> data.frame()
    
    pbs_locus_plot <- pbs_locus_plot %>%
      dplyr::select("chrom", "pos", "p", "snp", "LD","LDgrp","phe_id","var_id")
    
    
    fnf_ld <- leftjoin_fnf |>   mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))
    fnf_ld$LDgrp <- addNA(fnf_ld$LDgrp)
    
    fnf_locus_plot <- fnf_ld |>  dplyr::rename(chrom ="var_chr",
                                               pos = "var_from",
                                               p="nom_pval",
                                               LD ="R2",
                                               snp ="rsID") |> data.frame()
    
    
    fnf_locus_plot <- fnf_locus_plot %>%
      dplyr::select("chrom", "pos", "p", "snp", "LD","LDgrp","phe_id","var_id")
    
    #------------------------------------------------------------------------
    #plotting
    
    pageCreate(width = 4.5, height =8.5 , showGuides = FALSE)
    
    region_pg <- pgParams(assembly = "hg38",chrom = chrom,
                          chromstart = minregion,
                          chromend = maxregion,
                          x = 0.55, width = 3.75)
    pbs_ylim <- ceiling(max(log10(pbs_locus_plot$p)*-1)) + 2
    fnf_ylim <- ceiling(max(log10(fnf_locus_plot$p)*-1)) + 2
    ylim_pg <- max(pbs_ylim,fnf_ylim)
    
    pbs_locus <- plotManhattan(data = pbs_locus_plot ,
                               params = region_pg,
                               range = c(0, ylim_pg),
                               fill = colorby("LDgrp",
                                              palette = colorRampPalette(c("#262C74",
                                                                           "#98CDED",
                                                                           "#499A53",
                                                                           "#EEA741",
                                                                           "#DD3931",
                                                                           "grey"))),
                               y = 0.5, height = 1.5,
                               snpHighlights = data.frame(snp =rsID,
                                                          pch = c(24),
                                                          cex = c(0.75),
                                                          col = c("black")))
    annoYaxis(plot = pbs_locus, at = seq(0, ylim_pg, 2),
              axisLine = TRUE, fontsize = 8)
    plotText(
      label = "-log10(p-value)", x = 0.24, y = 1.25, rot = 90,
      fontsize = 8, just = "center",
      default.units = "inches", fontfamily = "Helvetica"
    )
    plotText(label = "PBS sQTL", x = 0.65, y = 0.5, just = c("left", "top"),
             fontfamily = "Helvetica", fontsize = 11)
    
    
    plotLegend(legend = c("0.8 - 1.0",
                          "0.6 - 0.8",
                          "0.4 - 0.6",
                          "0.2 - 0.4",
                          "0.0 - 0.2"),
               fill = c("#DD3931", "#EEA741", "#499A53","#98CDED","#262C74"),
               x = 3.5, y = 0.5, width = 0.1, height = 0.4, border = FALSE,
               fontsize = 6)
    
    #fnf plot
    fnf_locus <- plotManhattan(data = fnf_locus_plot ,
                               params = region_pg,
                               range = c(0, ylim_pg),
                               fill = colorby("LDgrp",
                                              palette = colorRampPalette(c("#262C74",
                                                                           "#98CDED",
                                                                           "#499A53",
                                                                           "#EEA741",
                                                                           "#DD3931",
                                                                           "grey"))),
                               y = 2.2, height = 1.5,
                               snpHighlights = data.frame(snp =rsID,
                                                          pch = c(24),
                                                          cex = c(0.75),
                                                          col = c("black")))
    
    
    annoYaxis(plot = fnf_locus, at = seq(0, ylim_pg, 2),
              axisLine = TRUE, fontsize = 8)
    plotText(
      label = "-log10(p-value)", x = 0.24, y = 2.95, rot = 90,
      fontsize = 8, just = "center",
      default.units = "inches", fontfamily = "Helvetica"
    )
    plotText(label = "FN-f sQTL", x = 0.6, y = 2.1, just = c("left", "top"),
             fontfamily = "Helvetica", fontsize = 11)
    
    
    grid.points(x = 3.5, y = 2.4, default.units = "native", pch = 24,
                size = unit(0.75, "char"))
    plotText(label = rsID,
             fontsize = 8, fontfamily = "Helvetica",
             just = "left", x = 3.6, y = 2.4)
    
    
    RNA_signals <- plotMultiSignal(data = list(ctl_signal,
                                               fnf_signal),
                                   params = region_pg,
                                   y = 3.85, height = 0.74, linecolor = c(yl_gn_bu[3],yl_gn_bu[6]), 
                                   fill = c(yl_gn_bu[3],yl_gn_bu[6]),
                                   default.units = "inches",
                                   gapdistance = 0.02)
    plotText(label = "PBS",
             fontsize = 7, x = 0.5, y =3.9 , just = "left", fontfamily = "Helvetica")
    plotText(label = "FN-f",
             fontsize = 7, x = 0.5, y =4.25, just = "left", fontfamily = "Helvetica")
    
    plotText(
      label = "RNA", x = 0.35, y = 4.2, rot = 90,
      fontsize = 8, just = "center",
      default.units = "inches", fontfamily = "Helvetica"
    )
    
    gene_hl <- pbs_specific_re_subset$SYMBOL[i]
    plotgenes <- plotGenes(params = region_pg, y = 4.7,
                           height = 0.5,
                           geneHighlights = data.frame("gene" = gene_hl,
                                                       "color" = "#37a7db"))
    annoGenomeLabel(plot = plotgenes, params = region_pg, fontsize = 8,y=5.25)
    
    
    annoHighlight(
      plot = fnf_locus,
      chrom = pbs_specific_re_subset$phe_chr[i],
      chromstart = pbs_specific_re_subset$phe_from[i],
      chromend = pbs_specific_re_subset$phe_to[i],
      y = 3.85, height = 1.4, just = c("left", "top"),
      default.units = "inches"
    )
    
    #-------------------------------------------------------------------------
    
    
    plotText(label = "shared highConf",  x = 0.5,
             y = 0.25,
             fontsize = 12, fontfamily = "Helvetica",
             just="left",
             fontcolor = "black",
             fontface = "bold")
    plotText(label = paste0("rank:",pbs_specific_re_subset$rank[i]),  x = 1.75,
             y = 0.25,
             fontsize = 10, fontfamily = "Helvetica",
             just="center",
             fontcolor = "black")
    
    plotText(label =paste0("Interaction_pvalue:",round(as.double(pbs_specific_re_subset$interaction_pval[i]),digits=4) ), 
             x = 0.5,
             y=5.5,
             fontsize = 8, fontfamily = "Helvetica",
             just="left",
             fontcolor = "black")
    
    plotText(label =  paste0("PBS-Beta:",round(as.double(pbs_specific_re_subset$PBS_beta[i]),digits=4)), 
             x = 0.5,
             y=5.65,
             fontsize = 8, fontfamily = "Helvetica",
             just="left",
             fontcolor = "black")
    plotText(label =  paste0("FN-f-Beta:",round(as.double(pbs_specific_re_subset$FNF_beta[i]),digits=4)), 
             x = 0.5,
             y=5.8,
             fontsize = 8, fontfamily = "Helvetica",
             just="left",
             fontcolor = "black")
    
    plotGG(plot =  boxplot_fnf_specific, x = 0, y = 6, height = 2.5, width = 3)
    
    plotText(label =intronID, 
             x = 0.5,
             y=5.95,
             fontsize =8, fontfamily = "Helvetica",
             just="left",
             fontcolor = "black")
    
  }
  dev.off()
}



