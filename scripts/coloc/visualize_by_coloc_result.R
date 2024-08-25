#FNF

coloc_fnf_sig <- coloc_results_fnf.df[coloc_results_fnf.df$coloc_H4 > 0.7,]
gene_name <- "GCAT"


test_condtion_ma <- coloc_fnf_sig
Intron_test <- coloc_fnf_sig$phe_id |>unique()

intron_test <- Intron_test[1]

for(intron_test in Intron_test){
psi_matrix <- psi_ratio[rownames(psi_ratio) %in% intron_test, ]
psi_data <- data.frame(sampleID = names(psi_matrix),
                       psi = as.double(psi_matrix))



variantID <- test_condtion_ma |> dplyr::filter(phe_id == intron_test) |> dplyr::pull("sQTL_var_id")
rsID <- test_condtion_ma |> dplyr::filter(phe_id == intron_test) |> dplyr::pull("sQTL_rsID")
minor_allele <- test_condtion_ma |> dplyr::filter(phe_id == intron_test) |> dplyr::pull("sQTL_minor_allele")

alleles <- unlist(strsplit(variantID, ":"))
ref_allele <- alleles[3]
alt_allele <- alleles[4]
var_pos <-  alleles[2] |> as.numeric()

protective_allele <- ifelse(minor_allele == ref_allele, alt_allele, ref_allele)
minregion <- var_pos - 100000
maxregion  <-  var_pos + 100000
geno_data <- all_geno_transpose_df[rownames(all_geno_transpose_df) %in% variantID, ] %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
colnames(geno_data)[2] <- "genotype"

meta_data_fi <- geno_data %>%
  left_join(meta_catl_simple,by = c("sampleID"="ID")) %>%
  left_join(psi_data,by = c("sampleID"="sampleID")) %>%
  dplyr::select(c("sampleID","Condition","psi","genotype"))

meta_combined_all <-  meta_data_fi %>%
  mutate(
    genotype = factor(genotype, levels = c("0", "1", "2")),
    Condition = ifelse(Condition == "CTL", "PBS", ifelse(Condition == "FNF", "FN-f", NA)),
    Condition = factor(Condition, levels = c("PBS", "FN-f")),
    genotype = case_when(
      genotype == "2" ~ paste(minor_allele, minor_allele, sep = "/"),
      genotype == "1" ~ paste(minor_allele, protective_allele, sep = "/"),
      genotype == "0" ~ paste(protective_allele, protective_allele, sep = "/")
    )
  )  %>%
  mutate(
    genotype = factor(genotype, levels = c(
      paste(protective_allele, protective_allele, sep = "/"),
      paste(minor_allele, protective_allele, sep = "/"),
      paste(minor_allele, minor_allele, sep = "/")
    ))
  )
maxrange <- range(meta_combined_all$psi)
yaxis_range_min <- min(0,maxrange[1])
yaxis_range_max <- round(maxrange[2],1)

genotype_plot <- ggplot(meta_combined_all, aes(x = genotype, y = psi, fill = Condition)) +
  geom_boxplot(outlier.shape = NA,
               linewidth = 0.25, alpha = 0.7) +
  #geom_jitter(width = 0.1, color = "grey40", size = 0.25)+
  geom_point(color = "grey40", position = position_jitterdodge(),size = 0.25) +
  geom_smooth(aes(group = Condition, color = Condition), method = "lm", se = FALSE, linetype = "solid") +
  labs(x = "Genotype", y = "Intron usage") +
  scale_fill_manual(values = c('#1e87a5','#FFB81C'))+
  scale_color_manual(values = c('#1e87a5','#FFB81C'))+
  scale_y_continuous( limits = c(yaxis_range_min, yaxis_range_max), breaks = seq(yaxis_range_min, yaxis_range_max, length.out=5)) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
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
        strip.text.x.top = element_text(size = 8, margin = margin(b = 5)),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        legend.position = c(0.1,0.9),
        panel.spacing.y = unit(0.5, "cm"))


chrom <- response_pbs_results |> dplyr::filter(SYMBOL == gene_name) |> dplyr::pull("phe_chr")
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

pageCreate(width = 4.5, height =10 , showGuides = FALSE)

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


grid.points(x = 3.5, y = 2.55, default.units = "native", pch = 24,
            size = unit(0.75, "char"))
plotText(label = rsID,
         fontsize = 8, fontfamily = "Helvetica",
         just = "left", x = 3.6, y = 2.4)
plotText(label = variantID,
         fontsize = 8, fontfamily = "Helvetica",
         just = "left", x = 3.6, y = 2.55)


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
#plotText(label = "OA",
#         fontsize = 7, x = 0.5, y =4.65, just = "left", fontfamily = "Helvetica")
plotText(
  label = "RNA", x = 0.35, y = 4.2, rot = 90,
  fontsize = 8, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)

gene_hl <- gene_name
#plotgenes <- plotGenes(params = region_pg, y = 4.7,
#                       height = 0.5,
#                       geneHighlights = data.frame("gene" = gene_hl,
#                                                   "color" = "#37a7db"))

trans_highlights <- data.frame(
  transcript = c("ENST00000259939.4","ENST00000486622.1"),
  color = c("grey40", "#37a7db")
)

plotTranscripts(
  params = region_pg,
  assembly = "hg38",
  height = 0.5,
  y=4.7,
  just = c("left", "top"), default.units = "inches",
  transcriptHighlights=trans_highlights,
  labels = "both",spaceHeight =0
)
annoGenomeLabel(plot = plotgenes, params = region_pg, fontsize = 8,y=5.25)

hl_start <- response_pbs_results |> dplyr::filter(SYMBOL == gene_name) |> dplyr::pull("phe_from")
hl_end <- response_pbs_results |> dplyr::filter(SYMBOL == gene_name) |> dplyr::pull("phe_to")
anno_highlight_plot <- annoHighlight(
  plot = fnf_locus,
  chrom =chrom,
  chromstart = hl_start,
  chromend = hl_end,
  y = 3.85, height = 1.4, just = c("left", "top"),
  default.units = "inches"
)

# zoom in 
annoZoomLines(
  plot = fnf_locus,   
  chrom =chrom,
  chromstart = hl_start,
  chromend = hl_end,
  y0 = hl_start, x1 = c(3, 4.5), y1 = hl_end, extend = c(5.5, 5.5),
  default.units = "inches",
  lty = 3
)

plotGG(plot =  genotype_plot, x = 0.1, y = 6.5, height = 2.5, width = 3)

plotText(label =paste0("Interaction_pvalue:", round(response_pbs_results |> dplyr::filter(SYMBOL == gene_name) |> dplyr::pull("interaction_pval"),digits=4) ), 
         x = 0.5,
         y=6,
         fontsize = 8, fontfamily = "Helvetica",
         just="left",
         fontcolor = "black")

plotText(label =  paste0("PBS-Beta:",round(response_pbs_results |> dplyr::filter(SYMBOL == gene_name) |> dplyr::pull("PBS_beta"),digits=4)), 
         x = 0.5,
         y=6.2,
         fontsize = 8, fontfamily = "Helvetica",
         just="left",
         fontcolor = "black")
plotText(label =  paste0("PBS-pvalue:",response_pbs_results |> dplyr::filter(SYMBOL == gene_name) |> dplyr::pull("PBS_p")), 
         x = 1.5,
         y=6.2,
         fontsize = 8, fontfamily = "Helvetica",
         just="left",
         fontcolor = "black")

plotText(label =  paste0("FN-f-Beta:",round(response_pbs_results |> dplyr::filter(SYMBOL == gene_name) |> dplyr::pull("FNF_beta"),digits=4)), 
         x = 0.5,
         y=6.35,
         fontsize = 8, fontfamily = "Helvetica",
         just="left",
         fontcolor = "black")
plotText(label =  paste0("FN-f-pvalue",response_pbs_results |> dplyr::filter(SYMBOL == gene_name) |> dplyr::pull("FNF_p")), 
         x = 1.5,
         y=6.35,
         fontsize = 8, fontfamily = "Helvetica",
         just="left",
         fontcolor = "black")

hl_start <- response_pbs_results |> dplyr::filter(SYMBOL == gene_name) |> dplyr::pull("phe_to") -1000
hl_end <-  response_pbs_results |> dplyr::filter(SYMBOL == gene_name) |> dplyr::pull("phe_to") +1000

highlight_pg <- pgParams(assembly = "hg38",chrom = chrom,
                         chromstart = hl_start,
                         chromend = hl_end,
                         x = 3, width = 1.5)

RNA_signals <- plotMultiSignal(data = list(ctl_signal,
                                           fnf_signal),
                               params = highlight_pg,
                               y = 5.5, height = 1.05, linecolor = c(yl_gn_bu[3],yl_gn_bu[6]), 
                               fill = c(yl_gn_bu[3],yl_gn_bu[6]),
                               default.units = "inches",
                               gapdistance = 0.1,binCap = FALSE)

plot_trans_hl <- plotTranscripts(
  params = highlight_pg,
  assembly = "hg38",
  height = 0.5,
  y=6.6,
  just = c("left", "top"), default.units = "inches",
  transcriptHighlights=trans_highlights,
  labels = "both",spaceHeight =0
)
annoGenomeLabel(plot = plot_trans_hl, params = highlight_pg, fontsize = 8,y=7.25,scale = "Mb")


# This is transcript level. 
degene_transcript <- load("output/quant/differential_transcript_expression_dds.rda")
vsd_trans <- fread("output/quant/normalized_vst_transcript.txt")
vsd_gene <- fread("output/quant/normalized_vst_gene.txt")
vsd_rnf144b <- vsd_trans |> dplyr::filter(ENST == "ENST00000486622.1")
vsd_rnf144b_gene <- vsd_gene |> dplyr::filter(ENSG == "ENSG00000137393.10")

vsd_rnf144b_gene_long <- vsd_rnf144b_gene %>%
  pivot_longer(cols = starts_with("AM"),
               names_to = "sampleID",
               values_to = "gene_expression")

vsd_long <- vsd_rnf144b %>%
  pivot_longer(cols = starts_with("AM"),
               names_to = "sampleID",
               values_to = "value")

# Separate groups and calculate means
mean_CTL <- vsd_long %>%
  filter(grepl("CTL", sample)) %>%
  summarise(mean_CTL = mean(value, na.rm = TRUE))

mean_FNF <- vsd_long %>%
  filter(grepl("FNF", sample)) %>%
  summarise(mean_FNF = mean(value, na.rm = TRUE))

# Combine the results
result <- bind_cols(mean_CTL, mean_FNF)


meta_trans <- left_join(meta_combined_all, vsd_long, by="sampleID")

meta_trans_gene <- left_join(meta_trans, vsd_rnf144b_gene_long, by="sampleID")

ggplot(meta_trans_gene, aes(x = genotype, y = value, fill = Condition)) +
  geom_boxplot(outlier.shape = NA,
               linewidth = 0.25, alpha = 0.7) +
  #geom_jitter(width = 0.1, color = "grey40", size = 0.25)+
  geom_point(color = "grey40", position = position_jitterdodge(),size = 0.25) +
  geom_smooth(aes(group = Condition, color = Condition), method = "lm", se = FALSE, linetype = "solid") +
  labs(x = "Genotype", y = "RNF144B_transcript(ENST00000486622.1) normalized counts") +
  scale_fill_manual(values = c('#1e87a5','#FFB81C'))+
  scale_color_manual(values = c('#1e87a5','#FFB81C'))+
  #scale_y_continuous( limits = c(yaxis_range_min, yaxis_range_max), breaks = seq(yaxis_range_min, yaxis_range_max, length.out=5)) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
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
        strip.text.x.top = element_text(size = 8, margin = margin(b = 5)),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        #legend.position = "none",
        panel.spacing.y = unit(0.5, "cm"))





