plotText(label =   paste0("\u03B2 = ", round(fnf_results$FNF_beta, 3)), 
         x = 2.5,
         y=4.75,
         fontsize = 7, fontfamily = "Helvetica",
         just="left",
         fontcolor = "black")


plotText(label =  paste0("pvalue = ", round(fnf_results$PBS_p, 3)), 
         x = 0.6,
         y=4.85,
         fontsize = 7, fontfamily = "Helvetica",
         just="left",
         fontcolor = "black")


plotText(label =   paste0("pvalue = ", round(fnf_results$FNF_p, 3)), 
         x = 2.5,
         y=4.85,
         fontsize = 7, fontfamily = "Helvetica",
         just="left",
         fontcolor = "black")

Plot_specific_gene <- plotGenes(chrom =chrom,
                                chromstart = gcat_start,
                                chromend = gcat_end,
                                y = 7.1,
                                height = 0.3,
                                x=0.55,
                                width=3.5,
                                just = c("left", "top"), default.units = "inches",
                                geneHighlights = data.frame("gene" = gene_hl,
                                                            "color" = "#37a7db"))

#trans_highlights <- data.frame(
#  transcript = c("ENST00000259939.4","ENST00000486622.1"),
#  color = c("grey40", "#37a7db")
#)

#plotTranscripts(
#  params = region_pg,
#  assembly = "hg38",
#  height = 0.5,
#  y=4.7,
#  just = c("left", "top"), default.units = "inches",
#transcriptHighlights=trans_highlights,
#  labels = "both",spaceHeight =0
#)


#hl_start <- fnf_results$phe_from
#hl_end <- fnf_results$phe_to

gcat_start <- 37807934
gcat_end <- 37816897
plotTranscripts(
  chrom =chrom,
  chromstart = gcat_start,
  chromend = gcat_end,
  assembly = "hg38",
  height = 1.5,
  y=5.5,
  x=0.55,
  width=3.5,
  just = c("left", "top"), default.units = "inches",
  #transcriptHighlights=trans_highlights,
  labels = "both",spaceHeight =0
)




# zoom in 
#annoZoomLines(
#  plot = fnf_locus,   
#  chrom =chrom,
#  chromstart = hl_start,
#  chromend = hl_end,
#  y0 = hl_start, x1 = c(3, 4.5), y1 = hl_end, extend = c(5.5, 5.5),
#  default.units = "inches",
#  lty = 3
#)



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

hl_start <- fnf_results$phe_from -50
hl_end <-  fnf_results$phe_from + 50

highlight_pg <- pgParams(assembly = "hg38",chrom = chrom,
                         chromstart = gcat_start,
                         chromend = gcat_end,
                         x = 0.55, width = 3.5)

RNA_signals <- plotMultiSignal(data = list(ctl_signal,
                                           fnf_signal),
                               params = highlight_pg,
                               chrom = chrom, 
                               y = 7.7, height = 1.05, linecolor = c(yl_gn_bu[7],yl_gn_bu[7]), 
                               fill = c(yl_gn_bu[7],yl_gn_bu[7]),
                               default.units = "inches",
                               gapdistance = 0.1,binCap = FALSE)

anno_highlight_plot <- annoHighlight(
  plot = RNA_signals,
  chrom = chrom, 
  chromstart = hl_start,
  chromend = hl_end,
  y = 5.5, height = 1.5, just = c("left", "top"),
  default.units = "inches"
)

#plot_trans_hl <- plotTranscripts(
#  params = highlight_pg,
#  assembly = "hg38",
#  height = 0.5,
#  y=6.6,
#  just = c("left", "top"), default.units = "inches",
#  transcriptHighlights=trans_highlights,
#  labels = "both",spaceHeight =0
#)
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



meta_combined_all_clusters %>%
  group_by(Condition, genotype) %>%
  summarize(mean_psi = mean(psi, na.rm = TRUE),
            .groups = 'drop')



