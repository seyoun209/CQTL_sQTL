# Supplementary Figure 2 - Differential PBS vs OA
## Author: Seyoun Byun
##Date:06.17.2024
#-------------------------------------------------------------------------------
setwd("/work/users/s/e/seyoun/CQTL_sQTL")
library(magrittr)
library(data.table)
library(dplyr)
library(colorspace)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(RColorBrewer)
library(plotgardener)
library(grid)
library(ggrepel)
library(leafcutter)
library(gridExtra)
library(grid)
library(rrvgo)
library(stringr)
library(colorspace)
library(gprofiler2)
source("scripts/utils/utils.R")

#-------------------------------------------------------------------------------
#Heatmap
meta_cqtl <- fread("output/clu_fnf/meta_cqtl")
meta_cqtl_subset <- meta_cqtl[!meta_cqtl$Donor %in% config$samples_to_omit, ]
meta_ctl_oa <- meta_cqtl_subset %>% dplyr::filter(Condition %in% c('CTL','OA')) %>%
  dplyr::select("ID","Age","Predicted_Ancestry","Sex", "Condition")
rownames(meta_ctl_oa) <- meta_ctl_oa$ID
#meta_ctl_oa <- meta_ctl_oa %>% dplyr::select(-"ID")

#Subset delta PSI 15% diff and FDR < 0.05
ratios_oa_qc <- read.table("output/clu_oa/ratio_oa.txt",sep='\t',header=T) |> as.data.frame()
rownames(ratios_oa_qc) <- ratios_oa_qc$Junction 
ratios_oa_qc <- ratios_oa_qc[,-1]
OA_deltapsiCalc_df <- calculate_delta_psi(ratios_oa_qc, "CTL", "OA")
introns_oa_pval_include <- join_introns_deltapsi_fdr(OA_deltapsiCalc_df,introns_oa,"./output/clu_oa/ctlvsoa_ds_cluster_significance.txt")
#save(introns_oa_pval_include, file= "./output/clu_oa/introns_oa_joinAll")
load("./output/clu_oa/introns_oa_joinAll")
introns_oa_sig <- introns_oa_pval_include %>% dplyr::filter(p.adjust < 0.05)
introns_oa_subset_sig <- introns_oa_sig[abs(introns_oa_sig$deltapsi_batch) > 0.15,]
psi_oa_sig_15 <- ratios_oa_qc[which(rownames(ratios_oa_qc) %in% introns_oa_subset_sig$phe_id ),]



cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_subset_norm <- t(apply(psi_oa_sig_15, 1, cal_z_score))

#z-score
brks <- seq(min(data_subset_norm, na.rm = TRUE), max(data_subset_norm, na.rm = TRUE), length = 51)
brks <- seq(-2,2,length.out=50) 

# Define the age group intervals
age_breaks <- c(30, 40, 50, 60, 70, 80, Inf)

# Labels for the age groups
age_labels <- c("30-39", "40-49", "50-59", "60-69", "70-79", "80+")



# Create age groups using the cut function
age_groups <- cut(meta_ctl_oa$Age, breaks = age_breaks, labels = age_labels, right = FALSE)
meta_ctl_oa$Age_range <- age_groups

age_palette <- c("#ffcba4","#cca283","#997a62","#806652","#665142","#332921")

# Create a named vector of colors for the age groups
age_colors <- setNames(age_palette, age_labels)


my_colour = list(
  Condition  = c("CTL" = "#9FCCE4", "OA" = "#FAB394"),
  Sex  = c( "M" =  "#0075B0", "F" = "#F8B7CD"),
  Predicted_Ancestry= c("AMR"="#F5BC9F", "EUR"="#FAF1D2" ,"AFR"="#86CBB7","SAS"="#6B7EA4"),
  Age_range = age_colors)

pheatmap::pheatmap(data_subset_norm,
                   annotation_col = meta_ctl_oa,
                   annotation_colors = my_colour,
                   show_rownames = FALSE,
                   clustering_method = "complete",
                   show_colnames = FALSE,
                   treeheight_row = 0,
                   treeheight_col = 0,
                   fontsize = 5,
                   width = 1,
                   height = 1,
                   #cutree_rows = 2,
                   #cutree_cols = 2,
                   #labels_row = labels,
                   color = colorRampPalette(c("#097EA4", "black", "#BFA527"))(50),
                   breaks = brks
)


heatmapColors <- colorRampPalette(c("#097EA4", "black", "#BFA527"))(5)

ID <- colnames(data_subset_norm)
meta_ctl_oa <- meta_ctl_oa %>% arrange(factor(ID, levels=ID)) %>%
  arrange(Condition)

sample_info <- meta_ctl_oa[, c("Age_range","Predicted_Ancestry","Sex", "Condition")]
colnames(sample_info) <- c("Age","Ancestry","Sex","Condition")

my_colour = list(
  Condition  = c("CTL" = "#B8B8B8", "OA" = "grey20"),
  Sex  = c("F" = "#DD8492", "M" = "#4788BA"),
  Ancestry= c("AFR"="#C74A53", "AMR"="#EFB06E","EUR"= "#5BAD58", "SAS"="#177F97"),
  Age = age_colors)



colAnn <- HeatmapAnnotation(df = sample_info,
                            col = my_colour,
                            gap = unit(0.5, 'mm'),
                            annotation_name_gp= gpar(fontsize = 7),
                            simple_anno_size = unit(3, "mm"),
                            annotation_legend_param  = list(
                              Age = list(
                                title = "Age",
                                title_position = "leftcenter",
                                title_gp = gpar(fontsize = 7),
                                labels_gp = gpar(fontsize = 6),
                                grid_width = unit(2, "mm"),
                                grid_height = unit(1, "mm"),
                                at = c("30-39","40-49","50-59","60-69","70-79","80+"),
                                labels = c("30","40","50","60","70","80+"),
                                ncol = 1
                              ),
                              Condition  = list(
                                title = "Condition",
                                title_gp = gpar(fontsize = 7),
                                labels_gp = gpar(fontsize = 6),
                                title_position = "leftcenter",
                                grid_width = unit(2, "mm"),
                                grid_height = unit(1, "mm"),
                                at = c("CTL", "FNF"),
                                labels = c("PBS","OA"),
                                ncol = 1
                              ),
                              Sex = list(
                                title = "Sex",
                                title_position = "leftcenter",
                                title_gp = gpar(fontsize = 7),
                                labels_gp = gpar(fontsize = 6),
                                grid_width = unit(2, "mm"),
                                grid_height = unit(1, "mm"),
                                at = c("M","F"),
                                labels = c("Male","Female"),
                                ncol = 1
                              ),
                              Ancestry= list(
                                title = "Ancestry",
                                title_position = "leftcenter",
                                title_gp = gpar(fontsize = 7),
                                labels_gp = gpar(fontsize = 6),
                                grid_width = unit(2, "mm"),
                                grid_height = unit(1, "mm"),
                                at = c("AFR","AMR","EUR","SAS"),
                                labels = c("AFR","AMR","EUR","SAS"),
                                ncol = 1
                                
                              )))


column_order <- meta_ctl_oa$ID
data_subset_norm_ordered <- data_subset_norm[, column_order,drop = FALSE]

hmap <- Heatmap(
  data_subset_norm_ordered,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,  # Disable column clustering
  column_order = column_order,  # Specify column order
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  row_dend_reorder = FALSE,
  column_dend_reorder = FALSE,
  top_annotation = colAnn,
  col = colorRamp2(seq(-2, 2), heatmapColors)
)



heatmapLegend <- Legend(at = c(-2, 2),
                        col_fun = colorRamp2(breaks = seq(-2, 2),
                                             colors = heatmapColors),
                        border = NA,
                        title_gp = gpar(fontsize = 0),
                        labels_gp = gpar(fontfamily = "Helvetica",
                                         fontsize = 8),
                        legend_width = unit(4.325, "in"),
                        grid_height = unit(0.11, "in"),
                        direction = "horizontal")

heatmapGrob <- grid.grabExpr(draw(hmap,
                                  show_annotation_legend = FALSE,
                                  show_heatmap_legend = FALSE,
                                  background = "transparent"))
heatmapLegendGrob <- grid.grabExpr(draw(heatmapLegend))

save(heatmapGrob, file = "output/results_plots/Supplementary_figures/supp_fig2_heatmapOA.rda")
save(heatmapLegendGrob, file = "output/results_plots/Supplementary_figures/supp_fig2_heatmapOALegendGrob.rda")

#Venn diagram with the same Directional and the high confidence. 
library(VennDiagram) 
library(ggvenn)
load("./output/clu_fnf/introns_fnf_joinAll")

introns_fnf_sig <- introns_fnf_pval_include %>% dplyr::filter(p.adjust < 0.05) %>%
  mutate(loc = str_split(phe_id, ":") %>%
           map_chr(~ paste(.x[1:3], collapse = ":")))
introns_oa_sig <- introns_oa_sig %>%
  mutate(loc = str_split(phe_id, ":") %>%
           map_chr(~ paste(.x[1:3], collapse = ":")))


oa_psi_up <- introns_oa_sig %>% dplyr::filter(deltapsi_batch > 0) 
oa_psi_down <- introns_oa_sig %>% dplyr::filter(deltapsi_batch < 0) 
fnf_psi_up <- introns_fnf_sig %>% dplyr::filter(deltapsi_batch > 0)
fnf_psi_down <- introns_fnf_sig %>% dplyr::filter(deltapsi_batch < 0) 


fnf_psi_lookup <- setNames(introns_fnf_sig$deltapsi_batch, introns_fnf_sig$loc)
oa_psi_lookup <- setNames(introns_oa_sig$deltapsi_batch, introns_oa_sig$loc)

common_loc <- intersect(names(fnf_psi_lookup), names(oa_psi_lookup))
nonmatching_loc <- common_loc[sign(fnf_psi_lookup[common_loc]) != sign(oa_psi_lookup[common_loc])]


fnf_psi_same_direction_filtered <- introns_fnf_sig %>% dplyr::filter(!loc %in% nonmatching_loc)
oa_psi_same_direction_filtered <- introns_oa_sig %>% dplyr::filter(!loc %in% nonmatching_loc)

overlaps_intron_junctions_up <- list("PBS vs. FN-f" = fnf_psi_up$loc, "PBS vs. OA" = oa_psi_up$loc)
overlaps_intron_junctions_down <- list("PBS vs. FN-f" = fnf_psi_down$loc, "PBS vs. OA" = oa_psi_down$loc)
overlaps_intron_junctions_all <- list("PBS vs. FN-f" = fnf_psi_same_direction_filtered$loc, "PBS vs. OA" = oa_psi_same_direction_filtered$loc)

overlaps_intersect_match <- intersect(fnf_psi_same_direction_filtered$loc,oa_psi_same_direction_filtered$loc)
fnf_psi_subset_matches <- fnf_psi_same_direction_filtered %>% dplyr::filter(loc %in% overlaps_intersect_match)
oa_psi_subset_matches <- oa_psi_same_direction_filtered %>% dplyr::filter(loc %in% overlaps_intersect_match)
#ggvenn(
#  overlaps_intron_junctions_all,
#  fill_color =   c("#FDCDAC", "#CBD5E8"),
#  stroke_size = 0.5, set_name_size = 4,show_outside = "none",
#)


overlaps_venn <- ggvenn(overlaps_intron_junctions_all, 
       fill_color = c("#BFDDFF","#C8F0BF"), 
       stroke_color = NA, 
       auto_scale = TRUE, 
       show_percentage = FALSE, 
       text_size = 3, 
       set_name_size = 0) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"))
save(overlaps_venn, file = "output/results_plots/Supplementary_figures/Supp_fig2_overlaps_venn.rda")
fnf_psi_subset_matches_group <- fnf_psi_subset_matches %>% mutate(group=ifelse(deltapsi_batch >= 0.15,"high-confidence","low-confidence"))
highconf_fnf <- fnf_psi_subset_matches_group[which(fnf_psi_subset_matches_group$group == "high-confidence" ),]
oa_psi_subset_matches_group <- oa_psi_subset_matches %>% mutate(group=ifelse(deltapsi_batch >= 0.15,"high-confidence","low-confidence"))

highconf_oa <- oa_psi_subset_matches_group[which(oa_psi_subset_matches_group$group == "high-confidence" ),]

highconf_genes_oa <- oa_psi_subset_matches_group %>% dplyr::filter(genes %in% intersect(highconf_fnf$genes,highconf_oa$genes)) %>% 
  dplyr::filter(group == "high-confidence")  %>% pull(genes)
highconf_genes_fnf <- fnf_psi_subset_matches_group %>% dplyr::filter(genes %in% intersect(highconf_fnf$genes,highconf_oa$genes)) %>% 
  dplyr::filter(group == "high-confidence")  %>% pull(genes)


genes_with_commas_fnf <- highconf_genes_fnf %>%
  tibble(genes = .) %>%
  dplyr::filter(str_detect(genes, ",")) %>%
  pull(genes)

genes_with_commas <- highconf_genes_oa %>%
  tibble(genes = .) %>%
  dplyr::filter(str_detect(genes, ",")) %>%
  pull(genes)

chart_fnf <- table(highconf_fnf$group) |> as.data.frame() %>% 
  mutate(group="PBS vs. FN-f") %>%
  dplyr::rename(type = Var1)
chart_oa <- table(highconf_oa$group) |> as.data.frame() %>% 
  mutate(group="PBS vs. OA") %>%
  dplyr::rename(type = Var1)

#overlaps_both_test <- highconf_fnf %>% dplyr::filter(loc %in% intersect(highconf_fnf$loc,highconf_oa$loc)) |> nrow()
combined_chart <- rbind(chart_fnf, chart_oa)
combined_chart$group <- factor(combined_chart$group, levels = c("PBS vs. OA","PBS vs. FN-f"))

highConf_barPlot <- ggplot(combined_chart, aes(x = group, y = Freq, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge())+
  scale_fill_manual(values = c(lighten("#9edb91", 0.1),lighten("#649fe3",0.1))) +
  scale_y_continuous( name = "Total number of High confidence intron junction",
                      limits = c(0, 200), breaks = seq(0, 200, 50)) +
  annotate(geom = "richtext", label = paste("n=",combined_chart$Freq[1]) ,
           x = 2, y = 195, family = "Helvetica", size = 3, fill = NA, label.color = NA) +
  annotate(geom = "richtext", label = paste("n=",combined_chart$Freq[2]) ,
           x = 1, y = 90, family = "Helvetica", size = 3, fill = NA, label.color = NA) +
  coord_flip(clip = "off")+
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 0.25),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_blank(),
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
save(highConf_barPlot , file = "output/results_plots/Supplementary_figures/Supp_fig2_highConf_barPlot.rda")


#venndiagram high confidence
overlaps_genes_high_conf <- list("PBS vs. FN-f" = highconf_fnf$loc, "PBS vs. OA" = highconf_oa$loc)
overlaps_highconf_Venn <-  ggvenn(overlaps_genes_high_conf, 
                                  fill_color = c("#649fe3","#9edb91"), 
                                  stroke_color = NA, 
                                  auto_scale = TRUE, 
                                  show_percentage = FALSE, 
                                  text_size = 3, 
                                  set_name_size = 0) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"))
save(overlaps_highconf_Venn, file = "output/results_plots/Supplementary_figures/Supp_fig2_highconf_venn.rda")
ggsave(filename = "output/results_plots/Supplementary_figures/Supp_fig2_highconf_venn.pdf",
       plot = overlaps_highconf_Venn, width = 5, height = 5, units = "in")

#-------------------------------------------------------------------------------
# plotting Supp Fig2 DS PBS vs. OA
#-------------------------------------------------------------------------------

pdf(file = "output/results_plots/Supplementary_figures/supp_fig2_DSOA.pdf",   # The directory you want to save the file in
    width = 9.5, # The width of the plot in inches
    height = 7.2)
pageCreate(width = 9.5, height = 7.2, showGuides = TRUE)
plotText("A", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")



load(file="output/results_plots/Supplementary_figures/supp_fig2_heatmapOA.rda")
load(file="output/results_plots/Supplementary_figures/supp_fig2_heatmapOALegendGrob.rda")

plotGG(plot = heatmapGrob, x = 0.25, y = 0.2, height = 4.5, width = 5)

# Colorbar
plotGG(plot = heatmapLegendGrob, x = 0.35, y = 4.7,
       width = 4.325, height = 0.11)


# Colorbar title
plotText(label = "Relative Percent Spliced In (PSI)", fontfamily = "Helvetica",
         fontsize = 8, x = 2.45, y = 4.8, just = "top")

y=0.3
x=5.1
# Age Legend
for (i in 1:length(age_colors)) { 
  plotRect(x = unit(5, "in") + unit(3*i, "mm"),
           y = 0.3, width = unit(3, "mm"), 
           height = unit(1, "mm"), linecolor = NA, fill = age_colors[i],
           just = c("top" ,"left"))
}

ageText <- c("30","40","50","60","70","80+")

for(i in 1:length(ageText)){
  plotText(label = ageText[i], x = unit(5.05, "in") + unit(3.1*(i), "mm"),
           y =y+0.07 ,
           fontsize = 5, fontfamily = "Helvetica")
}


# Race Legend
Ancestry= c("AFR"="#C74A53", "AMR"="#EFB06E","EUR"= "#5BAD58", "SAS"="#177F97")

for (i in 1:4) {
  plotRect(x = unit(5.05, "in") + unit(4 * i, "mm"),
           y = y+0.07*2, width = unit(4, "mm"),
           height = unit(1, "mm"), linecolor = NA,
           fill = Ancestry[i], just = c("top", "left"))
}


raceText <- c("AFR","AMR","EUR","SAS")

for(i in 1:4){
  plotText(label = raceText[i], x = unit(5.13, "in") + unit(4*(i), "mm"),
           y = y+0.07*2+0.08,
           fontsize = 4, fontfamily = "Helvetica")
}

# Sex legend
plotRect(x = unit(5.1, "in") + unit(3*6, "mm"), 
         y = y+0.07*4, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = my_colour[['Sex']]['M'],
         just = "right")
plotRect(x = unit(5.1, "in") + unit(3*5, "mm"), 
         y = y+0.07*4, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = my_colour[['Sex']]['F'],
         just = "right")
plotText(label = "M", x = unit(x, "in") + unit(3*5.5, "mm"),
         y = y+0.07*5,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "F", x = unit(x, "in") + unit(3*4.5, "mm"),
         y =y+0.07*5,fontsize = 5, fontfamily = "Helvetica")

# Condition legend
plotRect(x = unit(x, "in") + unit(3*6, "mm"), 
         y =  y+0.07*6, width = unit(4, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = "#4A4A4A",
         just = "right")
plotRect(x = unit(x, "in") + unit(3*4.7, "mm"), 
         y =y+0.07*6, width = unit(4, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = "#B8B8B8",
         just = "right")
plotText(label = "OA", x = unit(x, "in") + unit(3*5.4, "mm"),
         y =  y+0.07*7,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "PBS", x = unit(x, "in") + unit(3*4, "mm"),
         y =y+0.07*7,
         fontsize = 5, fontfamily = "Helvetica")


#-------------------------------------------------------------------------------
plotText("B", x = 6, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")


dchart <- table(introns_oa_subset_sig$verdict) |> as.data.frame()
dchart$Percentage <- round((dchart$Freq / sum(dchart$Freq)),4)*100
dchart$ymax <- cumsum(dchart$Freq / sum(dchart$Freq))
dchart$ymin <- c(0, head(dchart$ymax, n=-1))

# Compute label position
dchart$labelPosition <- (dchart$ymax + dchart$ymin) / 2

# Compute a good label
lab.pos = cumsum(dchart$Percentage)-.5*(dchart$Percentage)
colnames(dchart)[1] <- "Verdict"
pieColors <- c("annotated"="lightgrey",
               "cryptic_fiveprime"="#96bce3" ,
               "cryptic_threeprime"="#3f85cc",
               "cryptic_unanchored"="#085099",
               "novel annotated pair"="#043363")
pieVerdict <-ggplot(data = dchart, aes(x = "", y = Percentage, fill = Verdict)) +
  geom_bar(stat = "identity") +
  coord_polar("y", direction = -1,start=55) +
  theme_void() +
  scale_fill_manual(values = pieColors) +
  theme(legend.position = "none", plot.background = element_rect(fill = 'transparent', color = NA))

plotGG(plot = pieVerdict, x = 6.25, y = 0.3, height = 1.3, width = 1.3)

sp_text <- c('annotated',"cryptic_5'","cryptic_3'","cryptic_unanchored","novel annotated pair")
plotText(label = sp_text[1], x = unit(7.8, "in"),
         y = 0.45,
         fontsize = 7, fontfamily = "Helvetica",
         just="left",
         fontcolor = "lightgrey",
         fontface = "bold")
plotText(label = "90.8%", x = unit(7.75, "in"),
         y = 0.45,
         fontsize = 7, fontfamily = "Helvetica",
         just="right",
         fontcolor = "lightgrey",
         fontface = "bold")


plotText(label = sp_text[2], x = unit(7.8, "in"),
         y = 0.45+0.15,
         fontsize = 7, fontfamily = "Helvetica",
         just="left",
         fontcolor = "#96bce3",
         fontface = "bold")
plotText(label = "3.0%", x = unit(7.75, "in"),
         y = 0.45+0.15,
         fontsize = 7, fontfamily = "Helvetica",
         just="right",
         fontcolor = "#96bce3",
         fontface = "bold")


plotText(label = sp_text[3], x = unit(7.8, "in"),
         y = 0.45+0.15*2,
         fontsize = 7, fontfamily = "Helvetica",
         just="left",
         fontcolor = "#3f85cc",
         fontface = "bold")
plotText(label = "4.5%", x = unit(7.75, "in"),
         y = 0.45+0.15*2,
         fontsize = 7, fontfamily = "Helvetica",
         just="right",
         fontcolor = "#3f85cc",
         fontface = "bold")


plotText(label =  sp_text[4], x = unit(7.8, "in"),
         y =  0.45+0.15*3,
         fontsize = 7, fontfamily = "Helvetica",
         just="left",
         fontcolor = "#085099",
         fontface = "bold")
plotText(label = "0.2%", x = unit(7.75, "in"),
         y =  0.45+0.15*3,
         fontsize = 7, fontfamily = "Helvetica",
         just="right",
         fontcolor = "#085099",
         fontface = "bold")

plotText(label =  sp_text[5], x = unit(7.8, "in"),
         y =  0.45+0.15*4,
         fontsize = 7, fontfamily = "Helvetica",
         just="left",
         fontcolor = "#043363",
         fontface = "bold")
plotText(label = "1.5%", x = unit(7.75, "in"),
         y =  0.45+0.15*4,
         fontsize = 7, fontfamily = "Helvetica",
         just="right",
         fontcolor = "#043363",
         fontface = "bold")

#------------------------------------------------------------------------------
#Plot C
#------------------------------------------------------------------------------
plotText("C", x = 6, y = 1.6, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

load(file="output/results_plots/Figure1_differntial_splicing/OApathway_barplots.rda")

plotGG(OApathway_barplots, x = 6.15, y = 1.8, width = 3.1, height = 2.3)

load(file="output/results_plots/Figure1_differntial_splicing/OAGO_barplots.rda")
plotGG(OAGO_barplots, x = 6.15, y = 4.15, width = 3.1, height =  2.3)

#------------------------------------------------------------------------------
#Plot D
#------------------------------------------------------------------------------

plotText("D", x = 0.25, y = 5, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
load("output/results_plots/Supplementary_figures/Supp_fig2_overlaps_venn.rda")
plotGG(overlaps_venn, x = -0.6, y = 5.2, width = 4.5, height = 2.0)

plotText(label = "PBS vs. FN-f", x = 1.1,
         y = 5.6,
         fontsize = 7, fontfamily = "Helvetica",
         just="right",
         fontcolor = "#F2BC40",fontface = "bold")

plotText(label = "PBS vs. OA", x = 2.3,
         y = 5.85,
         fontsize = 7, fontfamily = "Helvetica",
         just="left",
         fontcolor = "#2057A7",fontface = "bold")

plotText(label = "Differential intron junctions", x = 1.7,
         y = 5.25,
         fontsize = 7, fontfamily = "Helvetica",
         just="center",
         fontcolor = "black",fontface = "bold")

load("output/results_plots/Supplementary_figures/Supp_fig2_highconf_venn.rda")
plotGG(overlaps_highconf_Venn, x = 2, y =  5.2, width = 4.5, height =2.0)


plotSegments(
  x0 = 2.1, y0 = 6.2, x1 = 2.1, y1 = 5.5,
  default.units = "inches",
  lwd = 1, lty = 2,
  linecolor = "#505050"
)

plotSegments(
  x0 = 2.1, y0 = 5.5, x1 = 4.5, y1 = 5.5,
  default.units = "inches",
  lwd = 1, lty = 2,
  linecolor = "#505050"
)

plotSegments(
  x0 = 4.5, y0 = 5.5, x1 = 4.5, y1 = 5.6,
  default.units = "inches",
  lwd = 1, lty = 2,
  linecolor = "#505050"
)

plotText(label = "High confidence intron junctions", x = 4.5,
         y = 5.65,
         fontsize = 7, fontfamily = "Helvetica",
         just="center",
         fontcolor = "black",fontface = "bold")

dev.off()
