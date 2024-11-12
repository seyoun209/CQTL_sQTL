setwd("/work/users/s/e/seyoun/CQTL_sQTL")

#Venndiagrm-redo (Oct 19th 2024)
load("output/clu_fnf/introns_fnf_joinAll")
introns_fnf_sig <- introns_fnf_pval_include %>% dplyr::filter(p.adjust < 0.05)
introns_fnf_subset_sig <- introns_fnf_sig[abs(introns_fnf_sig$deltapsi_batch) > 0.15,]
introns_fnf_all_sig <- introns_fnf_sig %>% dplyr::filter(abs(deltapsi_batch) > 0.15)
load("./output/clu_oa/introns_oa_joinAll")
introns_oa_sig <- introns_oa_pval_include %>% dplyr::filter(p.adjust < 0.05)
introns_oa_subset_sig <- introns_oa_sig[abs(introns_oa_sig$deltapsi_batch) > 0.15,]

overlaps_genes_all <- list("FN-f/PBS" = introns_fnf_subset_sig$ensemblID, "OA/PBS" = introns_oa_subset_sig$ensemblID)
overlaps_venn <- ggvenn(overlaps_genes_all,
                        fill_color = c("#BFDDFF","#C8F0BF"),
                        stroke_color = NA,
                        auto_scale = TRUE,
                        show_percentage = FALSE,
                        set_name_size = 0,
                        text_size = 2) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"))
save(overlaps_venn, file = "output/results_plots/Supplementary_figures/overall_venn_figure1d.rda")
ggsave(filename = "output/results_plots/Supplementary_figures/overall_venn_figure1d.pdf", 
       plot =overlaps_venn,width = 4, height = 4, units = "in")

#Figure1 new version of the generation
pdf(file = "output/results_plots/Figure1_differntial_splicing/figure1_v3.pdf",   # The directory you want to save the file in
    width = 11.8, # The width of the plot in inches
    height = 8)

pageCreate(width = 11.8, height = 8, showGuides = FALSE)
plotText("a", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
load(file="output/results_plots/Figure1_differntial_splicing/heatmapGrob.rda")
load(file="output/results_plots/Figure1_differntial_splicing/heatmapLegendGrob.rda")

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
age_labels <- c("30-39", "40-49", "50-59", "60-69", "70-79", "80+")
age_palette <- c("#ffcba4","#cca283","#997a62","#806652","#665142","#332921")
# Create a named vector of colors for the age groups
age_colors <- setNames(age_palette, age_labels)


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

my_colour = list(
  Condition  = c("CTL" = "#9FCCE4", "FNF" = "#FAB394"),
  Sex  = c( "M" =  "#0075B0", "F" = "#F8B7CD"),
  Predicted_Ancestry= c("AMR"="#F5BC9F", "EUR"="#FAF1D2" ,"AFR"="#86CBB7","SAS"="#6B7EA4"),
  Age_range = age_colors)


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
plotText(label = "FNF", x = unit(x, "in") + unit(3*5.4, "mm"),
         y =  y+0.07*7,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "PBS", x = unit(x, "in") + unit(3*4, "mm"),
         y =y+0.07*7,
         fontsize = 5, fontfamily = "Helvetica")


#------------------------------------------------------------------------------
#Plot B

plotText("b", x = 0.1, y = 5.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
for (i in 1:5) {
  if(i == 4){
    plotSegments(
      x0 = unit(0.6, "in"),
      y0 = 5.4+(0.15 * 4),
      x1 = unit(1.6, "in")+unit(7, "mm"),
      y1 = 5.4+(0.15 * 4),
      default.units = "inches",
      lwd = 1.5, lty = 1,"darkgrey",
      just ="left"
    )
    plotRect(x = unit(0.6, "in") + unit(7, "mm"),
             y = 5.4+(0.15 * i), width = unit(3, "mm"),
             height = unit(2, "mm"), linecolor = NA, fill = "#2473b5",
             just = "left")
    plotRect(x = unit(1.6, "in"),
             y = 5.4+(0.15 * i), width = unit(3, "mm"),
             height = unit(2, "mm"), linecolor = NA, fill = "#2473b5",
             just = "right")
  } else if(i == 5){
    plotRect(x = unit(0.6, "in"),
             y = 5.4+(0.15 * i), width = unit(7, "mm"),
             height = unit(2, "mm"), linecolor = NA, fill = "#434e57",
             just = "left")
    plotRect(x = unit(1.6, "in"),
             y = 5.4+(0.15 * i), width = unit(7, "mm"),
             height = unit(2, "mm"), linecolor = NA, fill = "#96bce3",
             just = "left")
    plotSegments(
      x0 = unit(0.6, "in") + unit(7, "mm"),
      y0 = 5.4+(0.15 * i),
      x1 = 1.6,
      y1 = 5.4+(0.15 * i),
      default.units = "inches",
      lwd = 1.5, lty = 1,"darkgrey",
      just ="left"
    )
  } else {
    plotRect(x = unit(0.6, "in"),
             y = 5.4+(0.15 * i), width = unit(7, "mm"),
             height = unit(2, "mm"), linecolor = NA, fill = "#434e57",
             just = "left")
    plotRect(x = unit(1.6, "in"),
             y = 5.4+(0.15 * i), width = unit(7, "mm"),
             height = unit(2, "mm"), linecolor = NA, fill = "#434e57",
             just = "left")
    plotSegments(
      x0 = unit(0.6, "in") + unit(7, "mm"),
      y0 = 5.4+(0.15 * i),
      x1 = 1.6,
      y1 = 5.4+ (0.15 * i),
      default.units = "inches",
      lwd = 1.5, lty = 1,"darkgrey",
      just ="left"
    )
  }
}

# Updated text annotations
sp_text <- c('annotated',"cryptic_5'","cryptic_3'","cryptic_unanchored","novel annotated pair")
for (i in 1:5) {
  plotText(label = sp_text[i], x = unit(1.9, "in"),
           y = 5.4+(0.15 * i),
           fontsize = 7, fontfamily = "Helvetica",
           just="left")
}

# Annotation legend
anno_text <- c(': annotated SS',": unannotated SS",": unpaired SS")
annocolor <- c("#434e57" ,"#2473b5","#96bce3")
for (i in 1:3) {
  plotRect(x = unit(-0.1+(0.75*i),"in") ,  # Adjust x here to align with panel b
           y = 6.3, width = unit(1.5, "mm"),
           height = unit(1.5, "mm"), linecolor = NA, fill = annocolor[i],
           just = "left")
}

# Adjust legend labels
plotText(label = 'annotated SS', 
         x = 0.75,  # Adjust x to align with rectangles
         y = 6.3,
         fontsize = 5, fontfamily = "Helvetica",
         just = "left", fontcolor = "#434e57")

plotText(label = "unannotated SS", 
         x = 0.75+ 0.77,  # Adjust x for second label
         y = 6.3,
         fontsize = 5, fontfamily = "Helvetica",
         just = "left", fontcolor = "#2473b5")

plotText(label = "unpaired SS", 
         x = 0.71 + 0.77 * 2,  # Adjust x for third label
         y = 6.3,
         fontsize = 5, fontfamily = "Helvetica",
         just = "left", fontcolor = "#96bce3")


## dchart
dchart <- table(introns_fnf_all_sig$verdict) |> as.data.frame()
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

plotGG(plot = pieVerdict, x = 0.5, y = 6.6, height = 1.3, width = 1.3)

plotText(label = sp_text[1], x = unit(2.2, "in"),
         y = 6.8,
         fontsize = 7, fontfamily = "Helvetica",
         just="left",
         fontcolor = "grey50",
         fontface = "bold")
plotText(label = "90.86%", x = unit(2.15, "in"),
         y = 6.8,
         fontsize = 7, fontfamily = "Helvetica",
         just="right",
         fontcolor = "grey50",
         fontface = "bold")


plotText(label = sp_text[2], x = unit(2.2, "in"),
         y = 6.8+0.15,
         fontsize = 7, fontfamily = "Helvetica",
         just="left",
         fontcolor = "#96bce3",
         fontface = "bold")
plotText(label = "6.88%", x = unit(2.15, "in"),
         y = 6.8+0.15,
         fontsize = 7, fontfamily = "Helvetica",
         just="right",
         fontcolor = "#96bce3",
         fontface = "bold")


plotText(label = sp_text[3], x = unit(2.2, "in"),
         y = 6.8+0.15*2,
         fontsize = 7, fontfamily = "Helvetica",
         just="left",
         fontcolor = "#3f85cc",
         fontface = "bold")
plotText(label = "1.03%", x = unit(2.15, "in"),
         y = 6.8+0.15*2,
         fontsize = 7, fontfamily = "Helvetica",
         just="right",
         fontcolor = "#3f85cc",
         fontface = "bold")


plotText(label =  sp_text[4], x = unit(2.2, "in"),
         y =  6.8+0.15*3,
         fontsize = 7, fontfamily = "Helvetica",
         just="left",
         fontcolor = "#085099",
         fontface = "bold")
plotText(label = "0.62%", x = unit(2.15, "in"),
         y =  6.8+0.15*3,
         fontsize = 7, fontfamily = "Helvetica",
         just="right",
         fontcolor = "#085099",
         fontface = "bold")

plotText(label =  sp_text[5], x = unit(2.2, "in"),
         y =  6.8+0.15*4,
         fontsize = 7, fontfamily = "Helvetica",
         just="left",
         fontcolor = "#043363",
         fontface = "bold")
plotText(label = "0.62%", x = unit(2.15, "in"),
         y =  6.8+0.15*4,
         fontsize = 7, fontfamily = "Helvetica",
         just="right",
         fontcolor = "#043363",
         fontface = "bold")

#------------------------------------------------------------------------------
#Plot C
plotText("c", x = 6, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
library(png)

PPI_test <- readPNG("output/results_plots/Figure1_differntial_splicing/CTL_FNF_deltapsi15_added_name.png")
plotRaster(image= PPI_test,
           x=6, y=0.25,width=5.7,height=4.5,
           just = c("top","left"))

#------------------------------------------------------------------------------
#Plot d
plotText("d", x = 3.5, y = 5.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
load(file = "output/results_plots/Supplementary_figures/overall_venn_figure1d.rda")
plotGG(overlaps_venn, x = 2.8, y = 5.7, width = 4.5, height = 1.85)

plotText(label = "Alternative splice intron junctions\n(padj < 0.05 & |ΔPSI| > 0.15)", x = 5.05,
         y = 5.5,
         fontsize = 7, fontfamily = "Helvetica",
         just="center",
         fontcolor = "black",fontface = "bold")

plotText(label = "FN-f/PBS", x = 4.2,
         y = 6.1,
         fontsize = 7, fontfamily = "Helvetica",
         just="right",
         fontcolor = "#2057A7",fontface = "bold")

plotText(label = "OA/PBS", x = 5.85,
         y = 6.2,
         fontsize = 7, fontfamily = "Helvetica",
         just="left",
         fontcolor = "#0E8F4C",fontface = "bold")

#Listing genes

genes_overlaps <- intersect(introns_fnf_subset_sig$gene, introns_oa_subset_sig$gene)
ensg_gene <- genes_overlaps[grep("^ENSG", genes_overlaps)]
non_ensg_genes <- genes_overlaps[!genes_overlaps %in% ensg_gene]

# Sort non-ENSG genes alphabetically
sorted_genes <- sort(non_ensg_genes)

# Append the ENSG gene at the end
sorted_genes <- c(sorted_genes, ensg_gene)


x_pos <- 6.3
y_start <- 5.25  # Adjusted start point for writing

# Number of genes to display per column
genes_per_column <- 19

# Vertical spacing between genes
y_spacing <- 0.15

# Horizontal spacing between columns
x_spacing <- 0.65

# Loop through all the genes and plot them
for (i in 1:length(sorted_genes)) {
  # Determine the current column
  col_num <- ceiling(i / genes_per_column)

  # Determine the current row in the column
  row_num <- (i - 1) %% genes_per_column + 1

  # Adjust x and y positions based on the column and row
  x_current <- x_pos + (col_num - 1) * x_spacing
  y_current <- y_start + (row_num - 1) * y_spacing  # Adjusted to start from y_start and go downwards

  # Plot the gene text at the calculated position
  plotText(label = sorted_genes[i],
           x = x_current,
           y = y_current,
           fontsize = 7, fontfamily = "Helvetica",
           just = "left",
           fontcolor = "black")
}


# load("output/results_plots/Supplementary_figures/Supp_fig2_overlaps_venn.rda")
# plotGG(overlaps_venn, x = 3, y = 5.2, width = 4.5, height = 1.85)
# 
# plotText(label = "FN-f/PBS", x = 4.4,
#          y = 5.6,
#          fontsize = 7, fontfamily = "Helvetica",
#          just="right",
#          fontcolor = "#2057A7",fontface = "bold")
# 
# plotText(label = "OA/PBS", x = 5.8,
#          y = 5.7,
#          fontsize = 7, fontfamily = "Helvetica",
#          just="left",
#          fontcolor = "#0E8F4C",fontface = "bold")
# 
# 
# plotText(label = "Alternative splice intron junctions (padj < 0.05)", x = 5.25,
#          y = 5.3,
#          fontsize = 7, fontfamily = "Helvetica",
#          just="center",
#          fontcolor = "black",fontface = "bold")
# 
# load("output/results_plots/Supplementary_figures/Supp_fig2_highconf_venn.rda")
# plotGG(overlaps_highconf_Venn, x = 6, y =  5.8, width = 2, height =0.9)
# 
# 
# 
# plotSegments(
#   x0 = 5.65, y0 = 6.08, x1 = 5.65, y1 = 5.55,
#   default.units = "inches",
#   lwd = 1, lty = 2,
#   linecolor = "#505050"
# )
# 
# plotSegments(
#   x0 = 5.65, y0 = 5.55, x1 = 7.1, y1 = 5.55,
#   default.units = "inches",
#   lwd = 1, lty = 2,
#   linecolor = "#505050"
# )
# 
# plotSegments(
#   x0 = 7.1, y0 = 5.55, x1 = 7.1, y1 = 5.7,
#   default.units = "inches",
#   lwd = 1, lty = 2,
#   linecolor = "#505050"
# )
# 
# 
# plotText(label = "High confidence splice intron junctions", 
#          x = 7.1,
#          y = 5.72,
#          fontsize = 7, fontfamily = "Helvetica",
#          just = "center",
#          fontcolor = "black", fontface = "bold")
# 
# # Second line of the text (slightly below the first line)
# plotText(label = "(padj < 0.05 , |ΔPSI| > 0.15)", 
#          x = 7.1,
#          y =5.82,  # Adjust y to position the second line right below the first
#          fontsize = 7, fontfamily = "Helvetica",
#          just = "center",
#          fontcolor = "black")
# 
# plotText(label = "FN-f/PBS", x = 6.45,
#          y = 6.0,
#          fontsize = 7, fontfamily = "Helvetica",
#          just="center",
#          fontcolor = "#2057A7",fontface = "bold")
# 
# plotText(label = "OA/PBS", x = 7.5,
#          y = 6.2,
#          fontsize = 7, fontfamily = "Helvetica",
#          just="left",
#          fontcolor = "#0E8F4C",fontface = "bold")
# 
# 
# #Listing genes
# 
# genes_overlaps <- intersect(highconf_fnf$gene, highconf_oa$gene)
# ensg_gene <- genes_overlaps[grep("^ENSG", genes_overlaps)]
# non_ensg_genes <- genes_overlaps[!genes_overlaps %in% ensg_gene]
# 
# # Sort non-ENSG genes alphabetically
# sorted_genes <- sort(non_ensg_genes)
# 
# # Append the ENSG gene at the end
# sorted_genes <- c(sorted_genes, ensg_gene)
# 
# 
# x_pos <- 4
# y_start <- 7.15  # Adjusted start point for writing
# 
# # Number of genes to display per column
# genes_per_column <- 6
# 
# # Vertical spacing between genes
# y_spacing <- 0.15
# 
# # Horizontal spacing between columns
# x_spacing <- 0.65
# 
# # Loop through all the genes and plot them
# for (i in 1:length(sorted_genes)) {
#   # Determine the current column
#   col_num <- ceiling(i / genes_per_column)
#   
#   # Determine the current row in the column
#   row_num <- (i - 1) %% genes_per_column + 1
#   
#   # Adjust x and y positions based on the column and row
#   x_current <- x_pos + (col_num - 1) * x_spacing
#   y_current <- y_start + (row_num - 1) * y_spacing  # Adjusted to start from y_start and go downwards
#   
#   # Plot the gene text at the calculated position
#   plotText(label = sorted_genes[i], 
#            x = x_current, 
#            y = y_current,  
#            fontsize = 7, fontfamily = "Helvetica",
#            just = "left",
#            fontcolor = "black")
# }

#-------------------------------------------------------------------------------

plotText("e", x = 9, y = 5.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")


load(file="output/results_plots/Supplementary_figures/Fig1E_OA_boxplot.rda")
plotGG(oa_boxplots, x = 9.2, y = 5.5, width = 2.8, height =2.4 )

dev.off()