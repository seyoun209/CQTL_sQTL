# supp_fig4 for the sQTL
## Author: Seyoun Byun
## Date: 06.26.2024
#-------------------------------------------------------------------------------
library(plotgardener)
library(ggplot2)

pdf(file = "output/results_plots/Supplementary_figures/supp_fig4.pdf",   # The directory you want to save the file in
    width = 11.5, # The width of the plot in inches
    height = 3.75)

pageCreate(width = 11.5, height = 3.75, default.units = "inches",showGuides =FALSE)
plotText("B", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
load("output/results_plots/sqtl_plots/counts_dotplot.rda")
plotGG(counts_dotPlot, x = 0.3, y = 0.5, width = 4, height = 3)

plotText("C", x = 4.4, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
load("output/results_plots/sqtl_plots/venn_introns.rda")
load("output/results_plots/sqtl_plots/venn_sGenes.rda")
load("output/results_plots/sqtl_plots/enn_snp_introns_pairs.rda")

plotGG(sIntrons_venn, x = 3.5, y = 0, width = 4, height = 2.75)
plotGG(sGene_venn, x = 5.1, y = 0, width = 4, height = 2.75)
plotGG(snp_introns_pairs_venn, x = 4.3, y = 1.5, width = 4, height = 2.75)

plotText("D", x = 8, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
load( "output/results_plots/sqtl_plots/counts_rank_vennplot.rda")
plotGG(counts_rank_barplot,x = 8.3, y = 0.5, width = 3, height = 2.5)

dev.off()


