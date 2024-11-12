#supp_Figure4_version2_for the bottom part

# Function
create_venn_diagram <- function(total_PBS, total_FNF, overlap, only_PBS, only_FNF, max_count) {
  ggplot() +
    geom_circle(aes(x0 = -0.6, y0 = 0, r = sqrt(total_PBS/max_count)), 
                fill = "#BFDDFF", color = NA, alpha = 0.5) +
    geom_circle(aes(x0 = 0.6, y0 = 0, r = sqrt(total_FNF/max_count)), 
                fill = "#FFDDA2", color = NA, alpha = 0.5) +
    geom_text(aes(x = -0.8, y = 0, label = only_PBS), size = 3) +
    geom_text(aes(x = 0.8, y = 0, label = only_FNF), size = 3) +
    geom_text(aes(x = 0.1, y = 0, label = overlap), size = 3) +
    coord_fixed() +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    )
}

#-------------------------------------------------------------------------------
pdf(file = "./results_plots/sqtl_plots/fig4_supp_bottom_half_only_venn.pdf", width = 7.5, height = 2)
pageCreate(width = 7.5, height= 2 , default.units = "inches", showGuides = FALSE)


#-------------------------------------------------------------------------------

plotText("c", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

#splice intron-junction
plotText(label =  "sIntron-junctions", 
         x = 1.15,
         y= 0.4,
         fontsize = 8, fontfamily = "Helvetica",
         just="center", fontface = "bold",
         fontcolor = "black")

total_PBS <- length(pbs_sGene_introns)
total_FNF <- length(fnf_sGene_introns)
overlap <- length(intersect(pbs_sGene_introns, fnf_sGene_introns))
only_PBS <- total_PBS - overlap
only_FNF <- total_FNF - overlap
max_count <- max(total_PBS, total_FNF)

sIntrons_venn <- create_venn_diagram(total_PBS, total_FNF, overlap, only_PBS, only_FNF, max_count)

plotGG(plot =  sIntrons_venn, x =0, y = 0.25, height = 2, width = 2.25)

plotText(label =  "PBS", 
         x = 0.6,
         y= 0.6,
         fontsize = 8, fontfamily = "Helvetica",
         just="center",
         fontcolor = "#2057A7")

plotText(label =  "Fn-f", 
         x =1.5,
         y= 0.6,
         fontsize = 8, fontfamily = "Helvetica",
         just="center",
         fontcolor = "#F2BC40")

#cluster counts
total_PBS <- length(pbs_sGene_clusters)
total_FNF <- length(fnf_sGene_clusters)
overlap <- length(intersect(pbs_sGene_clusters, fnf_sGene_clusters))
only_PBS <- total_PBS - overlap
only_FNF <- total_FNF - overlap
max_count <- max(total_PBS, total_FNF)

sClusters_venn <- create_venn_diagram(total_PBS, total_FNF, overlap, only_PBS, only_FNF, max_count)
plotGG(plot =  sClusters_venn, x = 1.75, y =0.25 , height = 2, width = 2.25)
plotText(label = "sClusters", 
         x = 2.9,
         y=0.4,
         fontsize = 8, fontfamily = "Helvetica",
         just="center",fontface = "bold",
         fontcolor = "black")

plotText(label =  "PBS", 
         x = 2.35,
         y= 0.6,
         fontsize = 8, fontfamily = "Helvetica",
         just="center",
         fontcolor = "#2057A7")

plotText(label =  "Fn-f", 
         x =3.15,
         y= 0.6,
         fontsize = 8, fontfamily = "Helvetica",
         just="left",
         fontcolor = "#F2BC40")
#-------------------------------------------------------------------------------
plotText("d", x = 3.75, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

total_PBS <- sum(resqtl_intron_tibble$PBS)
total_FNF <- sum(resqtl_intron_tibble$FNF)
overlap <- sum(resqtl_intron_tibble$PBS & resqtl_intron_tibble$FNF)
only_PBS <- total_PBS - overlap
only_FNF <- total_FNF - overlap
max_count <- max(total_PBS, total_FNF)

resQTL_intron_venn <- create_venn_diagram(total_PBS, total_FNF, overlap, only_PBS, only_FNF, max_count)

plotGG(plot =  resQTL_intron_venn, x = 3.75, y =0.25 , height = 2, width = 2.25)
plotText(label = "re-sIntron-junctions", 
         x = 5,
         y=0.4,
         fontsize = 8, fontfamily = "Helvetica",
         just="center",fontface = "bold",
         fontcolor = "black")

plotText(label =  "PBS", 
         x = 4.35,
         y= 0.6,
         fontsize = 8, fontfamily = "Helvetica",
         just="center",
         fontcolor = "#2057A7")

plotText(label =  "Fn-f", 
         x =5.25,
         y= 0.6,
         fontsize = 8, fontfamily = "Helvetica",
         just="left",
         fontcolor = "#F2BC40")
# Cluster for response----------------------------------------------------------
total_PBS <- sum(resqtl_clusters$PBS)
total_FNF <- sum(resqtl_clusters$FNF)
overlap <- sum(resqtl_clusters$PBS & resqtl_clusters$FNF)
only_PBS <- total_PBS - overlap
only_FNF <- total_FNF - overlap
max_count <- max(total_PBS, total_FNF)

resQTL_cluster_venn <- create_venn_diagram(total_PBS, total_FNF, overlap, only_PBS, only_FNF, max_count)

plotGG(plot =  resQTL_cluster_venn, x = 5.5, y =0.25 , height = 2, width = 2.25)
plotText(label = "re-sClusters", 
         x = 6.75,
         y=0.4,
         fontsize = 8, fontfamily = "Helvetica",
         just="center",fontface = "bold",
         fontcolor = "black")

plotText(label =  "PBS", 
         x = 6.1,
         y= 0.6,
         fontsize = 8, fontfamily = "Helvetica",
         just="center",
         fontcolor = "#2057A7")

plotText(label =  "Fn-f", 
         x =6.9,
         y= 0.6,
         fontsize = 8, fontfamily = "Helvetica",
         just="left",
         fontcolor = "#F2BC40")

dev.off()
