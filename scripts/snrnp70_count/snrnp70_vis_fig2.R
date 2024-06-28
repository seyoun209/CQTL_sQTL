# Figure2_SNRNP70
## Author: Seyoun Byun
## Date: 06.22.2024
#-------------------------------------------------------------------------------
setwd("/work/users/s/e/seyoun/CQTL_sQTL/output")
library(dplyr)
library(plotgardener)
library(RColorBrewer)
source("../scripts/utils/utils.R")
source("../scripts/utils/grid_boxplot.R")

yl_gn_bu <- brewer.pal(n = 9, name = "YlGnBu")
"#FFFFD9" "#EDF8B1" "#C7E9B4" "#7FCDBB" "#41B6C4" "#1D91C0" "#225EA8" "#253494" "#081D58"

#Signals pbs fnf oa

dir_merged <- "/work/users/s/e/seyoun/CQTL_sQTL/output/signals/merged_norm/"
ctl_signal <- paste0(dir_merged,"CTL_norm.bw")
fnf_signal <- paste0(dir_merged,"FNF_norm.bw")
oa_signal <- paste0(dir_merged,"OA_norm.bw")

dir_edited <- "/work/users/s/e/seyoun/crispr/02.test_seq/signals/norm_merged/"
wd_signal <- paste0(dir_merged,"WD.bw")
kd_signal <- paste0(dir_merged,"KD.bw")
#-------------------------------------------------------------------------------
#Making a boxplot for Fig2D
#-------------------------------------------------------------------------------
#Comparing with PBS vs. fnf
library(org.Hs.eg.db)
degenes_fnf_shrink <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/quant/de_genes_results.csv") |>
  dplyr::mutate(gene_id = gsub("\\..*$", "", gene_id))
annots <- AnnotationDbi::select(org.Hs.eg.db, degenes_fnf_shrink$gene_id,
                 columns=c("SYMBOL"), keytype="ENSEMBL")
annots <- dplyr::rename(annots, gene_id = ENSEMBL)
annots_unique <- annots %>% distinct(gene_id, .keep_all = TRUE)
degenes_fnf_shrink_symbol <- left_join(degenes_fnf_shrink, annots_unique, by = "gene_id")
degenes_fnf_shrink_symbol <- degenes_fnf_shrink_symbol %>% dplyr::select("gene_id","SYMBOL","log2FoldChange","padj")

#edited 
edited_dir <- "/work/users/s/e/seyoun/crispr/02.test_seq/"
load(paste0(edited_dir,"condition_de/differential_expression_dds.rda"))
de_gene_edited  <- fread(paste0(edited_dir,"condition_de/de_genes_results_with_symbol.csv"))
down_l2fc2_padj01_edited_degenes <- fread(paste0(edited_dir,"condition_de/downsig_deGenes_pval01_l2fc2.csv")) |>
  mutate(dir = ifelse(log2FoldChange > 0, "up", "down"))  |> 
  dplyr::rename(log2FoldChange_edited = log2FoldChange, dir_edited = dir)

up_l2fc2_padj01_edited_degenes <- fread(paste0(edited_dir,"condition_de/upsig_deGenes_pval01_l2fc2.csv")) |>
  mutate(dir = ifelse(log2FoldChange > 0, "up", "down"))  |> 
  dplyr::rename(log2FoldChange_edited = log2FoldChange, dir_edited = dir)

up_fnf_edited_subset <- degenes_fnf_shrink_symbol %>% dplyr::filter(gene_id %in% up_l2fc2_padj01_edited_degenes$gene_id)
down_fnf_edited_subset <- degenes_fnf_shrink_symbol %>% dplyr::filter(gene_id %in% down_l2fc2_padj01_edited_degenes$gene_id)


fnf_all_edited_subset <- bind_rows(up_fnf_edited_subset |> mutate(group = "Up in Het-KD"),
                               down_fnf_edited_subset |> mutate(group = "Down in Het-KD")) |>
  mutate(group = factor(group, levels = c("Up in Het-KD", "Down in Het-KD")))


up_test_kd <- wilcox.test(x = up_fnf_edited_subset$log2FoldChange, mu=0,
                          alternative = "greater")

down_test_kd <- wilcox.test(x = down_fnf_edited_subset$log2FoldChange,mu=0,
                            alternative = "less")



fnf_all_edited_subset <- fnf_all_edited_subset |>
  mutate(highlight = case_when(SYMBOL %in% c("MMP12","CXCL8","MMP13","IL1A","IL11","CXCL2") ~ "up",
                               SYMBOL %in% c("SEMA3E","PEG3","FGL2","HLF","CSRNP3") ~ "down")) |>
  mutate(highlight = factor(highlight, levels = c("up", "down")))  %>%
  dplyr::select(-gene_id) %>%
  dplyr::rename(symbol = SYMBOL)

library(gridtext)
library(colorspace)
library(grid)
edited_boxplot_plot(data = fnf_all_edited_subset, up_wilcox_test = up_test_kd,
               down_wilcox_test = down_test_kd,
               x = unit(6, "native"),
               y = unit(5, "native"), width = 3.5, height = 2.5)



#-------------------------------------------------------------------------------
# Figure2 SNRNP70
pdf(file = "results_plots/Figure2_SNRNP70/figure2.pdf",   # The directory you want to save the file in
    width = 9.75, # The width of the plot in inches
    height = 8.0)
pageCreate(width = 9.75, height = 8.0, default.units = "inches",showGuides = FALSE)
plotText("A", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
chr.nm = "chr19"
minpos = 49085450
maxpos = 49108604
p.s <-  pgParams(chrom = chr.nm, chromstart = minpos, chromend = maxpos,
                 assembly = "hg38",
                 resolution = 10000,
                 width =  4.0,
                 x = 0.5,
                 y = 0.35)

trans_highlights_nmd <- data.frame(
  transcript = c("ENST00000598441.6","ENST00000401730.5","ENST00000221448.9","ENST00000595231.5","ENST00000601065.5","ENST00000599687.1","ENST00000597936.5","ENST00000544278.2","ENST00000601411.1","ENST00000598984.5"),
  color = c('grey50', "grey50",'white',"white","white","white","white","white","white","white")
)

nmdplot <- plotTranscripts(
  params = p.s,
  assembly = "hg38",
  height = 1.5,
  y=0.55,
  just = c("left", "top"), default.units = "inches",
  transcriptHighlights=trans_highlights_nmd,
  labels = NULL,spaceHeight =0
)



height=1
snrnp70_signals <- plotMultiSignal(data = list(ctl_signal,
                                               fnf_signal,
                                               oa_signal),
                                   params = p.s,
                                   y = 0.35, height = 1.11, linecolor = c(yl_gn_bu[3],yl_gn_bu[6],yl_gn_bu[7]), 
                                   fill = c(yl_gn_bu[3],yl_gn_bu[6],yl_gn_bu[7]),
                                   #label = c("PBS", "FNF","OA"),
                                   default.units = "inches",
                                   gapdistance = 0.02)
height=0.35
plotText(label = "PBS", fontcolor = yl_gn_bu[3],
         fontsize = 7, x = 0.5, y = 0.4, just = "left", fontfamily = "Helvetica")
plotText(label = "FN-f", fontcolor = yl_gn_bu[6],
         fontsize = 7, x = 0.5, y = 0.42+height, just = "left", fontfamily = "Helvetica")
plotText(label = "OA", fontcolor = yl_gn_bu[7],
         fontsize = 7, x = 0.5, y = 0.44+height*2, just = "left", fontfamily = "Helvetica")



plotText(label= "ex1", x=0.55, y=2.1, fontcolor = "black",fontsize = 8,just="center")
plotText(label= "ex2", x=0.73, y=2.1, fontcolor = "black",fontsize = 8,just="center")
plotText(label= "ex3", x=1.33, y=2.1, fontcolor = "black",fontsize = 8,just="right")
plotText(label= "ex4", x=1.35, y=2.1, fontcolor = "black",fontsize = 8,just="left")
plotText(label= "ex5", x=2.76, y=2.1, fontcolor = "black",fontsize = 8,just="right")
plotText(label= "ex6", x=2.77, y=2.1, fontcolor = "black",fontsize = 8,just="left")
plotText(label= "ex7", x=3.265, y=2.1, fontcolor = "black",fontsize = 8,just="center")
plotText(label= "altex8", x=3.4, y=2.1, fontcolor = "#E07653",fontsize = 8,just="left")
plotText(label= "ex8", x=3.8, y=2.1, fontcolor = "black",fontsize = 8,just="center")
plotText(label= "ex9", x=4.36, y=2.1, fontcolor = "black",fontsize = 8,just="right")
plotText(label= "ex10", x=4.367, y=2.1, fontcolor = "black",fontsize = 8,just="left")
annoGenomeLabel(plot=p.s,y=2.2,x=0.5, scale="bp",fontsize = 8)

# adding the Count plot 


load(file="results_plots/Figure2_SNRNP70/snrnp70_combinedplot.rda")
#load(file="results_plots/Figure2_SNRNP70/snrnp70_combinedplot_r43.rda")
plotGG(combined_plot, x = 0.35, y = 2.4, width = 4.35, height = 3)

#-------------------------------------------------------------------------------
#B Wildtype vs Knockdown alt exon 8 
#-------------------------------------------------------------------------------
plotText("B", x = 5, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

trans_highlights_nmd <- data.frame(
  transcript = c("ENST00000598441.6","ENST00000401730.5","ENST00000221448.9","ENST00000595231.5","ENST00000601065.5","ENST00000599687.1","ENST00000597936.5","ENST00000544278.2","ENST00000601411.1","ENST00000598984.5"),
  color = c('white', "grey50",'white',"white","white","white","white","white","white","white")
)

ps_edited <-  pgParams(chrom = chr.nm, chromstart = minpos, chromend = maxpos,
                 assembly = "hg38",
                 resolution = 10000,
                 width =  4,
                 x = 5.5,
                 y = 0.35)

targed_enst <- plotTranscripts(
  params = ps_edited,
  assembly = "hg38",
  height = 1.2,
  y=0,
  just = c("left", "top"), default.units = "inches",
  transcriptHighlights=trans_highlights_nmd,
  labels = NULL,spaceHeight =0,
  boxHeight= unit(0.5, "cm")
)

plotText(label= "ex1", x=5.47, y=1.15, fontcolor = "black",fontsize = 8,just="center")
plotText(label= "ex2", x=5.62, y=1.15, fontcolor = "black",fontsize = 8,just="left")
plotText(label= "ex3", x=6.37, y=1.15, fontcolor = "black",fontsize = 8,just="right")
plotText(label= "ex4", x=6.375, y=1.15, fontcolor = "black",fontsize = 8,just="left")
plotText(label= "ex5", x=7.75, y=1.15, fontcolor = "black",fontsize = 8,just="right")
plotText(label= "ex6", x=7.755, y=1.15, fontcolor = "black",fontsize = 8,just="left")
plotText(label= "ex7", x=8.25, y=1.15, fontcolor = "black",fontsize = 8,just="center")
plotText(label= "altex8", x=8.395, y=1.15, fontcolor = "#E07653",fontsize = 8,just="left")
plotText(label= "ex8", x=8.8, y=1.15, fontcolor = "black",fontsize =  8,just="center")
plotText(label= "ex9", x=9.4, y=1.15, fontcolor = "black",fontsize =  8,just="right")
plotText(label= "ex10", x=9.405, y=1.15, fontcolor = "black",fontsize =  8,just="left")
annoGenomeLabel(plot=ps_edited,y=1.4,x=5.5, scale="bp",fontsize = 8)

#load(file="/work/users/s/e/seyoun/crispr/02.test_seq/plots/combined_plot_snrnp70_edited.rda")
load(file="/work/users/s/e/seyoun/crispr/02.test_seq/plots/combined_plot_snrnp70_edited_r44.rda")

plotGG(combined_plot, x = 5.35, y = 1.5, width = 4.35, height = 3)


#-------------------------------------------------------------------------------
##plot C - GO, KEGG
#-------------------------------------------------------------------------------
plotText("C", x = 0.1, y = 5.35, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

load(file="/work/users/s/e/seyoun/crispr/02.test_seq/plots/GO_barplots.rda")
load(file="/work/users/s/e/seyoun/crispr/02.test_seq/plots/pwathway_barplots.rda")

plotGG(pwathway_barplots, x = 0.315, y = 5.45, width = 2.5, height = 2.5)
plotGG(GO_barplots, x = 2.6, y = 5.45, width =4.0 , height = 2.5)


#-------------------------------------------------------------------------------
##plot D - boxplot
#-------------------------------------------------------------------------------
plotText("D", x = 6.5, y = 4.65, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

edited_boxplot_fnf_subset <- edited_boxplot_plot(data = fnf_all_edited_subset, up_wilcox_test = up_test_kd,
                    down_wilcox_test = down_test_kd,
                    x = unit(6.9, "native"),
                    y = unit(5, "native"), width = 4.4, height = 2.7)

dev.off()
