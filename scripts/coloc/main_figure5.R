setwd("/work/users/s/e/seyoun/CQTL_sQTL/")
library(plotgardener)
library(RColorBrewer)
library(grid)

response_pbs_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_pbs_resInv.rds")
response_fnf_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_fnf_resInv.rds")

pbs_norm_qtl_pc5_list <- readRDS("output/nominal_1mb/nominal_pbs/pbs_norm_qtl_pc5_list.rds")
fnf_norm_qtl_pc4_list <- readRDS("output/nominal_1mb/nominal_fnf/fnf_norm_qtl_pc4_list.rds")

load("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/freq/maf_id_add.rds") #maf_id_add  is the name
maf_subset_rsID <- maf_id_add |> dplyr::select(-c("minor_allele","MAF"))


#Plotting for per all the coloc FNF genes
coloc_results_fnf <- fread("output/coloc/coloc_result_table_fnf_prob_calc_.txt")
coloc_results_pbs <- fread("output/coloc/coloc_result_table_pbs_prob_calc_.txt")
sig_fnf_coloc <- coloc_results_fnf |> dplyr::filter(FNF_coloc_H4 >= 0.7) |> distinct(phe_id,sQTL_rsID, .keep_all = TRUE)
sig_pbs_coloc <- coloc_results_pbs |> dplyr::filter(PBS_coloc_H4 >= 0.7) |> distinct(phe_id,sQTL_rsID, .keep_all = TRUE)

#Locus zoom plot 

plot_gwas_qtl_locusZoom <- function(gene_name, primary_dataset, x_start, y_start, width, height, zoom_range = 200000) {
  # Validate input
  if (!primary_dataset %in% c("PBS", "FNF")) {
    stop("Primary dataset must be either 'PBS' or 'FNF'")
  }
  
  # Select the appropriate data-set
  highConf_resQtL <- if(primary_dataset == "PBS") response_pbs_results else response_fnf_results
  coloc_results <- if(primary_dataset == "PBS") sig_pbs_coloc else sig_fnf_coloc
  
  # Filter data for the specified gene
  test_boxplotInfo <- highConf_resQtL |> dplyr::filter(SYMBOL %in% gene_name ) 
  if (nrow(test_boxplotInfo) != 1) {
    max_coloc_phe_id <- coloc_results |> dplyr::filter(gene %in% gene_name) |> 
      slice_max(!!sym(paste0(primary_dataset, "_coloc_H4"))) |> 
      pull(phe_id)
    test_boxplotInfo <- test_boxplotInfo |> dplyr::filter(phe_id == max_coloc_phe_id)
  }
  
  # Extract necessary information
  chrom <- test_boxplotInfo$phe_chr
  variantID <- unique(test_boxplotInfo$var_id)
  rsID <- test_boxplotInfo$rsID
  intronID <- test_boxplotInfo$phe_id
  qtl_pos <-  test_boxplotInfo$var_from
  
  # Double-check if coloc result exist
  coloc_results_subset <- coloc_results |> dplyr::filter(phe_id   %in% intronID)
  
  if(isEmpty(coloc_results_subset)){
    print("no-coloc-result exist")
  }
  
  
  subtype <- coloc_results_subset$subtype
  gwas_rsID <- coloc_results_subset$GWAS_rsID
  gwas_testPos <- str_split(coloc_results_subset$GWAS_pos,":",simplify = TRUE)[2] |> as.double()

  PP4_label_pbs <- as.numeric(coloc_results_subset$PBS_coloc_H4)
  PP4_label_fnf  <- as.numeric(coloc_results_subset$FNF_coloc_H4)
  
  # Load LD data
  ld_file <- paste0("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld/", chrom, "/", variantID, ".ld")
  ld_calc <- fread(ld_file) |> dplyr::select(c("SNP_B","R2")) |>
    dplyr::rename("var_id" = "SNP_B")
  
  # Load and process QTL data for both PBS and FNF
  pbs_norm_qtl <- pbs_norm_qtl_pc5_list[[chrom]]
  fnf_norm_qtl <- fnf_norm_qtl_pc4_list[[chrom]]
  
  pbs_qtl_region <- pbs_norm_qtl |> dplyr::filter(phe_id %in% intronID)
  fnf_qtl_region <- fnf_norm_qtl |> dplyr::filter(phe_id %in% intronID)
  
  pbs_qtl_pval_rsid <- left_join(pbs_qtl_region, maf_subset_rsID, by="var_id")
  fnf_qtl_pval_rsid <- left_join(fnf_qtl_region, maf_subset_rsID, by="var_id")
  
  leftjoin_pbs <- inner_join(pbs_qtl_pval_rsid, ld_calc, by = c("var_id" = "var_id"))
  leftjoin_fnf <- inner_join(fnf_qtl_pval_rsid, ld_calc, by = c("var_id" = "var_id"))
  
  #GWAS data subsetting --------------------------------------------------------
  
  gwas_data <- fread(paste0("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg38_LD/", subtype, "/leads/EUR_", 
                            subtype, "_leads_ld_final.csv"))
  
  gwas_data_subset <- gwas_data |> dplyr::filter(rsID  == gwas_rsID) |>
    dplyr::filter(!grepl("NA", `ldbuddy_CHR:hg38POS`)) |>
    dplyr::select("ldbuddy_CHR:hg38POS", "ldbuddy_R2", "CHR:hg38POS", "ldbuddy_rsID","rsID","CHR:hg38POS", "p")|>
    dplyr::mutate(
      chr_ldbuddy = paste0("chr", str_split(`ldbuddy_CHR:hg38POS`, ":", simplify = TRUE)[, 1]),
      pos_ldbuddy = as.numeric(str_split(`ldbuddy_CHR:hg38POS`, ":", simplify = TRUE)[, 2])
    ) |> dplyr::rename(gwas_rsID = "rsID") |>
    dplyr::rename(R2 = "ldbuddy_R2", var_chr = "chr_ldbuddy", var_from ="pos_ldbuddy",rsID="ldbuddy_rsID",nom_pval ="p")
  
  
  prepare_locus_plot <- function(data) {
    data |> 
      mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1))) |>
      mutate(LDgrp = addNA(LDgrp)) |>
      dplyr::rename(chrom = "var_chr",
                    pos = "var_from",
                    p = "nom_pval",
                    LD = "R2",
                    snp = "rsID") |>
      dplyr::select("chrom", "pos", "p", "snp", "LD", "LDgrp") |> 
      mutate(LDgrp = factor(LDgrp, levels = c(NA,"(0.8,1]","(0.6,0.8]", "(0.4,0.6]","(0.2,0.4]","(0,0.2]"), ordered = TRUE)) |>
      arrange(desc(LDgrp)) |>  # Sort by LDgrp so that higher LD values are at the end
      data.frame()
  }
  
  # Prepare plotting data for both PBS and FNF
  pbs_locus_plot <- prepare_locus_plot(leftjoin_pbs)
  fnf_locus_plot <- prepare_locus_plot(leftjoin_fnf)
  gwas_locus_plot <- prepare_locus_plot(gwas_data_subset)
  
  # Set plot region with adjustable zoom
  min_pos_gwas <- gwas_testPos - zoom_range
  max_pos_gwas <- gwas_testPos + zoom_range
  
  min_pos_qtl <- qtl_pos - zoom_range
  max_pos_qtl <- qtl_pos + zoom_range
  
  minregion <- min(min_pos_qtl,min_pos_gwas)
  maxregion <- max(max_pos_qtl,max_pos_gwas)
  
  # Set up plot parameters
  region_pg <- pgParams(assembly = "hg38", chrom = chrom,
                        chromstart = minregion,
                        chromend = maxregion,
                        x = x_start, y = y_start, width = width, height = height)
  
  # Calculate y-axis limits
  pbs_ylim <- ceiling(max(log10(pbs_locus_plot$p)*-1)) + 2
  fnf_ylim <- ceiling(max(log10(fnf_locus_plot$p)*-1)) + 2
  gwas_ylim <- ceiling(max(log10(gwas_locus_plot$p)*-1)) + 2
  ylim_pg <- max(pbs_ylim, fnf_ylim)
  
  # Function to create Manhattan plot
  create_manhattan_plot <- function(data, y_offset, label, PP4_label = NULL, label_color) {
    plot_height <- (height - 0.3) / 3  # Reduce height to make room for space between plots
    locus_plot <- plotManhattan(data = data,
                                params = region_pg,
                                range = c(0, ylim_pg),
                                fill = colorby("LDgrp",
                                               palette = colorRampPalette(c("#DD3931",  
                                                                            "#EEA741", 
                                                                            "#499A53",  
                                                                            "#98CDED",  
                                                                            "#262C74" 
                                               ))),
                                y = y_start + y_offset, height = plot_height,
                                snpHighlights = data.frame(snp = c(gwas_rsID, rsID),
                                                           pch = c(24, 23),
                                                           cex = c(0.75, 0.75),
                                                           col = c("black", "black")))
    
    # Add y-axis
    annoYaxis(plot = locus_plot, at = seq(0, ylim_pg, 2),
              axisLine = TRUE, fontsize = 6)
    
    # Add label
    plotText(label = label, 
             x = x_start + 0.05, y = y_start + y_offset, 
             just = c("left", "top"),
             fontfamily = "Helvetica", fontsize = 7, fontcolor = label_color )
    
    # Add PP4 label with error handling
    plotText(label = paste0("PP4=", round(PP4_label, digits = 3)),
               x = x_start + 0.05, y = y_start + y_offset + 0.15, 
               just = c("left", "top"),
               fontfamily = "Helvetica", fontsize = 6)
    
    return(locus_plot)
  }
  
  
GWAS_manhattan_plot <- function(data, y_offset, label) {
    plot_height <- (height - 0.3) / 3  # Reduce height to make room for space between plots
    locus_plot <- plotManhattan(data = data,
                                params = region_pg,
                                range = c(0, gwas_ylim),
                                fill = colorby("LDgrp",
                                               palette = colorRampPalette(c("#DD3931",  
                                                                            "#EEA741", 
                                                                            "#499A53",  
                                                                            "#98CDED",  
                                                                            "#262C74" 
                                               ))),
                                y = y_start + y_offset, height = plot_height,
                                snpHighlights = data.frame(snp = c(gwas_rsID, rsID),
                                                           pch = c(24, 23),
                                                           cex = c(0.75, 0.75),
                                                           col = c("black", "black")))
    
    
    # Add y-axis
    annoYaxis(plot = locus_plot, at = seq(0, gwas_ylim, 2),
              axisLine = TRUE, fontsize = 6)
    
    
    # Add label
    plotText(label =  label, 
             x = x_start + 0.05, y = y_start + y_offset, 
             just = c("left", "top"),
             fontfamily = "Helvetica", fontsize = 7)
    
    return(locus_plot)
  }
  
  # Create both PBS and FNF plots with space between them
  gwas_plot <- GWAS_manhattan_plot(gwas_locus_plot, 0, paste0("GWAS-", subtype))
  pbs_plot <- create_manhattan_plot(pbs_locus_plot, (height - 0.3) / 3+0.15 , "PBS sQTL",PP4_label_pbs ,"#0067B9")
  fnf_plot <- create_manhattan_plot(fnf_locus_plot, (height - 0.3) / 3 *2 +0.3, "FN-f sQTL" ,PP4_label_fnf ,  '#FCCE52')
  

  
  
  # Add y-axis label
  plotText(label = "-log10(p-value)", 
           x = x_start - 0.3, y = y_start + height/2, 
           rot = 90, fontsize = 6, just = "center",
           default.units = "inches", fontfamily = "Helvetica")
  
  # Add legend (adjusted position)
  legend_x <- x_start + width-0.8
  legend_y <- y_start # Moved up to accommodate space between plots
  plotLegend(legend = c("0.8 - 1.0", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "0.0 - 0.2"),
             #title="R²",
             fill = c("#DD3931", "#EEA741", "#499A53", "#98CDED", "#262C74"),
             x = legend_x, y = legend_y, width = 0.05, height = 0.3, border = FALSE,
             fontsize = 6
             )
  
  
  # Add rsID label
  
  grid.points(x = x_start + width -0.6 , y = y_start+(height - 0.3) / 3+0.15 , default.units = "native", pch = 24,
              size = unit(0.5, "char"))
  grid.points(x = x_start + width -0.6 , y = y_start+(height - 0.3) / 3+0.4 , default.units = "native", pch = 23,
              size = unit(0.5, "char"))
 
   plotText(label =  paste(gwas_rsID, "(GWAS)",sep="\n"), 
           x =x_start + width-0.45 , y =y_start+(height - 0.3) / 3+0.1, 
           just = c("left", "top"),
           fontfamily = "Helvetica", fontsize = 7,
           lineheight = 0.8)
   
   plotText(label = paste(rsID, "(sQTL)",sep="\n"), 
            x =x_start + width-0.45 , y =y_start+(height - 0.3) / 3+0.3, 
            just = c("left", "top"),
            fontfamily = "Helvetica", fontsize = 7,
            lineheight = 0.8)
   
   plotText(label = intronID, 
            x =x_start + (width/2), y =y_start-0.25, 
            just = c("center", "top"),
            fontfamily = "Helvetica", fontsize = 8)

  
  # Add gene track (adjusted position)
  gene_hl <- test_boxplotInfo$SYMBOL
  plotgenes <- plotGenes(params = region_pg, 
                         y = y_start + height + 0.1,
                         height = 0.4,
                         geneHighlights = data.frame("gene" = gene_hl,
                                                     "color" = "#37a7db"), fontsize = 6,geneOrder=gene_hl)
  annoGenomeLabel(plot = plotgenes, params = region_pg, fontsize = 6, 
                  y = y_start + height + 0.55)
}

# plot gardener to keep both plot at the same time------------------------------
pdf(file = "output/results_plots/coloc/main_figure5.pdf",   # The directory you want to save the file in
    width = 10.5, # The width of the plot in inches
    height = 9)

pageCreate(width = 10.5, height =9 , default.units = "inches", showGuides = FALSE)


# scatter plots of PBS and FNF--------------------------------------------------

plotText("a", x = 0.1, y = 0.1, just = c("center", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plot_gwas_qtl_locusZoom("RNF144B", "PBS", x_start =0.7, y_start = 0.6, width = 2.5, height = 3,zoom_range = 100000)


plotText("b", x = 3.5, y = 0.1, just = c("center", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

plot_gwas_qtl_locusZoom("GCAT", "PBS",x_start =4.1, y_start = 0.6, width = 2.5, height = 3,zoom_range = 250000)

plotText("c", x = 6.9, y = 0.1, just = c("center", "top"), fontfamily = "Helvetica",
       fontsize = 14, fontface = "bold")

plot_gwas_qtl_locusZoom("PBRM1","FNF",x_start =7.5, y_start = 0.6, width = 2.5, height = 3,zoom_range = 200000)


plotText("d", x = 0.1, y = 4.5, just =c("center", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plot_gwas_qtl_locusZoom(gene_name="WWP2", primary_dataset = "PBS", x_start =0.7, y_start = 5.1, width = 2.5, 
                        height = 3,zoom_range = 200000)


plotText("e", x = 3.5, y = 4.5, just = c("center", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

plot_gwas_qtl_locusZoom("COLGALT2", "PBS",x_start =4.1, y_start = 5.1, width = 2.5, height = 3,zoom_range = 180000)

plotText("f", x = 6.9, y = 4.5, just =c("center", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

plot_gwas_qtl_locusZoom(gene_name="HMGN1", primary_dataset ="PBS",
                        x_start =7.5, y_start = 5.1, width = 2.5, height = 3,zoom_range = 120000)

dev.off()
