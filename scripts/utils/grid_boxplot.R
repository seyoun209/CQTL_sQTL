log2fcColors <- c("+" = "#F2BC40", "-" = "#2057A7")
# Function to make grid-based boxplots of log2FoldChanges with highlighted genes
edited_boxplot_plot <- function(data, up_wilcox_test, down_wilcox_test, x, y, width, height,
                            default.units = "inches", just = c("left", "top")) {
  
  ## Create a function to draw boxplot with distribution
  makeBoxplot <- function(vec,
                          fill,
                          col,
                          yscale,
                          xpos = 0.2,
                          width = 0.25,
                          alpha = 0.6){
    
    ## Calculate median, q1, q3, 1.5*IQR
    boxplot_median <- median(vec)
    boxplot_q1 <- quantile(vec, 0.25)
    boxplot_q3 <- quantile(vec, 0.75)
    boxplot_IQR <- boxplot_q3 - boxplot_q1
    
    ## Define viewport
    
    vp <- viewport(x = unit(xpos, "npc"),
                   y = unit(0.5, "npc"),
                   width = width,
                   height = 1, just = c("left", "center"),
                   yscale = yscale,
                   clip = "off",
                   name = "boxplot")
    
    # Create boxplot background box
    
    boxplot_box <- polygonGrob(
      # bottom left, top left, top right, bottom right
      x = c(unit(0, "npc"), unit(0, "npc"), unit(1, "npc"), unit(1, "npc")),
      y = c(unit(boxplot_q1, "native"), unit(boxplot_q3, "native"),
            unit(boxplot_q3, "native"), unit(boxplot_q1, "native")),
      default.units = "native",
      gp = gpar(fill = alpha(colour = fill, alpha = alpha),
                col = col, lwd = 1)
    )
    
    ## Boxplot distribution lines
    # Median
    median_line <- segmentsGrob(x0 = 0,
                                y0 = unit(boxplot_median, "native"),
                                x1 = 1, y1 = unit(boxplot_median, "native"),
                                gp = gpar(lwd = 3, lineend = "butt", col = col))
    
    # Q1
    q1_line <- segmentsGrob(x0 = 0,
                            y0 = unit(boxplot_q1, "native"),
                            x1 = 1, y1 = unit(boxplot_q1, "native"),
                            gp = gpar(lwd = 1, lineend = "square", col = col))
    
    
    # Q3
    q3_line <- segmentsGrob(x0 = 0,
                            y0 = unit(boxplot_q3, "native"),
                            x1 = 1, y1 = unit(boxplot_q3, "native"),
                            gp = gpar(lwd = 1, lineend = "square", col = col))
    
    
    # Q3 + 1.5*IQR
    q3_iqr_line <- segmentsGrob(x0 = 0.25,
                                y0 = unit(boxplot_q3 + 1.5*boxplot_IQR, "native"),
                                x1 = 0.75,
                                y1 = unit(boxplot_q3 + 1.5*boxplot_IQR, "native"),
                                gp = gpar(lwd = 2, lineend = "square", col = col))
    
    # Q1 - 1.5*IQR
    q1_iqr_line <- segmentsGrob(x0 = 0.25,
                                y0 = unit(boxplot_q1 - 1.5*boxplot_IQR, "native"),
                                x1 = 0.75,
                                y1 = unit(boxplot_q1 - 1.5*boxplot_IQR, "native"),
                                gp = gpar(lwd = 2, lineend = "square", col = col))
    
    # Vertical line
    v_line <- segmentsGrob(x0 = 0.5, y0 = unit(boxplot_q1 - 1.5*boxplot_IQR, "native"),
                           x1 = 0.5, y1 = unit(boxplot_q3 + 1.5*boxplot_IQR, "native"),
                           gp = gpar(lwd = 1, lineend = "square", col = col))
    
    boxplot_gtree <- gTree(vp = vp, children = gList(boxplot_box,
                                                     median_line,
                                                     q1_line,
                                                     q3_line,
                                                     q3_iqr_line,
                                                     q1_iqr_line,
                                                     v_line))
    return(boxplot_gtree)
  }
  
  ## Create a function to add jittered points
  makePoints <- function(vec, yscale,
                         xpos = 0.2,
                         xwidth = 0.05,
                         pch = 21, size = 0.5,
                         fill = "black", alpha = 1, col = NA,
                         jitter = TRUE, jwidth = 0.5){
    get_xval_diff <- function(val){
      if (as.numeric(val) < 0) {
        xval_diff <- 1 - abs(as.numeric(val))
      } else {
        xval_diff <- 1 + as.numeric(val)
      }
    }
    
    ## Toggle jittered xposition
    if (jitter){
      set.seed(123)
      xvals <- unit(runif(length(vec), -jwidth, jwidth), "npc")
    } else {
      xvals <- unit(rep(0, length(vec)), "npc")
    }
    
    ## Define viewport
    vp <- viewport(x = unit(xpos, "npc"), y = unit(0.5, "npc"), width = xwidth, height = 1,
                   just = c("left", "center"), yscale = yscale,
                   name = "points")
    ## Points
    points <- pointsGrob(x = xvals, y = unit(vec, "native"),
                         pch = pch, size = unit(size, "npc"),
                         gp = gpar(fill = fill, col = col, alpha = alpha))
    
    # Create gTree
    points_gtree <- gTree(vp = vp, children = gList(points))
    
    ## Return xvalues for labels
    xvals_diff <- unlist(lapply(xvals, get_xval_diff))
    names(xvals_diff) <- names(vec)
    
    return(list(points_gtree, xvals_diff))
  }
  
  yscale <- c(-8,16)
  
  plotVP <- viewport(x = x, y = y, width = width, height = height,
                     default.units = default.units, just = just,
                     yscale = yscale, name = "plot_area")
  
  oa_boxplot_gtree <- gTree(name = "OA_boxplots",
                            vp = plotVP)
  
  ## Y- axis
  ylabs <- seq(-8,16,2 )
  
  # Axis line
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree,
                              child = segmentsGrob(x0 = 0.01,
                                                   y0 = 0, x1 = 0.01, y1 = 1,
                                                   gp = gpar(lwd = 0.75)))
  
  # axis text y
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree,
                              child = textGrob(label = ylabs,
                                               x = -0.01,
                                               y = unit(ylabs, "native"),
                                               gp = gpar(fontfamily = "Helvetica",
                                                         fontsize = 6)))
  # axis title y 
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree,
                              child = richtext_grob(text =
                                                      "log~2~(fold change)<br> in response to FN-f",
                                                    x = -0.05, y = 0.5, rot = 90,
                                                    gp = gpar(fontfamily = "Helvetica",
                                                              fontsize = 6)))
  # 0 line
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree,
                              child = segmentsGrob(x0 = 0.01,
                                                   y0 = unit(0, "native"),
                                                   x1 = 0.575,
                                                   y1 = unit(0, "native"),
                                                   default.units = "npc",
                                                   gp = gpar(col = "grey25", lty = 2)))
  
  ## Calculate x-positioning for groups
  n <- 10
  xpos <- seq(1,2*n, 2)/(2*n)
  xpos <- xpos + 0.05 # shift xpos
  
  ## Define barwidth
  barwidth <- 0.08
  
  
  ## Jittered points
  # Up in OA, not genes of interest
  oa_up <- makePoints(vec = data |>
                        filter(group == "Up in Het-KD" & is.na(highlight)) |>
                        pull(log2FoldChange),
                      yscale = yscale,
                      xpos = xpos[2],
                      size = 0.25, fill = "grey80", jwidth = 0.8)
  
  ## Add gTree to plot gTree
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = oa_up[[1]])
  
  # Up in OA, genes of interest
  up_vec <- data |> filter(group == "Up in Het-KD" & highlight == "up")
  up_vec_data <- up_vec[["log2FoldChange"]]
  names(up_vec_data) <- up_vec[["symbol"]]
  
  oa_up_highlight <- makePoints(vec = up_vec_data,
                                yscale = yscale,
                                xpos = xpos[2],
                                size = 0.25, fill = darken(log2fcColors[["+"]], 0.3), jwidth = 0.8)
  
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = oa_up_highlight[[1]])
  

  # Down in OA, not genes of interest
  oa_down <- makePoints(vec = data |>
                          filter(group == "Down in Het-KD" & is.na(highlight)) |>
                          pull(log2FoldChange),
                        yscale = yscale,
                        xpos = xpos[4],
                        size = 0.25, fill = "grey80", jwidth = 0.8)
  
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = oa_down[[1]])
  
  # Down in OA, genes of interest
  down_vec <- data |> filter(group == "Down in Het-KD" & highlight == "down")
  down_vec_data <- down_vec[["log2FoldChange"]]
  names(down_vec_data) <- down_vec[["symbol"]]
  
  oa_down_highlight <- makePoints(vec = down_vec_data,
                                  yscale = yscale,
                                  xpos = xpos[4],
                                  size = 0.25,
                                  fill = darken(log2fcColors[["-"]], 0.3), jwidth = 0.8)
  
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = oa_down_highlight[[1]])
  
  ## Boxplots
  up_boxplot <- makeBoxplot(data |>
                              filter(group == "Up in Het-KD") |>
                              pull(log2FoldChange),
                            yscale = yscale,
                            fill = log2fcColors[["+"]],
                            col = darken(log2fcColors[["+"]], 0.3),
                            xpos = xpos[2]-0.0375, width = barwidth)
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = up_boxplot)
  down_boxplot <- makeBoxplot(data |>
                                filter(group == "Down in Het-KD") |>
                                pull(log2FoldChange),
                              yscale = yscale,
                              fill = log2fcColors[["-"]],
                              col = darken(log2fcColors[["-"]], 0.3),
                              xpos = xpos[4]-0.0375, width = barwidth)
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = down_boxplot)
  
  # Get genes of interest
  goi <- data |>
    filter(!is.na(highlight)) |>
    pull(symbol)
  
  ## Upregulated
  goi_LFC_up <- data |>
    filter(group == "Up in Het-KD" & symbol %in% goi) |>
    arrange(log2FoldChange)
  
  ypos2 <- seq(0.5, 0.8, length.out = 5)
  ## Draw segments connecting y-positions (left)
  left_seg1 <- segmentsGrob(x0 =  (xpos[2] - (0.5*barwidth)) - 0.055, y0 = ypos2,
                            x1 =  (xpos[2] - (0.5*barwidth)) - 0.04, y1 = ypos2,
                            gp = gpar(lty=3, col = darken(log2fcColors[["+"]], 0.3)))
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = left_seg1)
  
  left_seg2 <- segmentsGrob(x0 = (xpos[2] - (0.5*barwidth)) - 0.04, y0 = ypos2,
                            x1 = (xpos[2] - (0.5*barwidth)) - 0.02,
                            y1 = unit(goi_LFC_up |> pull(log2FoldChange), "native"),
                            gp = gpar(lty=3, col = darken(log2fcColors[["+"]], 0.3)))
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = left_seg2)
  
  left_seg3 <- segmentsGrob(x0 = (xpos[2] - (0.5*barwidth)) - 0.02,
                            y0 = unit(goi_LFC_up  |> arrange(match(symbol, names(oa_up_highlight[[2]]))) |> pull(log2FoldChange), "native"),
                            x1 = (xpos[2] - 0.5*barwidth) + oa_up_highlight[[2]]*barwidth*0.5,
                            y1  = unit(goi_LFC_up  |> arrange(match(symbol, names(oa_up_highlight[[2]]))) |> pull(log2FoldChange), "native"),
                            gp = gpar(lty = 3, col = darken(log2fcColors[["+"]], 0.3), lineend = "butt"))
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = left_seg3)
  
  ## Up gene labels
  upGenes <- textGrob(label = goi_LFC_up |> pull(symbol),
                      x = (xpos[2] - (0.5*barwidth)) - 0.06, y = ypos2, just = c("right", "center"),
                      gp = gpar(col = darken(log2fcColors[["+"]], 0.3), fontsize = 6,
                                fontfamily = "Helvetica"))
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = upGenes)
  
  ## Downregulated
  goi_LFC_down <- data |>
    filter(group == "Down in Het-KD" & symbol %in% goi) |>
    arrange(log2FoldChange)
  
  
  ypos2 <- seq(0.05, 0.35, length.out = 5)
  
  ## Draw segments connecting y-positions (right)
  right_seg1 <- segmentsGrob(x0 = (xpos[4] + (0.5*barwidth)) + 0.02,
                             y0 = unit(goi_LFC_down |>
                                         arrange(match(symbol, names(oa_down_highlight[[2]]))) |>
                                         pull(log2FoldChange), "native"),
                             x1 = (xpos[4] - 0.5*barwidth) + oa_down_highlight[[2]]*barwidth*0.55,
                             y1  = unit(goi_LFC_down |> arrange(match(symbol, names(oa_down_highlight[[2]])))
                                        |> pull(log2FoldChange), "native"),
                             gp = gpar(lty = 3, col = darken(log2fcColors[["-"]], 0.3),
                                       lineend = "butt"))
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = right_seg1)
  
  right_seg2 <- segmentsGrob(x0 = (xpos[4] + (0.5*barwidth)) + 0.02,
                             y0 = unit(goi_LFC_down |> pull(log2FoldChange), "native"),
                             x1 = (xpos[4] + (0.5*barwidth)) + 0.04,
                             y1 = ypos2,
                             gp = gpar(lty=3, col = darken(log2fcColors[["-"]], 0.3)))
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = right_seg2)
  
  right_seg3 <- segmentsGrob(x0 =  (xpos[4] + (0.5*barwidth)) + 0.04, y0 = ypos2,
                             x1 =  (xpos[4] + (0.5*barwidth)) + 0.055, y1 = ypos2,
                             gp = gpar(lty=3, col = darken(log2fcColors[["-"]], 0.3)))
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = right_seg3)
  
  
  ## Down gene labels
  downGenes <- textGrob(label = goi_LFC_down |> pull(symbol),
                        x = (xpos[4] + (0.5*barwidth)) + 0.06, y = ypos2, just = c("left", "center"),
                        gp = gpar(col = darken(log2fcColors[["-"]], 0.3), fontsize = 6,
                                  fontfamily = "Helvetica"))
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = downGenes)
  
  ## X-axis labels
  xaxis_labels <- textGrob(label = c("Up in Het-KD", "Down in Het-KD"),
                           x = xpos[c(2, 4)], y =-0.025, just = c("center", "top"),
                           gp = gpar(lineheight = 0.9, fontfamily = "Helvetica",
                                     fontsize = 8,
                                     col = c(darken(log2fcColors[["+"]], 0.3), darken(log2fcColors[["-"]], 0.3))))
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = xaxis_labels)
  
  ## Significance stars from Wilcox testing enrichment
  # Upregulated
  pval_up <- sprintf("p-value: %s", format.pval(up_wilcox_test$p.value, digits = 2))
  pval_down <- sprintf("p-value: %s", format.pval(down_wilcox_test$p.value, digits = 2))
  
  if (up_wilcox_test$p.value < 0.05){
    up_test_star <- textGrob(label = pval_up, x = xpos[2], y = 0.95,
                             gp = gpar(fontface = "bold", fontfamily = "Helvetica",
                                       fontsize = 6))
  } else {
    up_test_star <- textGrob(label = "ns", x = xpos[2], y = 0.95,
                             gp = gpar(fontface = "bold", fontfamily = "Helvetica"))
  }
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = up_test_star)
  # Downregulated
  if (down_wilcox_test$p.value < 0.05){
    down_test_star <- textGrob(label = pval_down, x = xpos[4], y = 0.95,
                               gp = gpar(fontface = "bold", fontfamily = "Helvetica",
                                         fontsize = 6))
  } else {
    down_test_star <- textGrob(label = "ns", x = xpos[4], y = 0.95,
                               gp = gpar(fontface = "bold", fontfamily = "Helvetica"))
  }
  oa_boxplot_gtree <- addGrob(oa_boxplot_gtree, child = down_test_star)
  
  
  grid.draw(oa_boxplot_gtree)
  
}

# Function to calculate sample-level l2fc for a gene in a count matrix
get_sample_l2fc <- function(gene, countMatrix){
  
  # Extract row of gene from countMatrix
  gene_counts <- countMatrix[gene,]
  
  # Convert to dataframe and extract donors/conditions into separate columns
  donor_gene_counts <- data.frame(gene_counts) |>
    rownames_to_column(var = "Sample") |>
    separate_wider_delim("Sample",
                         delim = "_",
                         names = c(NA, "Donor", NA, "Condition", NA, NA)) |>
    # Group by each donor and calculate l2FC
    group_by(Donor) |>
    summarize(log2FC =
                log2(gene_counts[Condition == "FNF"]/gene_counts[Condition == "CTL"])) |>
    ungroup() |>
    mutate(ENSEMBL = gene)
  
  return(donor_gene_counts)
}
