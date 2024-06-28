# Utils for the sqtl 
## Author: Seyoun Byun
## Date: 04.22.2024
## Edited:
#-------------------------------------------------------------------------------
library(rtracklayer)
library(AnnotationDbi)
setwd("/work/users/s/e/seyoun/CQTL_sQTL")

venn_font <- function(p, font){
  
  grep_grob <- function(gt, lab){
    which(sapply(gt, function(x) grepl(lab, x$name)))
  }
  
  p2 <- ggplot_gtable(ggplot_build(p))
  mygrobs <- p2$grobs
  # Break down grobs into panel, venn object, and text pieces
  panel_grob <- mygrobs[[grep_grob(mygrobs, "panel")]]
  venn_grob <- panel_grob$children[[grep_grob(panel_grob$children, "venn")]]
  text_grobs <- venn_grob$children[grep_grob(venn_grob$children, "text")]
  # Make both new font family
  text_grobs <- do.call(grid::gList, 
                        lapply(text_grobs, 
                               function(x) {x$gp$fontfamily <- font; 
                               x}))
  # Make titles bold
  text_grobs[[1]]$gp$fontface <- "bold"
  
  # Add grobs back
  venn_grob$children[grep_grob(venn_grob$children, "text")] <- text_grobs
  panel_grob$children[[grep_grob(panel_grob$children, "venn")]] <- venn_grob
  mygrobs[[grep_grob(mygrobs, "panel")]] <- panel_grob
  p2$grobs <- mygrobs
  grid::grid.newpage()
  grid::grid.draw(p2)
  
}


hg38_intron <- fread("/work/users/s/e/seyoun/CQTL_sQTL/tools/leafcutter/gencode_hg38_all_introns.bed.gz")
colnames(hg38_intron) <- c("chr","start","end","symbol","ensg","strand","enst","id","region","info")
#hg38_intron_filtered <- hg38_intron %>%
#  filter(!startsWith(symbol, "ENSG"))

hg38_intron_sub <- hg38_intron %>%
  mutate(genomicLoc = paste(chr, start, end, sep = ":")) %>%
  dplyr::select(genomicLoc, ensg)

hg38_intron_sub_select_first <- hg38_intron_sub %>%
  group_by(genomicLoc) %>%
  arrange(ensg) %>%  # Ensure that the data is sorted by ensg within each genomicLoc
  dplyr::slice(1) %>%       # Select the first ensg for each genomicLoc
  ungroup()


leafcutter_pheno <- fread("output/gtex_cluster/ctl_fnf.leafcutter.phenotype_groups.txt",header=F)
colnames(leafcutter_pheno) <- c("id_ensg","ensg")
leafcutter_pheno_mutated <- leafcutter_pheno %>%
  mutate(
    chr=gsub("\\..*$", "", sapply(strsplit(as.character(id_ensg), ":"), `[`, 1)),
    start=gsub("\\..*$", "", sapply(strsplit(as.character(id_ensg), ":"), `[`, 2)),
    end=gsub("\\..*$", "", sapply(strsplit(as.character(id_ensg), ":"), `[`, 3)),
    clusterID=gsub("\\..*$", "", sapply(strsplit(as.character(id_ensg), ":"), `[`, 4)),
    phe_id=sapply(strsplit(as.character(id_ensg), ":"), function(x) paste(x[1:4], collapse = ":")),
    genomicLoc=sapply(strsplit(as.character(id_ensg), ":"), function(x) paste(x[1:3], collapse = ":"))
  )
leafcutter_pheno_subset <- leafcutter_pheno_mutated %>%
  dplyr::select(ensg, clusterID) %>%
  group_by(clusterID) %>%
  arrange(ensg) %>%  # Ensure that the data is sorted by ensg within each genomicLoc
  dplyr::slice(1) %>%       # Select the first ensg for each genomicLoc
  ungroup()




gtf_path <- "/work/users/s/e/seyoun/Ref/genome/gencode.v45.annotation.gtf"
gtf_data <- import(gtf_path)
gene_info <- gtf_data[gtf_data$type == "gene"]

genes_df <- data.frame(
  gene_id = mcols(gtf_data)$gene_id,
  gene_name = mcols(gtf_data)$gene_name,
  stringsAsFactors = FALSE
)

genes_df_unique <- genes_df %>%
  dplyr::rename(ensg = gene_id) %>%
  dplyr::distinct(ensg, gene_name, .keep_all = TRUE)  


txdb <- loadDb("/work/users/s/e/seyoun/crispr/02.test_seq/gencode.v45.annotation.TxDb")
txdb_genes <- genes(txdb)


annotate_genes <- function(cond_pbs, hg38_intron_sub_select_first, leafcutter_pheno_subset, txdb_genes) {
  # Processing initial data
  conditional_pbs <- cond_pbs %>%
    mutate(
      clusterID = gsub("\\..*$", "", sapply(strsplit(as.character(phe_id), ":"), `[`, 4)),
      genomicLoc = sapply(strsplit(as.character(phe_id), ":"), function(x) paste(x[1:3], collapse = ":"))
    )
  
  # First annotation trial
  conditional_sGene_pbs_anno_v1 <- conditional_pbs %>%
    left_join(hg38_intron_sub_select_first, by = "genomicLoc") 
  conditional_sGene_pbs_NO_anno_v1 <- conditional_sGene_pbs_anno_v1[is.na(conditional_sGene_pbs_anno_v1$ensg),]
  
  # Second annotation trial
  conditional_sGene_pbs_anno_v2 <- conditional_sGene_pbs_NO_anno_v1 %>%
    left_join(leafcutter_pheno_subset, by = "clusterID") %>%
    mutate(ensg = coalesce(ensg.x, ensg.y)) %>%
    dplyr::select(-ensg.x, -ensg.y)
  conditional_sGene_pbs_NO_anno_v2 <- conditional_sGene_pbs_anno_v2[is.na(conditional_sGene_pbs_anno_v2$ensg),]
  
  # Third annotation trial with genomic ranges
  conditional_sGene_pbs_gr <- GRanges(
    seqnames = conditional_sGene_pbs_NO_anno_v2$phe_chr,
    ranges = IRanges(start = conditional_sGene_pbs_NO_anno_v2$phe_from-1, end = conditional_sGene_pbs_NO_anno_v2$phe_to)
  )
  overlaps <- mergeByOverlaps(conditional_sGene_pbs_gr, txdb_genes, type="within") 
  overlaps_df <- data.frame(
    gene_id = overlaps$gene_id,
    genomicLoc = paste0(seqnames(overlaps$conditional_sGene_pbs_gr), ":", start(overlaps$conditional_sGene_pbs_gr), ":", end(overlaps$conditional_sGene_pbs_gr))
  )
  overlaps_df_selected_first <- overlaps_df %>% 
    group_by(genomicLoc) %>%
    arrange(gene_id) %>%
    dplyr::slice(1) %>%
    ungroup()
  
  conditional_sGene_pbs_anno_v3 <- conditional_sGene_pbs_NO_anno_v2 %>%
    left_join(overlaps_df_selected_first, by = "genomicLoc") %>%
    mutate(ensg = coalesce(gene_id, ensg)) %>%
    dplyr::select(-gene_id) 
  
  # Merging all
  conditional_pbs_final <- rbind(
    conditional_sGene_pbs_anno_v1[!is.na(conditional_sGene_pbs_anno_v1$ensg),],
    conditional_sGene_pbs_anno_v2[!is.na(conditional_sGene_pbs_anno_v2$ensg),], 
    conditional_sGene_pbs_anno_v3
  )
  
  return(conditional_pbs_final)
}


#making a pie chart for the response QTL
create_pie_chart <- function(data, colors, group) {
  # Compute percentages and cumulative positions
  data$Percentage <- round((data$Freq / sum(data$Freq)), 3) * 100
  data$ymax <- cumsum(data$Freq / sum(data$Freq))
  data$ymin <- c(0, head(data$ymax, n = -1))
  data$labelPosition <- (data$ymax + data$ymin) / 2
  
  # Create the pie chart
  p <- ggplot(data, aes(x = "", y = Percentage, fill = group)) +
    geom_bar(stat = "identity") +
    theme_void() +
    scale_fill_manual(values = colors) +
    theme(legend.position = "none", plot.background = element_rect(fill = 'transparent', color = NA))
  
  if (group == "PBS") {
    p <- p + coord_polar("y", direction = -1, start = 0)
  } else if (group == "FNF") {
    p <- p + coord_polar("y", direction = 1, start = 0)
  }
  
  return(p)
}



#Dot from for the line plot 
get_y_for_median <- function(data, median_x, binwidth) {
  bins <- hist(data$abs_distance, breaks = seq(0, max(data$abs_distance) + binwidth, by = binwidth), plot = FALSE)
  bin_index <- which(bins$mids >= median_x)[1]
  return(bins$counts[bin_index])
}

