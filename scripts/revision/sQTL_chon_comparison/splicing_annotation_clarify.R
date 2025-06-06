setwd("/work/users/s/e/seyoun/CQTL_sQTL/output")

library(stringr)
library(data.table)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(GenomicFeatures)
library(leafcutter)
library(BSgenome)
library(rtracklayer)
library(scales)
library(ggplot2)
library(fuzzyjoin)
library(ggtext) 
library(plotgardener)
library(BSgenome.Hsapiens.UCSC.hg38)

#-----------------------------------------------------------------------
# Load the all unique junctions from Leafcutter output
#-----------------------------------------------------------------------
# Load the reference genome
# 1.Question is regarding that if I included all unique junctions from Leafcutter and how exactly define previously annotated junctions
pbs_fnf <- load("clu_fnf/PBSvsFNF.Rdata")
load("clu_fnf/introns_fnf_joinAll") #introns_fnf_pval_include

#Extract all unique junctions from Leafcutter introns

#introns_with_strand <- introns %>%
#  mutate(strand = str_extract(clusterID, "[+-]$"))
#unique_junctions <- introns_with_strand %>%
#  distinct(chr, start, end, strand)

#total_unique_junctions <- nrow(unique_junctions)
# Create GRanges object for junctions
#junction_gr <- GRanges(
#  seqnames = unique_junctions$chr,
#  ranges = IRanges(start = unique_junctions$start, end = unique_junctions$end),
#  strand = unique_junctions$strand
#)

# Load the reference files you have available
# Import BED files using read.table with flexible column specifications
read_bed_to_granges <- function(file_path) {
  # Read the BED file with column types explicitly set to character
  df <- read.table(file_path, 
                  sep="\t", 
                  header=FALSE, 
                  stringsAsFactors=FALSE,
                  colClasses="character",
                  quote="")
  
  # Assign column names based on BED format
  colnames <- c("chromosome", "start", "end", "gene_name", "gene_id", "strand", 
                "transcript_id", "intron_number", "transcript_type", "annotations")
  
  # Ensure we only use available columns
  colnames <- colnames[1:min(length(colnames), ncol(df))]
  colnames(df)[1:length(colnames)] <- colnames
  
  # Convert to GRanges
  gr <- GRanges(
    seqnames = df$chromosome,
    ranges = IRanges(start = as.numeric(df$start), end = as.numeric(df$end)),
    strand = df$strand
  )
  
  # Add all available metadata columns
  metadata_cols <- setdiff(colnames(df), c("chromosome", "start", "end", "strand"))
  mcols(gr) <- df[, metadata_cols, drop=FALSE]
  
  return(gr)
}

# Use the custom function to load the files
gencode_all_introns <- read_bed_to_granges("/work/users/s/e/seyoun/CQTL_sQTL/tools/leafcutter/gencode_hg38_all_introns.bed.gz")
gencode_fiveprime <- read_bed_to_granges("/work/users/s/e/seyoun/CQTL_sQTL/tools/leafcutter/gencode_hg38_fiveprime.bed.gz")
gencode_threeprime <- read_bed_to_granges("/work/users/s/e/seyoun/CQTL_sQTL/tools/leafcutter/gencode_hg38_threeprime.bed.gz")
genome <- BSgenome.Hsapiens.UCSC.hg38

# Extract all unique junctions from your counts or ratios data
all_junctions <- rownames(counts)
total_junctions <- length(all_junctions)

# Parse junction coordinates
junction_parts <- strsplit(all_junctions, ":")
junction_coords <- data.frame(
  chr = sapply(junction_parts, function(x) x[1]),
  start = as.numeric(sapply(junction_parts, function(x) x[2])),
  end = as.numeric(sapply(junction_parts, function(x) x[3])),
  cluster = sapply(junction_parts, function(x) x[4])
)

# Extract strand information
junction_coords$strand <- substr(sapply(junction_coords$cluster, 
                                        function(x) sub(".*_", "", x)), 
                                nchar(sapply(junction_coords$cluster, 
                                             function(x) sub(".*_", "", x))), 
                                nchar(sapply(junction_coords$cluster, 
                                             function(x) sub(".*_", "", x))))

# Create GRanges object
junction_gr <- GRanges(
  seqnames = junction_coords$chr,
  ranges = IRanges(start = junction_coords$start, end = junction_coords$end),
  strand = junction_coords$strand
)

# Define donor and acceptor sites
donors <- GRanges(
  seqnames = seqnames(junction_gr),
  ranges = IRanges(start = start(junction_gr), width = 1),
  strand = strand(junction_gr)
)

acceptors <- GRanges(
  seqnames = seqnames(junction_gr),
  ranges = IRanges(start = end(junction_gr), width = 1),
  strand = strand(junction_gr)
)


# Calculate donor and acceptor overlaps
donor_annotated <- countOverlaps(donors, gencode_fiveprime, ignore.strand=FALSE) > 0
acceptor_annotated <- countOverlaps(acceptors, gencode_threeprime, ignore.strand=FALSE) > 0

# Calculate different categories
fully_annotated <- donor_annotated & acceptor_annotated
cryptic_fiveprime <- !donor_annotated & acceptor_annotated
cryptic_threeprime <- donor_annotated & !acceptor_annotated
cryptic_unanchored <- !donor_annotated & !acceptor_annotated

# Calculate percentages
percent_fully_annotated <- sum(fully_annotated) / total_junctions * 100
percent_cryptic_fiveprime <- sum(cryptic_fiveprime) / total_junctions * 100
percent_cryptic_threeprime <- sum(cryptic_threeprime) / total_junctions * 100
percent_cryptic_unanchored <- sum(cryptic_unanchored) / total_junctions * 100

cat(sprintf("Total unique junctions: %d\n", total_junctions))
cat(sprintf("Fully annotated junctions: %d (%.2f%%)\n", 
            sum(fully_annotated), percent_fully_annotated))

# 2. Additional validation of unannotated junctions (Canonical splice junction that is unannotated) and potential oligo(dT)
# Get unannotated junctions
unannotated_gr <- junction_gr[!(fully_annotated)]

# Check for canonical splice sites
check_canonical_sites <- function(gr) {
  canonical <- logical(length(gr))
  
  for (i in seq_along(gr)) {
    chr <- as.character(seqnames(gr[i]))
    if (!chr %in% seqnames(genome)) next
    
    strand_i <- as.character(strand(gr[i]))
    
    if (strand_i == "+") {
      # For + strand, check GT at donor (5') and AG at acceptor (3')
      donor_seq <- getSeq(genome, chr, start(gr[i]), start(gr[i])+1)
      acceptor_seq <- getSeq(genome, chr, end(gr[i])-1, end(gr[i]))
      
      canonical[i] <- as.character(donor_seq) == "GT" && as.character(acceptor_seq) == "AG"
    } else {
      # For - strand, check CT at donor (3') and AC at acceptor (5')
      donor_seq <-  BSgenome::getSeq(genome, chr, end(gr[i])-1, end(gr[i]))
      acceptor_seq <-  BSgenome::getSeq(genome, chr, start(gr[i]), start(gr[i])+1)
      
      canonical[i] <- as.character(donor_seq) == "CT" && as.character(acceptor_seq) == "AC"
    }
  }
  
  return(canonical)
}

canonical_sites <- check_canonical_sites(unannotated_gr)
percent_canonical <- sum(canonical_sites, na.rm=TRUE) / length(unannotated_gr) * 100



# Check for potential oligo(dT) artifacts
check_oligo_dt <- function(gr, threeprime_refs, distance = 250) {
  near_threeprime <- logical(length(gr))
  
  for (i in seq_along(gr)) {
    chr_i <- as.character(seqnames(gr[i]))
    strand_i <- as.character(strand(gr[i]))
    
    relevant_sites <- threeprime_refs[
      seqnames(threeprime_refs) == chr_i & 
      strand(threeprime_refs) == strand_i
    ]
    
    if (length(relevant_sites) == 0) next
    
    # For + strand, check 3' end; for - strand, check 5' end
    if (strand_i == "+") {
      min_dist <- min(abs(end(gr[i]) - start(relevant_sites)))
    } else {
      min_dist <- min(abs(start(gr[i]) - start(relevant_sites)))
    }
    
    near_threeprime[i] <- min_dist < distance
  }
  
  return(near_threeprime)
}

oligo_dt_results <- check_oligo_dt(unannotated_gr, gencode_threeprime)
percent_oligo_dt <- sum(oligo_dt_results, na.rm=TRUE) / length(unannotated_gr) * 100

cat(sprintf("Unannotated junctions with canonical splice sites: %d (%.2f%%)\n", 
            sum(canonical_sites, na.rm=TRUE), percent_canonical))
cat(sprintf("Potential oligo(dT) artifacts: %d (%.2f%%)\n", 
            sum(oligo_dt_results, na.rm=TRUE), percent_oligo_dt))


save(unannotated_gr,canonical_sites,oligo_dt_results,
 file="/work/users/s/e/seyoun/CQTL_sQTL/output/revision/data/canonical.Rdata")



#------------------------------------------------------------------------
# read from STAR SJ.out.tab
#------------------------------------------------------------------------

sj_files <- list.files(path="align",pattern = "*.SJ.out.tab", full.names = TRUE)

# Exclude files containing "OA", "AM7352", or "AM7244"
sj_files_filtered <- sj_files[!(grepl("OA", sj_files) | 
                               grepl("AM7352", sj_files) | 
                               grepl("AM7244", sj_files))]

cat(sprintf("Filtered files (CTL and FNF): %d\n", length(sj_files_filtered)))
valid_chroms <- c(paste0("chr", 1:22), "chrX", "chrY")
sj_list <- lapply(sj_files, function(file) {
  dt <- fread(file)
  dt <- dt[V1 %in% valid_chroms, ]
  return(dt)
})

# Combine and remove duplicates
sj_all <- rbindlist(sj_list)
sj_all_unique <- unique(sj_all, by = c("V1", "V2", "V3", "V4","V5", "V6"))
cat(sprintf("Total unique junctions (chr1-22, X, Y): %d\n", nrow(sj_all)))


introns_dt <- as.data.table(introns)
# Extract strand from clusterID (+ or -)
introns_dt[, strand := sub(".*_(.)$", "\\1", clusterID)]
introns_dt[, V4 := ifelse(strand == "+", 1, 2)]  # Map to STAR: 1 = +, 2 = -

# Ensure matching column names
setnames(introns_dt, c("chr", "start", "end"), c("V1", "V2", "V3"))

# Convert LeafCutter introns to match STAR SJ coordinates
introns_dt_star_format <- copy(introns_dt)
# For + strand junctions: adjust start (+1)
# For - strand junctions: adjust end (-1)
# Correct adjustment for both strands
introns_dt_star_format <- copy(introns_dt)
introns_dt_star_format_cryptic <- introns_dt_star_format |> filter(verdict %in% c("cryptic_fiveprime", "cryptic_threeprime", "cryptic_unanchored"))
introns_dt_start_format_annotated <- introns_dt_star_format |> filter(verdict %in% c("annotated","novel annotated pair"))
# Fuzzy matching approach to account for small coordinate differences

leafcutter_star_comparison <- introns_dt_star_format_cryptic[, {
  # Find potential matches in STAR data with some tolerance
  star_matches <- sj_all_unique[
    V1 == .BY$V1 &                        # Same chromosome
    abs(V2 - .BY$V2) <= 2 &               # Start position within 2bp
    abs(V3 - .BY$V3) <= 2 &               # End position within 2bp  
    V4 == .BY$V4                          # Same strand
  ]
  
  if(nrow(star_matches) > 0) {
    # If multiple matches, take the closest one
    if(nrow(star_matches) > 1) {
      star_matches[, dist := abs(V2 - .BY$V2) + abs(V3 - .BY$V3)]
      star_matches <- star_matches[order(dist)][1]
    }
    # Return the match with original LeafCutter data
    cbind(.SD, star_matches[, .(STAR_V2=V2, STAR_V3=V3, V5, V6, V7, V8, V9)])
  } else {
    # No match found
    cbind(.SD, data.table(STAR_V2=NA_integer_, STAR_V3=NA_integer_, 
                         V5=NA_integer_, V6=NA_integer_, 
                         V7=NA_integer_, V8=NA_integer_, V9=NA_integer_))
  }
}, by=.(V1, V2, V3, V4)]

leafcutter_star_comparison_annotated <- introns_dt_start_format_annotated[, {
  # Find potential matches in STAR data with some tolerance
  star_matches <- sj_all_unique[
    V1 == .BY$V1 &                        # Same chromosome
    abs(V2 - .BY$V2) <= 2 &               # Start position within 2bp
    abs(V3 - .BY$V3) <= 2 &               # End position within 2bp  
    V4 == .BY$V4                          # Same strand
  ]
  
  if(nrow(star_matches) > 0) {
    # If multiple matches, take the closest one
    if(nrow(star_matches) > 1) {
      star_matches[, dist := abs(V2 - .BY$V2) + abs(V3 - .BY$V3)]
      star_matches <- star_matches[order(dist)][1]
    }
    # Return the match with original LeafCutter data
    cbind(.SD, star_matches[, .(STAR_V2=V2, STAR_V3=V3, V5, V6, V7, V8, V9)])
  } else {
    # No match found
    cbind(.SD, data.table(STAR_V2=NA_integer_, STAR_V3=NA_integer_, 
                         V5=NA_integer_, V6=NA_integer_, 
                         V7=NA_integer_, V8=NA_integer_, V9=NA_integer_))
  }
}, by=.(V1, V2, V3, V4)]

save(introns_dt_star_format_cryptic
sj_all_unique, leafcutter_star_comparison,
     file="/work/users/s/e/seyoun/CQTL_sQTL/output/revision/data/leafcutter_star_comparison.Rdata")


save(introns_dt_start_format_annotated, 
sj_all_unique, leafcutter_star_comparison_annotated,leafcutter_star_comparison_annotated_fi,
     file="/work/users/s/e/seyoun/CQTL_sQTL/output/revision/data/leafcutter_star_comparison_annotated.Rdata")

load("/work/users/s/e/seyoun/CQTL_sQTL/output/revision/data/leafcutter_star_comparison.Rdata")
load("/work/users/s/e/seyoun/CQTL_sQTL/output/revision/data/leafcutter_star_comparison_annotated.Rdata")

leafcutter_star_comparison_annotated_na <- leafcutter_star_comparison_annotated[is.na(V5)] |> dplyr::select(-c(STAR_V2,STAR_V3,V5,V6,V7,V8,V9))

leafcutter_star_na_reannot <- leafcutter_star_comparison_annotated_na[, {
  # Find potential matches in STAR data with some tolerance
  star_matches <- sj_all_unique[
    V1 == .BY$V1 &                        # Same chromosome
    abs(V2 - .BY$V2) <= 2 &               # Start position within 2bp
    abs(V3 - .BY$V3) <= 2                   # End position within 2bp
  ]
  
  if(nrow(star_matches) > 0) {
    # If multiple matches, take the closest one
    if(nrow(star_matches) > 1) {
      star_matches[, dist := abs(V2 - .BY$V2) + abs(V3 - .BY$V3)]
      star_matches <- star_matches[order(dist)][1]
    }
    # Return the match with original LeafCutter data
    cbind(.SD, star_matches[, .(STAR_V2=V2, STAR_V3=V3, V5, V6, V7, V8, V9)])
  } else {
    # No match found
    cbind(.SD, data.table(STAR_V2=NA_integer_, STAR_V3=NA_integer_, 
                         V5=NA_integer_, V6=NA_integer_, 
                         V7=NA_integer_, V8=NA_integer_, V9=NA_integer_))
  }
}, by=.(V1, V2, V3, V4)]

leafcutter_star_comparison_annotated_fi <- rbind(leafcutter_star_comparison_annotated[!is.na(V5)], leafcutter_star_na_reannot)


#-----------------------------------------------------------------------
# 2. Make a bar plot for total and also differential only 

## Total comparison

# Define mapping for V5 values to motif labels
motif_labels <- c(
  "0" = "non-canonical",
  "1" = "GT/AG",
  "2" = "CT/AC",
  "3" = "GC/AG",
  "4" = "CT/GC",
  "5" = "AT/AC",
  "6" = "GT/AT"
)

# Summarize cryptic data and filter out any NA
cryptic_motif_summary <- leafcutter_star_comparison %>% 
  group_by(V5) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(motif = factor(as.character(V5),
                        levels = names(motif_labels),
                        labels = motif_labels[names(motif_labels)]),
         Group = "Cryptic")



# Summarize annotated data and filter out any NA
annotated_motif_summary <- leafcutter_star_comparison_annotated_fi %>% 
  group_by(V5) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(motif = factor(as.character(V5),
                        levels = names(motif_labels),
                        labels = motif_labels[names(motif_labels)]),
         Group = "Annotated")

# Combine the summaries
combined_summary <- bind_rows(cryptic_motif_summary, annotated_motif_summary) %>%
  mutate(Group = factor(Group, levels = c("Annotated", "Cryptic")))

# Define full motif levels to enforce consistent x-axis
full_motifs <- c( "GT/AG", "CT/AC", "GC/AG", "CT/GC", "AT/AC", "GT/AT","non-canonical")



# Complete the data so each motif appears in both groups (fill missing with count = 0)
combined_summary_complete <- combined_summary %>%
  complete(motif = full_motifs, Group, fill = list(count = 0)) %>%
  group_by(Group) %>%
  mutate(total = sum(count),
         perc = count / total * 100) %>%
  ungroup()

# Create the barplot with percentage on the y-axis
oligo_barplot_full <- ggplot(combined_summary_complete, aes(x = motif, y = perc, fill = Group)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = sprintf("%.2f%%", perc)),
            position = position_dodge(width = 0.9),
            vjust = -0.5,
            size = 2.5,
            family = "Helvetica") +
  labs(
    x = "Intron Motif",
    y = "Percentage of splice intron junctions"
  ) +
  scale_y_continuous(
    name = "Percentage of splice intron junctions",
    limits = c(0, 60),
    breaks = seq(0, 60, by = 10),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(values = c("Cryptic" = "grey65", "Annotated" = "grey35")) +
  scale_x_discrete(
    limits = full_motifs
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(linewidth = 0.25),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black", linewidth = 0.25),
    axis.ticks.length.y = unit(-0.1, "cm"),
    axis.title.x = element_markdown(size = 8, family = "Helvetica", margin = margin(t = 5)),
    axis.title.y = element_markdown(size = 8, family = "Helvetica", margin = margin(r = 5)),
    text = element_text(family = "Helvetica"),
    axis.text.y = element_text(color = "black", size = 6),
    axis.text.x = element_text(color = "black", size = 6, margin = margin(t = 5), angle = 45, hjust = 1),
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position.inside = c(0.9, 0.9),
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 6),
    panel.spacing.y = unit(0.5, "cm"),
    plot.title = element_blank()
  )
save(oligo_barplot_full,
     file = "revision/plots/splicing_oligo_barplot_full.rda")

load("revision/plots/splicing_oligo_barplot_full.rda")


## Differentail version
introns_fnf_sig <- introns_fnf_pval_include |>  dplyr::filter(p.adjust < 0.05) |> 
  dplyr::filter(abs(deltapsi_batch) > 0.15) |> mutate(junction_id = paste(chr, start, end, sep=":"))

# Summarize cryptic data and filter out any NA
cryptic_motif_summary <- leafcutter_star_comparison %>% 
  mutate(junction_id = paste(V1, V2, V3, sep=":")) %>%
  filter(junction_id %in% introns_fnf_sig$junction_id) %>%
  group_by(V5) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(motif = factor(as.character(V5),
                        levels = names(motif_labels),
                        labels = motif_labels[names(motif_labels)]),
         Group = "Cryptic")



# Summarize annotated data and filter out any NA
annotated_motif_summary <- leafcutter_star_comparison_annotated_fi %>%
  mutate(junction_id = paste(V1, V2, V3, sep=":")) %>%
  filter(junction_id %in% introns_fnf_sig$junction_id) %>% 
  group_by(V5) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(motif = factor(as.character(V5),
                        levels = names(motif_labels),
                        labels = motif_labels[names(motif_labels)]),
         Group = "Annotated")  


# Combine the summaries
differential_combined_summary <- bind_rows(cryptic_motif_summary, annotated_motif_summary) %>%
  mutate(Group = factor(Group, levels = c("Annotated", "Cryptic")))

# Define full motif levels to enforce consistent x-axis
full_motifs <- c( "GT/AG", "CT/AC", "GC/AG", "CT/GC", "AT/AC", "GT/AT","non-canonical")


# Complete the data so each motif appears in both groups (fill missing with count = 0)
    differential_combined_summary_perct <- differential_combined_summary %>%
      complete(motif = full_motifs, Group, fill = list(count = 0)) %>%
      group_by(Group) %>%
      mutate(total = sum(count),
            perc = count / total * 100) %>%
      ungroup()

# Create the barplot with percentage on the y-axis
oligo_barplot_differential_only <- ggplot(differential_combined_summary_perct, aes(x = motif, y = perc, fill = Group)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = sprintf("%.2f%%", perc)),
            position = position_dodge(width = 0.9),
            vjust = -0.5,
            size = 2.5,
            family = "Helvetica") +
  labs(
    x = "Intron Motif",
    y = "Percentage of splice intron junctions",
    title = "Distribution of Intron Motifs"
  ) +
  scale_y_continuous(
    name = "Percentage of splice intron junctions",
    limits = c(0, 60),
    breaks = seq(0, 60, by = 10),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(values = c("Cryptic" = "grey75", "Annotated" = "grey35")) +
  scale_x_discrete(
    limits = full_motifs
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(linewidth = 0.25),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black", linewidth = 0.25),
    axis.ticks.length.y = unit(-0.1, "cm"),
    axis.title.x = element_markdown(size = 8, family = "Helvetica", margin = margin(t = 5)),
    axis.title.y = element_markdown(size = 8, family = "Helvetica", margin = margin(r = 5)),
    text = element_text(family = "Helvetica"),
    axis.text.y = element_text(color = "black", size = 6),
    axis.text.x = element_text(color = "black", size = 6, margin = margin(t = 5), hjust = 1),
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position.inside = c(0.8, 0.8),
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 6),
    panel.spacing.y = unit(0.5, "cm"),
    plot.title = element_blank()
  )
save(oligo_barplot_differential_only,
     file = "revision/plots/splicing_oligo_barplot_differential_only.rda")





# plot gardener to keep both plot at the same time------------------------------
pdf(file = "revision/plots/splicing_oligo_Barplot.pdf",   # The directory you want to save the file in
    width = 7.5, # The width of the plot in inches
    height = 6)

pageCreate(width = 7.5, height =6 , default.units = "inches", showGuides = FALSE)
load("revision/plots/splicing_oligo_barplot_differential_only.rda")
load("revision/plots/splicing_oligo_barplot_full.rda")

plotGG(oligo_barplot_full, x = 0.4, y = 0.5, width = 7, height = 2.5)
plotGG(oligo_barplot_differential_only, x = 0.4, y = 3.5, width = 7, height = 2.5)


dev.off()

# 3. Explanation for unexpected results with the SNRNP70



