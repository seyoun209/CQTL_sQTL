#RBPmap genomic coordinate 

# Function to convert rMATS coordinates to genomic coordinate format for RBPmap
convert_rmats_to_rbpmap <- function(rmats_data, event_type) {
  rbpmap_coords <- character(nrow(rmats_data))
  
  for(i in 1:nrow(rmats_data)) {
    chr <- rmats_data$Chr[i]
    strand <- rmats_data$Strand[i]
    
    if(event_type == "SE") {
      # For skipped exons: focus on upstream and downstream splice sites
      # These are the critical regions where splicing factors bind
      target_coords <- strsplit(rmats_data$exon_target[i], "-")[[1]]
      upstream_coords <- strsplit(rmats_data$exon_upstream[i], "-")[[1]]
      downstream_coords <- strsplit(rmats_data$exon_downstream[i], "-")[[1]]
      
      # Generate 3 coordinates - upstream 3' splice site, target exon with flanking regions, and downstream 5' splice site
      upstream_3ss <- paste0(chr, ":", as.numeric(upstream_coords[2]) - 20, "-", as.numeric(upstream_coords[2]) + 40, ":", strand)
      target_region <- paste0(chr, ":", as.numeric(target_coords[1]) - 50, "-", as.numeric(target_coords[2]) + 50, ":", strand)
      downstream_5ss <- paste0(chr, ":", as.numeric(downstream_coords[1]) - 40, "-", as.numeric(downstream_coords[1]) + 20, ":", strand)
      
      rbpmap_coords[i] <- paste(upstream_3ss, target_region, downstream_5ss, sep = ";")
      
    } else if(event_type == "RI") {
      # For retained introns: focus on both splice sites
      ri_coords <- strsplit(rmats_data$exon_ir[i], "-")[[1]]
      up_coords <- strsplit(rmats_data$exon_upstream[i], "-")[[1]]
      down_coords <- strsplit(rmats_data$exon_downstream[i], "-")[[1]]
      
      # 5' splice site and 3' splice site
      ss5 <- paste0(chr, ":", as.numeric(ri_coords[1]) - 40, "-", as.numeric(ri_coords[1]) + 60, ":", strand)
      ss3 <- paste0(chr, ":", as.numeric(ri_coords[2]) - 60, "-", as.numeric(ri_coords[2]) + 40, ":", strand)
      
      rbpmap_coords[i] <- paste(ss5, ss3, sep = ";")
      
    } else if(event_type == "A3SS") {
      # Alternative 3' splice sites
      long_coords <- strsplit(rmats_data$exon_long[i], "-")[[1]]
      short_coords <- strsplit(rmats_data$exon_short[i], "-")[[1]]
      
      # Focus on both competing 3' splice sites
      short_3ss <- paste0(chr, ":", as.numeric(short_coords[1]) - 60, "-", as.numeric(short_coords[1]) + 40, ":", strand)
      long_3ss <- paste0(chr, ":", as.numeric(long_coords[1]) - 60, "-", as.numeric(long_coords[1]) + 40, ":", strand)
      
      rbpmap_coords[i] <- paste(short_3ss, long_3ss, sep = ";")
      
    } else if(event_type == "A5SS") {
      # Alternative 5' splice sites
      long_coords <- strsplit(rmats_data$exon_long[i], "-")[[1]]
      short_coords <- strsplit(rmats_data$exon_short[i], "-")[[1]]
      
      # Focus on both competing 5' splice sites
      short_5ss <- paste0(chr, ":", as.numeric(short_coords[2]) - 40, "-", as.numeric(short_coords[2]) + 60, ":", strand)
      long_5ss <- paste0(chr, ":", as.numeric(long_coords[2]) - 40, "-", as.numeric(long_coords[2]) + 60, ":", strand)
      
      rbpmap_coords[i] <- paste(short_5ss, long_5ss, sep = ";")
      
    } else if(event_type == "MXE") {
      # Mutually exclusive exons
      exon1_coords <- strsplit(rmats_data$exon_1[i], "-")[[1]]
      exon2_coords <- strsplit(rmats_data$exon_2[i], "-")[[1]]
      up_coords <- strsplit(rmats_data$exon_upstream[i], "-")[[1]]
      down_coords <- strsplit(rmats_data$exon_downstream[i], "-")[[1]]
      
      # Focus on all four relevant splice sites
      up_3ss <- paste0(chr, ":", as.numeric(up_coords[2]) - 20, "-", as.numeric(up_coords[2]) + 40, ":", strand)
      exon1_5ss <- paste0(chr, ":", as.numeric(exon1_coords[2]) - 40, "-", as.numeric(exon1_coords[2]) + 60, ":", strand)
      exon2_3ss <- paste0(chr, ":", as.numeric(exon2_coords[1]) - 60, "-", as.numeric(exon2_coords[1]) + 40, ":", strand)
      down_5ss <- paste0(chr, ":", as.numeric(down_coords[1]) - 40, "-", as.numeric(down_coords[1]) + 20, ":", strand)
      
      rbpmap_coords[i] <- paste(up_3ss, exon1_5ss, exon2_3ss, down_5ss, sep = ";")
    }
  }
  
  # Expand the semicolon-separated coordinates into individual entries
  expanded_coords <- unlist(strsplit(rbpmap_coords, ";"))
  
  return(expanded_coords)
}

# Process all event types
se_coords <- convert_rmats_to_rbpmap(summary(sig_kd, 'SE'), "SE")
ri_coords <- convert_rmats_to_rbpmap(summary(sig_kd, 'RI'), "RI")
a3ss_coords <- convert_rmats_to_rbpmap(summary(sig_kd, 'A3SS'), "A3SS")
a5ss_coords <- convert_rmats_to_rbpmap(summary(sig_kd, 'A5SS'), "A5SS")
mxe_coords <- convert_rmats_to_rbpmap(summary(sig_kd, 'MXE'), "MXE")

# Combine all coordinates
all_coords <- c(se_coords, ri_coords, a3ss_coords, a5ss_coords, mxe_coords)

# Write coordinates to a file for RBPmap
write.table(all_coords, file="/work/users/s/e/seyoun/CQTL_sQTL/output/rmats_edited/rbpmap_coordinates.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)


#------------------------------------------------------------------------------
#rbp output cleaning up

path_rbp_re <- "/work/users/s/e/seyoun/CQTL_sQTL/output/rmats_edited/rbpmap/"
rbp_edited <- readLines(paste0(path_rbp_re,"all_prediction_snrnp70_high_conserved.txt"))

# Now map these coordinates to your original gene list
# A better approach to create the mapping between coordinates and genes
# Start with empty dataframes for each event type
se_mapping <- data.frame()
ri_mapping <- data.frame()
a3ss_mapping <- data.frame()
a5ss_mapping <- data.frame()
mxe_mapping <- data.frame()

# For SE events
if(length(se_coords) > 0) {
  # Get the original SE data
  se_data <- summary(sig_kd, 'SE')
  
  # Calculate how many coordinates were generated per SE event
  coords_per_event <- length(se_coords) / nrow(se_data)
  
  # Create mapping
  se_mapping <- data.frame(
    coord = se_coords,
    gene = rep(se_data$geneSymbol, each = coords_per_event),
    event_type = "SE"
  )
}

# For RI events
if(length(ri_coords) > 0) {
  ri_data <- summary(sig_kd, 'RI')
  coords_per_event <- length(ri_coords) / nrow(ri_data)
  ri_mapping <- data.frame(
    coord = ri_coords,
    gene = rep(ri_data$geneSymbol, each = coords_per_event),
    event_type = "RI"
  )
}

# For A3SS events
if(length(a3ss_coords) > 0) {
  a3ss_data <- summary(sig_kd, 'A3SS')
  coords_per_event <- length(a3ss_coords) / nrow(a3ss_data)
  a3ss_mapping <- data.frame(
    coord = a3ss_coords,
    gene = rep(a3ss_data$geneSymbol, each = coords_per_event),
    event_type = "A3SS"
  )
}

# For A5SS events
if(length(a5ss_coords) > 0) {
  a5ss_data <- summary(sig_kd, 'A5SS')
  coords_per_event <- length(a5ss_coords) / nrow(a5ss_data)
  a5ss_mapping <- data.frame(
    coord = a5ss_coords,
    gene = rep(a5ss_data$geneSymbol, each = coords_per_event),
    event_type = "A5SS"
  )
}

# For MXE events
if(length(mxe_coords) > 0) {
  mxe_data <- summary(sig_kd, 'MXE')
  coords_per_event <- length(mxe_coords) / nrow(mxe_data)
  mxe_mapping <- data.frame(
    coord = mxe_coords,
    gene = rep(mxe_data$geneSymbol, each = coords_per_event),
    event_type = "MXE"
  )
}

# Combine all mappings
all_coords_info <- rbind(√, ri_mapping, a3ss_mapping, a5ss_mapping, mxe_mapping)

#-------------------------------------------------------------------------------


# Initialize variables
coord_hits <- character(0)
binding_details <- list()  # Store detailed binding information
current_coord <- ""
current_details <- list()

# Process the file to extract coordinates with motifs and their details
for(i in 1:length(rbp_edited)) {
  line <- rbp_edited[i]
  
  # Check if the line contains a genomic coordinate
  if(grepl("^chr", line) && grepl("-", line) && !grepl("=", line)) {
    # Save previous coordinate if it had motifs
    if(length(current_details) > 0) {
      coord_hits <- c(coord_hits, current_coord)
      binding_details[[current_coord]] <- current_details
    }
    
    # Start processing a new coordinate
    current_coord <- line
    current_details <- list()
  }
  
  # Check if the line contains motif information (but not the header line)
  if(grepl("\\d+\\s+chr", line) && grepl("rwucaag", line)) {
    # Parse the line to extract all details
    # Expected format: position, coordinate, motif, k-mer, z-score, p-value
    parts <- strsplit(gsub("\\s+", " ", trimws(line)), " ")[[1]]
    
    if(length(parts) >= 6) {
      # Find the indices that likely contain the values we want
      pos_idx <- which(grepl("^\\d+$", parts))[1]
      if(!is.na(pos_idx) && pos_idx <= length(parts) - 5) {
        detail <- list(
          position = parts[pos_idx],
          genomic_coord = parts[pos_idx + 1],
          motif = parts[pos_idx + 2],
          k_mer = parts[pos_idx + 3],
          z_score = as.numeric(parts[pos_idx + 4]),
          p_value = as.numeric(parts[pos_idx + 5])
        )
        
        # Only include if Z-score is significant
        if(detail$z_score > 1.6) {
          current_details[[length(current_details) + 1]] <- detail
        }
      }
    }
  }
}

# Check the last coordinate
if(length(current_details) > 0) {
  coord_hits <- c(coord_hits, current_coord)
  binding_details[[current_coord]] <- current_details
}

# Create a comprehensive results dataframe
snrnp70_binding_results <- data.frame(
  gene = character(0),
  event_type = character(0),
  coordinate = character(0),
  binding_coord = character(0),
  position = character(0),
  genomic_coord = character(0),
  motif = character(0),
  k_mer = character(0),
  z_score = numeric(0),
  p_value = numeric(0),
  stringsAsFactors = FALSE
)

# For each gene/event in your original data, check if it has SNRNP70 binding
for(i in 1:nrow(all_coords_info)) {
  gene_coord <- all_coords_info$coord[i]
  gene_symbol <- all_coords_info$gene[i]
  event_type <- all_coords_info$event_type[i]
  
  # Check if this coordinate overlaps with any SNRNP70 binding site
  has_binding <- FALSE
  
  for(hit_coord in coord_hits) {
    if(coords_overlap(hit_coord, gene_coord)) {
      # This gene has binding site(s)
      has_binding <- TRUE
      
      # Add details for each binding site to the results
      for(detail in binding_details[[hit_coord]]) {
        snrnp70_binding_results <- rbind(snrnp70_binding_results, data.frame(
          gene = gene_symbol,
          event_type = event_type,
          coordinate = gene_coord,
          binding_coord = hit_coord,
          position = detail$position,
          genomic_coord = detail$genomic_coord,
          motif = detail$motif,
          k_mer = detail$k_mer,
          z_score = detail$z_score,
          p_value = detail$p_value,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # If no binding, still add a row but with NA values for binding details
  if(!has_binding) {
    snrnp70_binding_results <- rbind(snrnp70_binding_results, data.frame(
      gene = gene_symbol,
      event_type = event_type,
      coordinate = gene_coord,
      binding_coord = NA,
      position = NA,
      genomic_coord = NA,
      motif = NA,
      k_mer = NA,
      z_score = NA,
      p_value = NA,
      stringsAsFactors = FALSE
    ))
  }
}


snrnp70_binding_re_noNA <- snrnp70_binding_results |> filter(!is.na(p_value)) |> 
  distinct(coordinate,binding_coord,genomic_coord,motif,k_mer,z_score,p_value,.keep_all = TRUE)




#-------------------------------------------------------------------------------
#Boxplot visualization for the SNRNP70 binding site near the intronjunction

# Get unique genes with SNRNP70 binding from your results
snrnp70_binding_genes <- unique(snrnp70_binding_re_noNA$gene)

# Add SNRNP70 binding status to your plot_data
kd_comparedtoOA_plot_data_subset <- kd_comparedtoOA_plot_data |> filter(!geneSymbol_sig %in% snrnp70_binding_genes)

fnf_comparedtoKD_plot_data_subset <- fnf_comparedtoKD_plot_data |> filter(geneSymbol_sig %in% snrnp70_binding_genes)
