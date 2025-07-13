library(GenomicRanges)
library(data.table)

load("external_data/encode_rbp/rbp_prep_rbpOnly/pbs_rbp_elip.Rdata") #pbs_rbp_eclip 
load("external_data/encode_rbp/rbp_prep_rbpOnly/fnf_rbp_elip.Rdata") #fnf_rbp_eclip

pbs_rbp_eclip_subset <- pbs_rbp_eclip |> group_by(rbp) %>%
  dplyr::filter(observed > 5) |>
  ungroup()  |>
  arrange(desc(odd_med))

fnf_rbp_eclip_subset <- fnf_rbp_eclip |> group_by(rbp) %>%
  dplyr::filter(observed > 5) |>
  ungroup()  |>
  arrange(desc(odd_med))

combined_data_rbp_eclip_sig <- rbind( pbs_rbp_eclip_subset,fnf_rbp_eclip_subset )
all_sig_enriched_rbps  <- unique(combined_data_rbp_eclip_sig$rbp)

# First, subset the all sig enriched rbps bed file and make grange to find the overlaps with the sQTLs. 

#qtl
pbs_sQtl_bed <- fread("output/Enrichment/rbp/data_prep/significant_pbs_rank0.bed")
colnames(pbs_sQtl_bed) <- c("var_chr", "var_start", "var_end", "var_id", "phe_id", "phe_strd")
fnf_sQtl_bed <- fread("output/Enrichment/rbp/data_prep/significant_fnf_rank0.bed")
colnames(fnf_sQtl_bed) <-  c("var_chr", "var_start", "var_end", "var_id", "phe_id", "phe_strd")

# Convert sQTL data to GRanges objects
pbs_gr <- with(pbs_sQtl_bed, GRanges(var_chr, IRanges(var_start, var_end)))
fnf_gr <- with(fnf_sQtl_bed, GRanges(var_chr, IRanges(var_start, var_end)))

#rbp_bed

rbp_sig.list <- list()
for(i in all_sig_enriched_rbps){
  print(i)
  rbp_bed <- fread(paste0("external_data/encode_rbp/rbp_prep_rbpOnly/",i,".bed"))
  colnames(rbp_bed) <- c("chr","start","end","strand","rbp_nm","cell_line","rep_n","-log2FC","-log10pval")
  rbp_sig.list[[i]] <- rbp_bed
}


# Function to find overlaps and return a data frame
find_overlaps <- function(sqtl_gr, rbp_bed) {
  rbp_gr <- with(rbp_bed, GRanges(chr, IRanges(start, end)))
  overlaps <- findOverlaps(sqtl_gr, rbp_gr)
  data.frame(
    sqtl_index = queryHits(overlaps),
    rbp_index = subjectHits(overlaps)
  )
}


pbs_overlaps <- lapply(rbp_sig.list, function(rbp) find_overlaps(pbs_gr, rbp))
fnf_overlaps <- lapply(rbp_sig.list, function(rbp) find_overlaps(fnf_gr, rbp))

#chr16: 4624826 4690972 
#4683932	4688796

#wwp2
#chr16 69925485 69925980
#mgrn1 <- GRanges(seqnames = "chr16", IRanges(4661789,4661790))

#wwp2 <- GRanges(seqnames = "chr16", IRanges(69921901,69921902))

#find_overlaps_mgrn1 <- find_overlaps(mgrn1, rbp_all.df)
#find_overlaps_wwp2 <- find_overlaps(wwp2, rbp_all.df)

#Now, 
create_overlap_matrix <- function(overlaps, sqtl_bed) {
  rbp_names <- names(overlaps)
  matrix_data <- sapply(rbp_names, function(rbp) {
    overlap_data <- overlaps[[rbp]]
    sqtl_bed$phe_id[overlap_data$sqtl_index]
  })
  
  
  for (i in seq_along(rbp_names)) {
    psi_matrix <- ctl_fnf_ratio |> dplyr::filter(Junction %in% matrix_data[[i]])
  }
  
  psi_matrix
}

pbs_matrix <- create_overlap_matrix(pbs_overlaps, pbs_sQtl_bed)
fnf_matrix <- create_overlap_matrix(fnf_overlaps, fnf_sQtl_bed)


create_overlap_matrix <- function(qtl_rbp_overlaps, sqtl_bed,
                                  ctl_fnf_ratio, vsd_geneExp, response_results, 
                                  target_rbps = c("EFTUD2", "SRSF9", "PRPF8", "SF3B4", "DROSHA"), 
                                  r2_threshold = 0.2) {
  
  # Combine PBS and FNF overlaps
  all_overlaps <- c(qtl_rbp_overlaps)
  all_sqtl_bed <- rbind(sqtl_bed)
  
  # Get variants overlapping with target RBPs
  var_id_counts <- all_overlaps %>%
    keep(names(.) %in% target_rbps) %>%
    imap_dfr(~ tibble(
      rbp = .y,
      var_id = all_sqtl_bed$var_id[.x$sqtl_index],
      phe_id = all_sqtl_bed$phe_id[.x$sqtl_index]
    )) %>%
    group_by(var_id, phe_id) %>%
    summarise(
      rbp_count =  n_distinct(rbp),
      overlapping_rbps = list(unique(rbp)),
      .groups = "drop"
    ) %>%
    arrange(desc(rbp_count))
  
  # Get PSI data
  psi_data <- ctl_fnf_ratio %>%
    pivot_longer(cols = -Junction, names_to = "sample", values_to = "PSI") %>%
    filter(Junction %in% var_id_counts$phe_id)
  
  # Get RBP expression data
  rbp_exp_data <- vsd_geneExp %>%
    filter(SYMBOL %in% target_rbps) %>%
    pivot_longer(cols = -c(ENSG, gene_id, SYMBOL), names_to = "sample", values_to = "expression")
  
  # Combine PSI and RBP expression data
  combined_data <- psi_data %>%
    left_join(rbp_exp_data, by = "sample")
  
  # Calculate R-squared using linear regression for each variant-RBP pair
  correlations <- var_id_counts %>%
    rowwise() %>%
    mutate(
      correlations = list(
        map_dfr(overlapping_rbps, function(rbp) {
          data <- combined_data %>%
            filter(Junction == phe_id, SYMBOL == rbp)
          if(nrow(data) > 0) {
            model <- lm(PSI ~ expression, data = data)
            tibble(
              SYMBOL = rbp,
              r_squared = summary(model)$r.squared,
              p_value = summary(model)$coefficients[2,4]
            )
          } else {
            tibble(SYMBOL = rbp, r_squared = NA, p_value = NA)
          }
        })
      )
    ) %>%
    unnest(correlations) %>%
    filter(r_squared >= r2_threshold, !is.na(r_squared)) |> 
    dplyr::rename(rbp_name = SYMBOL)
  
  # Join with response results
  result <- correlations %>%
    left_join(response_results, by = c("var_id", "phe_id")) %>%
    arrange(desc(r_squared)) 
  
  return(result)
}

fnf_results <- create_overlap_matrix(fnf_overlaps, fnf_sQtl_bed, 
                                ctl_fnf_ratio, vsd_geneExp, response_fnf_results, 
                                target_rbps = c("EFTUD2", "SRSF9", "PRPF8", "SF3B4", "DROSHA","UTP3","NSUN2"), 
                                r2_threshold = 0)

save(fnf_results, file="output/Enrichment/fnf_results_rbp_overlapping.rda")

pbs_results <- create_overlap_matrix(pbs_overlaps, pbs_sQtl_bed, 
                                     ctl_fnf_ratio, vsd_geneExp, response_pbs_results, 
                                     target_rbps = c("EFTUD2", "SRSF9", "PRPF8", "SF3B4", "DROSHA","DDX55","TROVE2","SATU2","U2AF2","TRA2A","CSTF2","SRSF7","SND1","AATF"), 
                                     r2_threshold =0)
save(pbs_results, file="output/Enrichment/pbs_results_rbp_overlapping.rda")


#-------------------------------------------------------------------------------
#Trying to get all the significant RBP from 120 RBPs. 

#Getting combined to one that which rbp overlaps with
create_overlap_matrix <- function(qtl_rbp_overlaps, sqtl_bed,
                                  ctl_fnf_ratio, vsd_geneExp, response_results, 
                                  target_rbps = c("EFTUD2", "SRSF9", "PRPF8", "SF3B4", "DROSHA"), 
                                  r2_threshold = 0.2, p_value_threshold = 0.05) {
  
  # Combine PBS and FNF overlaps
  all_overlaps <- c(qtl_rbp_overlaps)
  all_sqtl_bed <- rbind(sqtl_bed)
  
  # Get variants overlapping with target RBPs
  var_id_counts <- all_overlaps %>%
    keep(names(.) %in% target_rbps) %>%
    imap_dfr(~ tibble(
      rbp = .y,
      var_id = all_sqtl_bed$var_id[.x$sqtl_index],
      phe_id = all_sqtl_bed$phe_id[.x$sqtl_index]
    )) %>%
    group_by(var_id, phe_id) %>%
    summarise(
      rbp_count = n_distinct(rbp),
      overlapping_rbps = list(unique(rbp)),
      .groups = "drop"
    ) %>%
    arrange(desc(rbp_count))
  
  # Get PSI data
  psi_data <- ctl_fnf_ratio %>%
    pivot_longer(cols = -Junction, names_to = "sample", values_to = "PSI") %>%
    filter(Junction %in% var_id_counts$phe_id)
  
  # Get RBP expression data
  rbp_exp_data <- vsd_geneExp %>%
    filter(SYMBOL %in% target_rbps) %>%
    pivot_longer(cols = -c(ENSG, gene_id, SYMBOL), names_to = "sample", values_to = "expression")
  
  # Combine PSI and RBP expression data
  combined_data <- psi_data %>%
    left_join(rbp_exp_data, by = "sample")
  
  # Calculate R-squared using linear regression for each variant-RBP pair
  correlations <- var_id_counts %>%
    rowwise() %>%
    mutate(
      correlations = list(
        map_dfr(overlapping_rbps, function(rbp) {
          data <- combined_data %>%
            filter(Junction == phe_id, SYMBOL == rbp)
          if(nrow(data) > 0) {
            model <- lm(PSI ~ expression, data = data)
            tibble(
              SYMBOL = rbp,
              r_squared = summary(model)$r.squared,
              p_value = summary(model)$coefficients[2,4]
            )
          } else {
            tibble(SYMBOL = rbp, r_squared = NA, p_value = NA)
          }
        })
      )
    ) %>%
    unnest(correlations) %>%
    filter(r_squared >= r2_threshold, p_value <= p_value_threshold, !is.na(r_squared)) %>% 
    dplyr::rename(rbp_name = SYMBOL)
  
  # Join with response results and format overlapping_rbps
  result <- correlations %>%
    left_join(response_results, by = c("var_id", "phe_id")) %>%
    group_by(var_id, phe_id) %>%
    summarise(
      rbp_count = n(),
      overlapping_rbps = paste(rbp_name, collapse = ", "),
      .groups = "drop"
    ) %>%
    left_join(response_results, by = c("var_id", "phe_id")) %>%
    arrange(desc(rbp_count))
  
  return(result)
}




#This will be all of them 
fnf_results_120rbp <- create_overlap_matrix(fnf_overlaps, fnf_sQtl_bed,  ctl_fnf_ratio, vsd_geneExp,  response_fnf_results, 
                                     target_rbps = all_sig_enriched_rbps,
                                     r2_threshold = 0.2, p_value_threshold = 0.05)

save(fnf_results_120rbp, file="output/Enrichment/fnf_results_120rbp_overlapping_.rda")

pbs_results_120rbp <- create_overlap_matrix(pbs_overlaps, pbs_sQtl_bed, 
                                     ctl_fnf_ratio, vsd_geneExp, response_pbs_results, 
                                     target_rbps = all_sig_enriched_rbps, 
                                     r2_threshold = 0.2, p_value_threshold = 0.05)
save(pbs_results_120rbp, file="output/Enrichment/pbs_results_120rbp_overlapping.rda")












psi_summary <- psi_with_meta %>%
  group_by(Junction, Condition) %>%
  summarize(mean_PSI = mean(PSI, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Condition, values_from = mean_PSI) %>%
  mutate(PSI_difference = FNF - CTL)


library(plotgardener)
pageCreate(width = 7.5, height = 2.1, default.units = "inches")

#mgrn1 chr16: 4624826 4690972 

region <- pgParams(
  chrom = "chr16",
  chromstart = 4624826, chromend = 4690972,
  assembly = "hg38",
  range = c(0, 45)
)



#-------------------------------------------------------------------------------



