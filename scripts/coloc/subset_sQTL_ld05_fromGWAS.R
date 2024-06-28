load(file="output/coloc/ld05_leadsnps/ld05_combined_leadsnps.rda") # name is the combined_ld_data
#load(file="output/coloc/ld0_leadsnps/ld0_combined_leadsnps.rda")# name is the combined_ld_data

OAsubtypes <- c("AllOA", "FingerOA", "HandOA", "HipOA", "KneeHipOA", "KneeOA", "THR", "ThumbOA", "TJR", "TKR")

# Output directory for CSV files
output_dir <- "output/coloc/gwas_qtl_ld0_subset"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

GWAS_subset_ld0 <- list()
QTL_subset0 <- list()

for (subtype in OAsubtypes) {
  gwas_data <- fread(paste0("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg38_LD/", subtype, "/leads/EUR_", subtype, "_leads_ld_final.csv"))
  #gwas_data <- fread(paste0("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg38_LD/", subtype, "/leads/ALL_", subtype, "_leads_ld_final.csv"))
  subtype_QTL_ld05sharedGWAS_subset <- data.table()
  subtype_GWAS_sharedQTL_subset <- data.table()
  for (chr in 1:22) {
    gwas_chr_data <- gwas_data %>% filter(chrom == chr)
    combined_ld_qtl <- combined_ld_data[[paste0("chr", chr)]]
    
    #Finding the leadsnp is in the GWAS data LD > 0.5
    gwas_lead_ldbuddy <- gwas_chr_data %>% 
      dplyr::filter(`ldbuddy_CHR:hg38POS` %in% combined_ld_qtl$leadsnp)  |> dplyr::filter(ldbuddy_R2 > 0.5) %>%
      dplyr::select(`CHR:hg38POS`,`ldbuddy_CHR:hg38POS`)
    
  
    qtl_ld05with_gwas <- combined_ld_qtl[which(combined_ld_qtl$leadsnp %in% gwas_lead_ldbuddy$`ldbuddy_CHR:hg38POS`),]
    
    if(nrow(qtl_ld05with_gwas[qtl_ld05with_gwas$ldbuddy %in%  gwas_lead_ldbuddy$`CHR:hg38POS`,] ) == 0){
      next
    }
    
    gwas_overlapswithQTL <- gwas_chr_data[which(gwas_chr_data$`CHR:hg38POS` %in% gwas_lead_ldbuddy$`CHR:hg38POS`),]
    
    
    subtype_QTL_ld05sharedGWAS_subset <- rbind(subtype_QTL_ld05sharedGWAS_subset, qtl_ld05with_gwas)
    subtype_GWAS_sharedQTL_subset <- rbind(subtype_GWAS_sharedQTL_subset, gwas_overlapswithQTL)
    
  }
  # Save results
  # Check if there is data to be saved before writing files
  if (nrow(subtype_QTL_ld05sharedGWAS_subset) > 0) {
    fwrite(subtype_QTL_ld05sharedGWAS_subset, paste0("output/coloc/gwas_qtl_ld05_subset/QTL_ld05_withGWAS_", subtype, ".csv"))
  }
  
  if (nrow(subtype_GWAS_sharedQTL_subset) > 0) {
    fwrite(subtype_GWAS_sharedQTL_subset, paste0("output/coloc/gwas_qtl_ld05_subset/GWAS_shared_withQTL_", subtype, ".csv"))
  }
  
  # Store results in lists
  GWAS_subset_ld0[[subtype]] <- subtype_GWAS_sharedQTL_subset
  QTL_subset0[[subtype]] <- subtype_QTL_ld05sharedGWAS_subset
  
  cat("Processed:", subtype, "\n")
}