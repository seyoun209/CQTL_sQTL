load(file="output/coloc/ld05_leadsnps/ld05_combined_leadsnps.rda") # name is the combined_ld_data
#load(file="output/coloc/ld0_leadsnps/ld0_combined_leadsnps.rda")# name is the combined_ld_data

load(file="output/01.qtltools_re/conditional_filterout_snps.txt") # conditional_filterout_snps

# Function to read and process BED files
read_bed_file <- function(file_path) {
  fread(file_path, col.names = c("chr", "start", "var_id", "ref", "alt", "rsID"))
}
maf_all <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/freq/cqtl.frq")
maf_subset <-  maf_all %>% dplyr::select(c("SNP","MAF","A1")) %>% 
  dplyr::rename("var_id" =SNP, "minor_allele" =A1) 
# Get list of BED files
bed_files <- list.files("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/rsid",pattern = "_rsid.bed",full.names = TRUE)
all_bed_data <- rbindlist(lapply(bed_files, read_bed_file))
all_bed_dataID <- all_bed_data |> dplyr::select(c("var_id","rsID"))

#merge rsIDs
maf_id_add <- left_join(maf_subset , all_bed_dataID,by="var_id") |> dplyr::select(c('var_id', 'rsID',"MAF"))
save(maf_id_add, file="output/geno/maf_withRsID")

OAsubtypes <- c("AllOA", "FingerOA", "HandOA", "HipOA", "KneeHipOA", "KneeOA", "THR", "ThumbOA", "TJR", "TKR")

# Output directory for CSV files
output_dir <- "output/coloc/gwas_qtl_ld05_subset"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

GWAS_subset_ld0 <- list()
QTL_subset0 <- list()
#OAsubtypes <- OAsubtypes[5]
for (subtype in OAsubtypes) {
  print(subtype)
  gwas_data <- fread(paste0("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg38_LD/", subtype, "/leads/EUR_", subtype, "_leads_ld_final.csv"))
  #gwas_data <- fread(paste0("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg38_LD/", subtype, "/leads/ALL_", subtype, "_leads_ld_final.csv"))
  subtype_QTL_ld05sharedGWAS_subset <- data.table()
  subtype_GWAS_sharedQTL_subset <- data.table()
  for (chr in unique(gwas_data$chrom)) {
    print(paste0("chr",chr))
    gwas_chr_data <- gwas_data %>% filter(chrom == chr)
    #combined_ld_qtl <- combined_ld_data[[paste0("chr", chr)]]
    ldfiles <- list.files(paste0("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld05/","chr",chr), pattern = ".ld")
    snpslist <- gsub(".ld", "",ldfiles) |> data.frame()
    colnames(snpslist)[1] <- c("var_id")
    snplist_rsid <- left_join(snpslist, maf_id_add,by="var_id")
    snps_lead <- snplist_rsid[!snplist_rsid$var_id  %in% conditional_filterout_snps,] 

    #Finding the leadsnp is in the GWAS data LD > 0.5
    gwas_lead_ldbuddy <- gwas_chr_data %>% 
      dplyr::filter(ldbuddy_rsID %in% snps_lead$rsID)  |> dplyr::filter(ldbuddy_R2 >= 0.5) %>%
      dplyr::select(rsID ,ldbuddy_rsID, `CHR:hg38POS`, `ldbuddy_CHR:hg38POS` )
    
    snps_leadqtl_overlaps <- snps_lead |> dplyr::filter(rsID %in% gwas_lead_ldbuddy$ldbuddy_rsID) |> dplyr::select(c("var_id","rsID"))
    if (is.na(snps_leadqtl_overlaps[1,1])){next}
    for(j in c(1:nrow(snps_leadqtl_overlaps))){
      
      qtl_ld0with_gwas <- fread(paste0("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld/","chr",chr,"/",snps_leadqtl_overlaps$var_id[j],".ld")) |>
        dplyr::select(c("SNP_B" ,"R2")) |>
        dplyr::mutate(lead_rsID = snps_leadqtl_overlaps$rsID[j],
                      lead_var_id=snps_leadqtl_overlaps$var_id[j] ) |> dplyr::rename("var_id" ="SNP_B")
      qtl_ld_0_overlapswithGWAS <- left_join(qtl_ld0with_gwas,maf_id_add,by="var_id") |>  dplyr::rename("ld_buddy_rsID" ="rsID")
      
      subtype_QTL_ld05sharedGWAS_subset <- rbind(subtype_QTL_ld05sharedGWAS_subset, qtl_ld_0_overlapswithGWAS)
    }
    
    
    gwas_overlapswithQTL <- gwas_chr_data |> dplyr::filter(`CHR:hg38POS` %in% gwas_lead_ldbuddy$`CHR:hg38POS`)
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

save(GWAS_subset_ld0, file ="output/coloc/gwas_qtl_ld05_subset/GWAS_subset_qtl")
save(QTL_subset0, file ="output/coloc/gwas_qtl_ld05_subset/QTL_ld05_subset")
