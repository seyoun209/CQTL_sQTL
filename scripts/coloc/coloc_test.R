OAsubtypes <-OAsubtypes[2]
for (subtype in OAsubtypes){
  coloc_results <- list()
  if (subtype == "HipOA") {
    # For HipOA, exclude files that contain KneeHipOA
    all_files <- list.files("./output/coloc/gwas_qtl_ld05_subset", full.names = TRUE)
    subtype_fi <- all_files[grep("_HipOA", all_files, fixed = TRUE)]
    subtype_fi <- subtype_fi[!grepl("KneeHipOA", subtype_fi, fixed = TRUE)]
  } else {
    # For all other subtypes, use the original pattern matching
    subtype_fi <- list.files("./output/coloc/gwas_qtl_ld05_subset", pattern = subtype, full.names = TRUE)
  }
  
  if(length(subtype_fi) == 0 ){
    print(paste0("no-subtype: ",subtype))
    next
  }
  #calcualte fraction cases
  
  test_fractionset <- case_control_sizes[case_control_sizes$OAsubtype == subtype,]
  fraction_cases <- calculate_case_fraction(test_fractionset)
  gwas_df <- fread(subtype_fi[1])
  
  gwas_leadSig <- gwas_df$`CHR:hg38POS` |> unique() |> as.character() 
  gwas_leadSig <- data.frame(gwas_leadSig) |> dplyr::rename("value" ="gwas_leadSig") |>
    dplyr::mutate(
      chr = paste0("chr",str_split(value, ":", simplify = TRUE)[,1]),
      pos = as.double(str_split(value, ":", simplify = TRUE)[,2])
    )
  
  gwas_data <- fread(paste0("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg38_LD/", subtype, "/leads/EUR_", subtype, "_leads_ld_final.csv"))
  
  #---------------------------------------------------------------------------
  # qtl
  
  qtl_df <- fread(subtype_fi[2])
  qtl_leadsig <- unique(qtl_df$lead_var_id ) |> tibble(value =`unique(qtl_df$lead_var_id)`) |> 
    dplyr::select(-c(`unique(qtl_df$lead_var_id)`)) %>%
    dplyr::mutate(
      chr = paste0(str_split(value, ":", simplify = TRUE)[,1]),
      pos = as.double(str_split(value, ":", simplify = TRUE)[,2])
    )
  
  for (chr in gwas_leadSig$chr){
    gwas_testID <- gwas_leadSig[gwas_leadSig[,2] == chr,"value"]
    
    file_path <- paste0("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg38_LD/",
                        subtype,
                        "/summary_stats/",
                        subtype,"_", chr, ".csv")
    
    gwas_summary_stat <- fread(file_path)
    
    if(isEmpty(which(gwas_summary_stat$`CHR:hg38POS` == gwas_testID))){
      next
    }
    
    gwas_summary_stat <- gwas_summary_stat[!is.na(gwas_summary_stat$hg38pos), ]
    
    
    gwas_summary_stat_maf <- gwas_summary_stat %>%
      mutate(MAF = ifelse(EAF < 0.5, EAF, 1 - EAF)) %>%
      filter(MAF != 0)
    
    gwas_data_subset_lds <- gwas_data %>% dplyr::filter(`CHR:hg38POS` %in% gwas_testID) |>
      dplyr::select("ldbuddy_CHR:hg38POS","ldbuddy_ref","ldbuddy_alt") |> dplyr::filter(!is.na(ldbuddy_ref))
    gwas_dataset2 <- inner_join(gwas_data_subset_lds, gwas_summary_stat_maf, by=c("ldbuddy_CHR:hg38POS"="CHR:hg38POS"))
    gwas_dataset2_fi <- gwas_dataset2 |> dplyr::rename("beta"="BETA",
                                                       "snp"="ldbuddy_CHR:hg38POS",
                                                       "MAF"="MAF",
                                                       "pvalues" = "p")
    
    qtl_leadsig_test <- qtl_leadsig[qtl_leadsig$chr == chr,"value"] 
    
    # Load only necessary columns and filter by relevant phe_id values
    pbs_norm_qtl_chr <- fread(paste0("output/nominal_1mb/nominal_pbs/pc5/", chr, ".pbs.cis"), 
                              select = c(1, 8, 9, 11,12,13,14)) 
    colnames(pbs_norm_qtl_chr) <- norm_header[c(1, 8, 9, 11,12,13,14)]
    
    min_pos <- min(gwas_dataset2_fi$hg38pos ) - 500000
    max_pos <- min(gwas_dataset2_fi$hg38pos ) + 500000
    
    pbs_norm_qtl_subset <- pbs_norm_qtl_chr %>%
      dplyr::filter(var_to >= min_pos &
                      var_to <= max_pos) 
    
    pbs_norm_qtl_subset <- pbs_norm_qtl_subset %>%
      dplyr::mutate(var_split = str_split(var_id, ":"),
                    var_allele1 = sapply(var_split, function(x) x[3]),
                    var_allele2 = sapply(var_split, function(x) x[4])) %>%
      dplyr::select(-var_split)
    
    setDT(pbs_norm_qtl_subset)
    setDT(gwas_dataset2_fi)
    # Create column for the gwas-daaset2
    # Generate both orientations of var_id for gwas_dataset2_fi
    gwas_dataset2_fi[, c("var_id1", "var_id2") := create_var_id(chrom, hg38pos, ldbuddy_ref, ldbuddy_alt)]
    unique_pbs_norm_qtl_varID <- unique(pbs_norm_qtl_subset$var_id) |> data.frame() |> 
      dplyr::rename(var_id ="unique.pbs_norm_qtl_subset.var_id.")
    # Perform the matching
    matched_id <- gwas_dataset2_fi[unique_pbs_norm_qtl_varID, on = .(var_id1 = var_id), nomatch = 0] |> dplyr::select(-"var_id2")
    need_to_flip_id <- gwas_dataset2_fi[unique_pbs_norm_qtl_varID, on = .(var_id2 = var_id), nomatch = 0] |> 
      dplyr::select(-"var_id1") |> dplyr::rename("var_id1" ="var_id2")
    
    for(k in qtl_leadsig_test$value){
      
      pheid_pbs <- response_pbs_results |> dplyr::filter(var_id %in% k) |> 
        dplyr::select("phe_id", "SYMBOL")
      
      

      # For non-matching rows, use the regular format
      non_matching <- gwas_dataset2_fi[!var_id1 %in% unique_pbs_norm_qtl_varID$var_id &
                                         !var_id2 %in% unique_pbs_norm_qtl_varID$var_id] |> dplyr::select(-"var_id2")


      # Combine the results
      gwas_dataset2_flip <- rbindlist(list(matched_id, need_to_flip_id,non_matching), fill = TRUE)
      
      
      for(p in pheid_pbs$phe_id){
        
        #pheid_fnf <- fnf_sig_qtl_cond_annot |> dplyr::filter(lead_sqtl %in% qtl_leadsig_test) |> dplyr::select("phe_id", "SYMBOL")
        #ld0_qtl_chr <- ld0_qtl[[chr]] %>% dplyr::filter(leadsnp == k) 
        #ld0_qtl_chr_range <- ld0_qtl_chr %>% dplyr::mutate( chr_ldbuddy = paste0("chr",str_split(ld0_qtl_chr$ldbuddy,":",simplify=TRUE)[,1]),
        #                                                    pos_ldbuddy = as.double(str_split(ld0_qtl_chr$ldbuddy,":",simplify=TRUE)[,2]))
        
        
        
        #-------------------------------------------------------------------------
        # QTL
        
        pbs_qtl_region <- pbs_norm_qtl_chr |> dplyr::filter(phe_id %in% p)
        pbs_qtl_region <- copy(pbs_qtl_region)
        pbs_qtl_region[, snpchr_loc := process_var_id(pbs_qtl_region$var_id)]
        
        pbs_qtl_region_fi <- inner_join(pbs_qtl_region, maf_id_add, by=c("var_id"="var_id"))
        
        qtl_min <- as.numeric(strsplit(k,":")[[1]][2]) - 500000
        qtl_max <- as.numeric(strsplit(k,":")[[1]][2]) + 500000
        
        gwas_dataset2_flip_region <- gwas_dataset2_flip %>%
          dplyr::filter(hg38pos >= qtl_min &
                          hg38pos <= qtl_max) 
        # Remove duplicates based on specific columns, e.g., 'var_id1'
        gwas_dataset2_final_range <- unique(gwas_dataset2_flip_region, by = "var_id1")
        
        gwas_dataset_list <- list(beta = gwas_dataset2_final_range$beta,
                                  snp = gwas_dataset2_final_range$var_id1,
                                  MAF = gwas_dataset2_final_range$MAF,
                                  pvalues= gwas_dataset2_final_range$pvalues,
                                  type="cc",
                                  s = fraction_cases$fraction_cases,
                                  N = fraction_cases$total_samples)
        
        pbs_dataset_list <- list(pvalues = pbs_qtl_region_fi$nom_pval,
                                 beta = pbs_qtl_region_fi$slope,
                                 N = n,
                                 MAF= pbs_qtl_region_fi$MAF,
                                 snp = pbs_qtl_region_fi$var_id,
                                 type = "quant")
        
        
        # coloc run
        coloc_result_name <- paste0(subtype,"_",gwas_testID,"_",k,"_",p)
        print(coloc_result_name)
        coloc_result <- coloc.abf(dataset1 = pbs_dataset_list, dataset2= gwas_dataset_list)
        
        # Store the result in the list
        coloc_results[[coloc_result_name]] <- coloc_result
        
      }
    }
  }
  save(coloc_results, file= paste0("output/coloc/coloc_results_pbs2_",subtype,".rda" ))
}