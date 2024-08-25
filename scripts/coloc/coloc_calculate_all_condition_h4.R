source("scripts/utils/utils.R")
library(coloc)
#PBS----------------------------------------------------------------------------
response_pbs_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_pbs_resInv.rds")
response_fnf_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_fnf_resInv.rds")

pbs_norm_qtl_pc5_list <- readRDS("output/nominal_1mb/nominal_pbs/pbs_norm_qtl_pc5_list.rds")
fnf_norm_qtl_pc4_list <- readRDS("output/nominal_1mb/nominal_fnf/fnf_norm_qtl_pc4_list.rds")

load("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/freq/maf_id_add.rds") #maf_id_add  is the name
maf_subset_rsID <- maf_id_add |> dplyr::select(-c("minor_allele","MAF"))

case_control_sizes <- read_csv("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg38_LD/Case_Control_sampleSizes.csv")

n <- 101
OAsubtypes <- c("AllOA", "FingerOA", "HandOA", "HipOA", "KneeHipOA", "KneeOA", "THR", "ThumbOA", "TJR", "TKR")
coloc_results.df <- data.frame()
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
  #calculate fraction cases
  
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
    gwas_testID <- gwas_leadSig[gwas_leadSig$chr == chr,"value"]
    gwas_testPos <- gwas_leadSig[gwas_leadSig$chr == chr,"pos"]
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
      dplyr::select("ldbuddy_CHR:hg38POS","ldbuddy_ref","ldbuddy_alt","rsID") |> dplyr::filter(!is.na(ldbuddy_ref))
    gwas_dataset2 <- inner_join(gwas_data_subset_lds, gwas_summary_stat_maf, by=c("ldbuddy_CHR:hg38POS"="CHR:hg38POS"))
    gwas_dataset2_fi <- gwas_dataset2 |> dplyr::rename("beta"="BETA",
                                                       "snp"="ldbuddy_CHR:hg38POS",
                                                       "MAF"="MAF",
                                                       "pvalues" = "p")
    
    qtl_leadsig_test <- qtl_leadsig[qtl_leadsig$chr == chr,"value"] 
    qtl_leadsig_test_pos <- qtl_leadsig[qtl_leadsig$chr == chr,"pos"]
    gwas_leadsig_test <- gwas_data  |>  dplyr::filter(ldbuddy_rsID %in% gwas_data_subset_lds$rsID[1]) %>%
      dplyr::slice(1)
    
    # Load only necessary columns and filter by relevant phe_id values
    pbs_norm_qtl_chr <- pbs_norm_qtl_pc5_list[[chr]]
    fnf_norm_qtl_chr <- fnf_norm_qtl_pc4_list[[chr]]
    
    pbs_norm_qtl_alleles <- pbs_norm_qtl_chr %>%
      dplyr::mutate(var_split = str_split(var_id, ":"),
                    var_allele1 = sapply(var_split, function(x) x[3]),
                    var_allele2 = sapply(var_split, function(x) x[4])) %>%
      dplyr::select(-var_split)
    
    fnf_norm_qtl_alleles <- fnf_norm_qtl_chr %>%
      dplyr::mutate(var_split = str_split(var_id, ":"),
                    var_allele1 = sapply(var_split, function(x) x[3]),
                    var_allele2 = sapply(var_split, function(x) x[4])) %>%
      dplyr::select(-var_split)
    
    
    setDT(pbs_norm_qtl_alleles)
    setDT(fnf_norm_qtl_alleles)
    setDT(gwas_dataset2_fi)
    # Create column for the gwas-data-set2
    # Generate both orientations of var_id for gwas-data-set2 fi
    gwas_dataset2_fi[, c("var_id1", "var_id2") := create_var_id(chrom, hg38pos, ldbuddy_ref, ldbuddy_alt)]
    
    unique_pbs_norm_qtl_varID <- unique(pbs_norm_qtl_alleles$var_id) |> data.frame() |> 
      dplyr::rename(var_id ="unique.pbs_norm_qtl_alleles.var_id.")
    
    unique_fnf_norm_qtl_varID <- unique(fnf_norm_qtl_alleles$var_id) |> data.frame() |> 
      dplyr::rename(var_id ="unique.fnf_norm_qtl_alleles.var_id.")
    # Perform the matching
    
    matched_id <- gwas_dataset2_fi[unique_pbs_norm_qtl_varID, on = .(var_id1 = var_id), nomatch = 0] |> dplyr::select(-"var_id2")
    need_to_flip_id <- gwas_dataset2_fi[unique_pbs_norm_qtl_varID, on = .(var_id2 = var_id), nomatch = 0] |> 
      dplyr::select(-"var_id1") |> dplyr::rename("var_id1" ="var_id2")
    
    for(k in qtl_leadsig_test$value){
      pheid_pbs <- response_pbs_results |> dplyr::filter(var_id %in% k) |> 
        dplyr::select("phe_id","SYMBOL","rsID","rank","var_id","minor_allele","MAF","interaction_pval")
      
      qtl_pos <- as.numeric(str_split(k,":",simplify =TRUE)[,2])
      
      
      min_pos_gwas <- gwas_testPos - 250000
      max_pos_gwas <- gwas_testPos + 250000
      
      min_pos_qtl <- qtl_pos - 250000
      max_pos_qtl <- qtl_pos + 250000
      
      min_region <- min(min_pos_qtl,min_pos_gwas)
      max_region <- max(max_pos_qtl,max_pos_gwas)
      
      
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
        
         pheid_pbs_filtered <- pheid_pbs |> dplyr::filter(phe_id == p)
  
  if(nrow(pheid_pbs_filtered) == 0) {
    print(paste("No PBS data for phenotype:", p))
    next  # Skip to the next iteration of the loop
  }
        
        
        #-------------------------------------------------------------------------
        # QTL
        
        pbs_qtl_region <- pbs_norm_qtl_alleles |> dplyr::filter(phe_id %in% pheid_pbs_filtered$phe_id)
        pbs_qtl_region <- copy(pbs_qtl_region)
        pbs_qtl_region[, snpchr_loc := process_var_id(pbs_qtl_region$var_id)]
        
        pbs_qtl_region_fi <- inner_join(pbs_qtl_region, maf_id_add, by=c("var_id"="var_id"))
        
        
        pbs_qtl_region_subset <- pbs_qtl_region_fi %>%
          dplyr::filter(var_from >= min_region &
                          var_from <= max_region) 
        

        
        
        # FNF analysis (using the same phenotype ID as PBS)
        fnf_qtl_region <- fnf_norm_qtl_alleles |> dplyr::filter(phe_id == p)
        fnf_qtl_region <- copy(fnf_qtl_region)
        fnf_qtl_region[, snpchr_loc := process_var_id(fnf_qtl_region$var_id)]
        
        fnf_qtl_region_fi <- inner_join(fnf_qtl_region, maf_id_add, by=c("var_id"="var_id"))
        
        fnf_qtl_region_subset <- fnf_qtl_region_fi %>%
          dplyr::filter(var_from >= min_region & var_from <= max_region) 
        
        
        
        # Gwas analysis
        
        gwas_dataset2_flip_region <- gwas_dataset2_flip %>%
          dplyr::filter(hg38pos >= min_region &
                          hg38pos <= max_region) 
        
        # Remove duplicates based on specific columns, e.g., 'var_id1'
        gwas_dataset2_final_range <- unique(gwas_dataset2_flip_region, by = "var_id1")
        
        gwas_dataset_list <- list(beta = gwas_dataset2_final_range$beta,
                                  snp = gwas_dataset2_final_range$var_id1,
                                  MAF = gwas_dataset2_final_range$MAF,
                                  pvalues= gwas_dataset2_final_range$pvalues,
                                  type="cc",
                                  s = fraction_cases$fraction_cases,
                                  N = fraction_cases$total_samples)
        
        pbs_dataset_list <- list(pvalues = pbs_qtl_region_subset$nom_pval,
                                 beta = pbs_qtl_region_subset$slope,
                                 N = n,
                                 MAF= pbs_qtl_region_subset$MAF,
                                 snp = pbs_qtl_region_subset$var_id,
                                 type = "quant")
        
        fnf_dataset_list <- list(pvalues = fnf_qtl_region_subset$nom_pval,
                                 beta = fnf_qtl_region_subset$slope,
                                 N = n,
                                 MAF= fnf_qtl_region_subset$MAF,
                                 snp = fnf_qtl_region_subset$var_id,
                                 type = "quant")
        
        
        # coloc run
        coloc_result_name <- paste0(subtype,"_",gwas_testID,"_",k,"_",p)
        print(coloc_result_name)
        
        coloc_result <- coloc.abf(dataset1 = pbs_dataset_list, dataset2 = gwas_dataset_list)
        coloc_result_fnf <- coloc.abf(dataset1 = fnf_dataset_list, dataset2 = gwas_dataset_list)
        
        coloc_summary <- c(coloc_result$summary[["PP.H3.abf"]], coloc_result$summary[["PP.H4.abf"]])
        coloc_summary_fnf <- c(coloc_result_fnf$summary[["PP.H3.abf"]], coloc_result_fnf$summary[["PP.H4.abf"]])
        
        
        # Before creating coloc_results_table, add these print statements:
        cat("Debugging information:\n")
        cat("subtype:", subtype, "\n")
        cat("phe_id:", pheid_pbs_filtered$phe_id, "\n")
        cat("gene:", pheid_pbs_filtered$SYMBOL, "\n")
        cat("GWAS_rsID:", gwas_leadsig_test$rsID, "\n")
        cat("GWAS_pos:", gwas_leadsig_test$`CHR:hg38POS`, "\n")
        cat("GWAS_ref:", gwas_leadsig_test$ldbuddy_ref, "\n")
        cat("GWAS_alt:", gwas_leadsig_test$ldbuddy_alt, "\n")
        cat("sQTL_rsID:", pheid_pbs_filtered$rsID, "\n")
        cat("sQTL_rank:", pheid_pbs_filtered$rank, "\n")
        cat("sQTL_var_id:", pheid_pbs_filtered$var_id, "\n")
        cat("sQTL_minor_allele:", pheid_pbs_filtered$minor_allele, "\n")
        cat("sQTL_MAF:", pheid_pbs_filtered$MAF, "\n")
        cat("sQTL_interaction_pval:", pheid_pbs_filtered$interaction_pval, "\n")
        cat("PBS_coloc_H3:", coloc_summary[1], "\n")
        cat("PBS_coloc_H4:", coloc_summary[2], "\n")
        cat("FNF_coloc_H3:", coloc_summary_fnf[1], "\n")
        cat("FNF_coloc_H4:", coloc_summary_fnf[2], "\n")
        
        
        coloc_results_table <- data.frame(
          subtype = subtype,
          phe_id = pheid_pbs_filtered$phe_id,
          gene = pheid_pbs_filtered$SYMBOL,
          GWAS_rsID = gwas_leadsig_test$rsID,
          GWAS_pos = gwas_leadsig_test$`CHR:hg38POS`,
          GWAS_ref = gwas_leadsig_test$ldbuddy_ref,
          GWAS_alt = gwas_leadsig_test$ldbuddy_alt,
          sQTL_rsID = pheid_pbs_filtered$rsID,
          sQTL_rank = pheid_pbs_filtered$rank,
          sQTL_var_id = pheid_pbs_filtered$var_id,
          sQTL_minor_allele = pheid_pbs_filtered$minor_allele,
          sQTL_MAF = pheid_pbs_filtered$MAF,
          sQTL_interaction_pval = pheid_pbs_filtered$interaction_pval,
          PBS_coloc_H3 = coloc_summary[1],
          PBS_coloc_H4 = coloc_summary[2],
          FNF_coloc_H3 = coloc_summary_fnf[1],
          FNF_coloc_H4 = coloc_summary_fnf[2]
        )
        
        
        coloc_results.df <- rbind(coloc_results.df,coloc_results_table)
        
        # Store the result in the list
        coloc_results[[coloc_result_name]] <- coloc_result
        
      }
    }
  }
  
  save(coloc_results, file= paste0("output/coloc/coloc_results_pbs_prob_calc_",subtype,".rda" ))
  
}

colnames(coloc_results.df) <- c("subtype","phe_id","gene","GWAS_rsID","GWAS_pos","GWAS_ref","GWAS_alt",
                                "sQTL_rsID","sQTL_rank","sQTL_var_id","sQTL_minor_allele","sQTL_MAF","sQTL_interaction_pval",
                                "PBS_coloc_H3","PBS_coloc_H4","FNF_coloc_H3","FNF_coloc_H4")
write.table(coloc_results.df, file="output/coloc/coloc_result_table_pbs_prob_calc_.txt",sep="\t",quote=F,col.names = TRUE, row.names = FALSE)






#-------------------------------------------------------------------------------
#-FNF---------------------------------------------------------------------------


n <- 101
OAsubtypes <- c("AllOA", "FingerOA", "HandOA", "HipOA", "KneeHipOA", "KneeOA", "THR", "ThumbOA", "TJR", "TKR")
coloc_results.df <- data.frame()
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
  #calculate fraction cases
  
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
# now gwas part. 
      
  for (chr in gwas_leadSig$chr){
    gwas_testID <- gwas_leadSig[gwas_leadSig$chr == chr,"value"]
    gwas_testPos <- gwas_leadSig[gwas_leadSig$chr == chr,"pos"]
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
      dplyr::select("ldbuddy_CHR:hg38POS","ldbuddy_ref","ldbuddy_alt","rsID") |> dplyr::filter(!is.na(ldbuddy_ref))
    gwas_dataset2 <- inner_join(gwas_data_subset_lds, gwas_summary_stat_maf, by=c("ldbuddy_CHR:hg38POS"="CHR:hg38POS"))
    gwas_dataset2_fi <- gwas_dataset2 |> dplyr::rename("beta"="BETA",
                                                       "snp"="ldbuddy_CHR:hg38POS",
                                                       "MAF"="MAF",
                                                       "pvalues" = "p")
    
    
    
    qtl_leadsig_test <- qtl_leadsig[qtl_leadsig$chr == chr,"value"] 
    qtl_leadsig_test_pos <- qtl_leadsig[qtl_leadsig$chr == chr,"pos"]
    gwas_leadsig_test <- gwas_data  |>  dplyr::filter(ldbuddy_rsID %in% gwas_data_subset_lds$rsID[1]) %>%
      dplyr::slice(1)
    
    # Load QTL data
    pbs_norm_qtl_chr <- pbs_norm_qtl_pc5_list[[chr]]
    fnf_norm_qtl_chr <- fnf_norm_qtl_pc4_list[[chr]]
    
    pbs_norm_qtl_alleles <- pbs_norm_qtl_chr %>%
      dplyr::mutate(var_split = str_split(var_id, ":"),
                    var_allele1 = sapply(var_split, function(x) x[3]),
                    var_allele2 = sapply(var_split, function(x) x[4])) %>%
      dplyr::select(-var_split)
    
    fnf_norm_qtl_alleles <- fnf_norm_qtl_chr %>%
      dplyr::mutate(var_split = str_split(var_id, ":"),
                    var_allele1 = sapply(var_split, function(x) x[3]),
                    var_allele2 = sapply(var_split, function(x) x[4])) %>%
      dplyr::select(-var_split)
    
    
    setDT(pbs_norm_qtl_alleles)
    setDT(fnf_norm_qtl_alleles)
    setDT(gwas_dataset2_fi)
    # Create column for the gwas-data-set2
    # Generate both orientations of var_id for gwas-data-set2 fi
    gwas_dataset2_fi[, c("var_id1", "var_id2") := create_var_id(chrom, hg38pos, ldbuddy_ref, ldbuddy_alt)]
    
    unique_pbs_norm_qtl_varID <- unique(pbs_norm_qtl_alleles$var_id) |> data.frame() |> 
      dplyr::rename(var_id ="unique.pbs_norm_qtl_alleles.var_id.")
    
    unique_fnf_norm_qtl_varID <- unique(fnf_norm_qtl_alleles$var_id) |> data.frame() |> 
      dplyr::rename(var_id ="unique.fnf_norm_qtl_alleles.var_id.")
    # Perform the matching
    
    matched_id <- gwas_dataset2_fi[unique_fnf_norm_qtl_varID, on = .(var_id1 = var_id), nomatch = 0] |> dplyr::select(-"var_id2")
    need_to_flip_id <- gwas_dataset2_fi[unique_fnf_norm_qtl_varID, on = .(var_id2 = var_id), nomatch = 0] |> 
      dplyr::select(-"var_id1") |> dplyr::rename("var_id1" ="var_id2")
    
    for(k in qtl_leadsig_test$value){
      
      pheid_fnf <- response_fnf_results |> dplyr::filter(var_id %in% k) |> 
        dplyr::select("phe_id","SYMBOL","rsID","rank","var_id","minor_allele","MAF","interaction_pval")
      
      qtl_pos <- as.numeric(str_split(k,":",simplify =TRUE)[,2])
      
      
      min_pos_gwas <- gwas_testPos - 250000
      max_pos_gwas <- gwas_testPos + 250000
      
      min_pos_qtl <- qtl_pos - 250000
      max_pos_qtl <- qtl_pos + 250000
      
      min_region <- min(min_pos_qtl,min_pos_gwas)
      max_region <- max(max_pos_qtl,max_pos_gwas)
      
      
      # For non-matching rows, use the regular format
      non_matching <- gwas_dataset2_fi[!var_id1 %in% unique_pbs_norm_qtl_varID$var_id &
                                         !var_id2 %in% unique_pbs_norm_qtl_varID$var_id] |> dplyr::select(-"var_id2")
      
      
      # Combine the results
      gwas_dataset2_flip <- rbindlist(list(matched_id, need_to_flip_id,non_matching), fill = TRUE)
      
      for(p in pheid_fnf$phe_id){
        
        pheid_fnf_filtered <- pheid_fnf |> dplyr::filter(phe_id == p)
        
        if(nrow(pheid_fnf_filtered) == 0) {
          print(paste("No FNF data for phenotype:", p))
          next  # Skip to the next iteration of the loop
        }
        
        #-------------------------------------------------------------------------
        # QTL
        
        # FNF analysis
        fnf_qtl_region <- fnf_norm_qtl_alleles |> dplyr::filter(phe_id %in% pheid_fnf_filtered$phe_id)
        fnf_qtl_region <- copy(fnf_qtl_region)
        fnf_qtl_region[, snpchr_loc := process_var_id(fnf_qtl_region$var_id)]
        
        fnf_qtl_region_fi <- inner_join(fnf_qtl_region, maf_id_add, by=c("var_id"="var_id"))
        
        fnf_qtl_region_subset <- fnf_qtl_region_fi %>%
          dplyr::filter(var_from >= min_region & var_from <= max_region) 
        
        # PBS analysis (using the same phenotype ID as FNF)
        pbs_qtl_region <- pbs_norm_qtl_chr |> dplyr::filter(phe_id == p)
        pbs_qtl_region <- copy(pbs_qtl_region)
        pbs_qtl_region[, snpchr_loc := process_var_id(pbs_qtl_region$var_id)]
        
        pbs_qtl_region_fi <- inner_join(pbs_qtl_region, maf_id_add, by=c("var_id"="var_id"))
        
        pbs_qtl_region_subset <- pbs_qtl_region_fi %>%
          dplyr::filter(var_from >= min_region & var_from <= max_region)
        
        
        
        # Gwas analysis
        
        gwas_dataset2_flip_region <- gwas_dataset2_flip %>%
          dplyr::filter(hg38pos >= min_region &
                          hg38pos <= max_region) 
        
        # Remove duplicates based on specific columns, e.g., 'var_id1'
        gwas_dataset2_final_range <- unique(gwas_dataset2_flip_region, by = "var_id1")

        
        gwas_dataset_list <- list(beta = gwas_dataset2_final_range$beta,
                                  snp = gwas_dataset2_final_range$var_id1,
                                  MAF = gwas_dataset2_final_range$MAF,
                                  pvalues= gwas_dataset2_final_range$pvalues,
                                  type="cc",
                                  s = fraction_cases$fraction_cases,
                                  N = fraction_cases$total_samples)
        
        # Prepare dataset lists
        fnf_dataset_list <- list(pvalues = fnf_qtl_region_subset$nom_pval,
                                 beta = fnf_qtl_region_subset$slope,
                                 N = n,
                                 MAF= fnf_qtl_region_subset$MAF,
                                 snp = fnf_qtl_region_subset$var_id,
                                 type = "quant")
        
        pbs_dataset_list <- list(pvalues = pbs_qtl_region_subset$nom_pval,
                                 beta = pbs_qtl_region_subset$slope,
                                 N = n,
                                 MAF= pbs_qtl_region_subset$MAF,
                                 snp = pbs_qtl_region_subset$var_id,
                                 type = "quant")
        
        # coloc run
        coloc_result_name <- paste0(subtype,"_",gwas_testID,"_",k,"_",p)
        print(coloc_result_name)
        
        coloc_result_fnf <- coloc.abf(dataset1 = fnf_dataset_list, dataset2 = gwas_dataset_list)
        coloc_result_pbs <- coloc.abf(dataset1 = pbs_dataset_list, dataset2 = gwas_dataset_list)
        
        coloc_summary_fnf <- c(coloc_result_fnf$summary[["PP.H3.abf"]], coloc_result_fnf$summary[["PP.H4.abf"]])
        coloc_summary_pbs <- c(coloc_result_pbs$summary[["PP.H3.abf"]], coloc_result_pbs$summary[["PP.H4.abf"]])
        
        # Debugging information
        cat("Debugging information:\n")
        cat("subtype:", subtype, "\n")
        cat("phe_id:", pheid_fnf_filtered$phe_id, "\n")
        cat("gene:", pheid_fnf_filtered$SYMBOL, "\n")
        cat("GWAS_rsID:", gwas_leadsig_test$rsID, "\n")
        cat("GWAS_pos:", gwas_leadsig_test$`CHR:hg38POS`, "\n")
        cat("GWAS_ref:", gwas_leadsig_test$ldbuddy_ref, "\n")
        cat("GWAS_alt:", gwas_leadsig_test$ldbuddy_alt, "\n")
        cat("sQTL_rsID:", pheid_fnf_filtered$rsID, "\n")
        cat("sQTL_rank:", pheid_fnf_filtered$rank, "\n")
        cat("sQTL_var_id:", pheid_fnf_filtered$var_id, "\n")
        cat("sQTL_minor_allele:", pheid_fnf_filtered$minor_allele, "\n")
        cat("sQTL_MAF:", pheid_fnf_filtered$MAF, "\n")
        cat("sQTL_interaction_pval:", pheid_fnf_filtered$interaction_pval, "\n")
        cat("PBS_coloc_H3:", coloc_summary_pbs[1], "\n")
        cat("PBS_coloc_H4:", coloc_summary_pbs[2], "\n")
        cat("FNF_coloc_H3:", coloc_summary_fnf[1], "\n")
        cat("FNF_coloc_H4:", coloc_summary_fnf[2], "\n")

        coloc_results_table <- data.frame(
          subtype = subtype,
          phe_id = pheid_fnf_filtered$phe_id,
          gene = pheid_fnf_filtered$SYMBOL,
          GWAS_rsID = gwas_leadsig_test$rsID,
          GWAS_pos = gwas_leadsig_test$`CHR:hg38POS`,
          GWAS_ref = gwas_leadsig_test$ldbuddy_ref,
          GWAS_alt = gwas_leadsig_test$ldbuddy_alt,
          sQTL_rsID = pheid_fnf_filtered$rsID,
          sQTL_rank = pheid_fnf_filtered$rank,
          sQTL_var_id = pheid_fnf_filtered$var_id,
          sQTL_minor_allele = pheid_fnf_filtered$minor_allele,
          sQTL_MAF = pheid_fnf_filtered$MAF,
          sQTL_interaction_pval = pheid_fnf_filtered$interaction_pval,
          PBS_coloc_H3 = coloc_summary_pbs[1],
          PBS_coloc_H4 = coloc_summary_pbs[2],
          FNF_coloc_H3 = coloc_summary_fnf[1],
          FNF_coloc_H4 = coloc_summary_fnf[2]
        )
        
        coloc_results.df <- rbind(coloc_results.df, coloc_results_table)
        
        # Store the result in the list
        coloc_results[[coloc_result_name]] <- coloc_result
      }
    }
  }
  
  save(coloc_results, file= paste0("output/coloc/coloc_results_fnf_prob_calc_",subtype,".rda" ))
}

colnames(coloc_results.df) <- c("subtype","phe_id","gene","GWAS_rsID","GWAS_pos","GWAS_ref","GWAS_alt",
                                "sQTL_rsID","sQTL_rank","sQTL_var_id","sQTL_minor_allele","sQTL_MAF","sQTL_interaction_pval",
                                "PBS_coloc_H3","PBS_coloc_H4", "FNF_coloc_H3","FNF_coloc_H4")
write.table(coloc_results.df, file="output/coloc/coloc_result_table_fnf_prob_calc_.txt",sep="\t",quote=F,col.names = TRUE, row.names = FALSE)
    