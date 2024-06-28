library(tidyverse)
library(tibble)
library(janitor)
library(vcfR)
library(car)
library(lme4)
library(data.table)
source("./utils.R")
setwd("/work/users/s/e/seyoun/CQTL_sQTL/output")



#-------------------------------------------------------------------------------
#geno_data 
pbsGeno_raw <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/pbs_geno/06.subset_sigSNps/recodeA_pbs.raw")
colnames(pbsGeno_raw) <- sub("_.*", "", colnames(pbsGeno_raw))

fnfGeno_raw <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/fnf_geno/06.subset_sigSNps/recodeA_fnf.raw")
colnames(fnfGeno_raw) <- sub("_.*", "", colnames(fnfGeno_raw))
all_geno_matrix <- bind_rows(pbsGeno_raw,fnfGeno_raw)
all_geno_transpose_df <- all_geno_matrix %>%
  column_to_rownames(var = "IID")  %>%
  dplyr::select(-FID, -MAT,-PAT, -SEX, -PHENOTYPE) %>% 
  as.data.frame() %>%  # Ensure it's still a data frame
  t()  

#Normalized PSI value 
ctl_fnf_ratio <- fread("clu_fnf/ratio_fnf.txt") |> as.data.frame()
rownames(ctl_fnf_ratio) <- ctl_fnf_ratio$Junction
ratios_fnf <- ctl_fnf_ratio[,-1]

#covariates
pc5_cov_all <- fread("gtex_cluster/qtltools_prep/covariates_PC5") %>%
  t() %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  `colnames<-`(.[1, ]) %>%
  .[-1, ] %>%
  cbind(id = rownames(.), .) %>%
  mutate(across(-id, as.numeric))

pc4_cov_all <- fread("gtex_cluster/qtltools_prep/covariates_PC4") %>%
  t() %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  `colnames<-`(.[1, ]) %>%
  .[-1, ] %>%
  cbind(id = rownames(.), .) %>%
  mutate(across(-id, as.numeric))

#meta---------------------------------------------------------------------------
meta_cqtl <- fread("clu_fnf/meta_cqtl")
meta_catl_simple <- meta_cqtl %>%
  dplyr::filter(Condition %in% c("CTL","FNF")) %>%
  dplyr::select(-Sex,-Age,-FragmentBatch,-RIN,-RNAextractionKitBatch,-RNAshippedDate,-Predicted_Ancestry) %>%
  dplyr::filter(!str_detect(Donor, 'AM7352|AM7244'))


header_cond <- c("phe_id","phe_chr","phe_from",
                 "phe_to","phe_strd","n_var_in_cis","dist_phe_var","var_id",
                 "var_chr","var_from","var_to","rank",
                 "fwd_pval","fwd_r_squared","fwd_slope","fwd_best_hit","fwd_sig",
                 "bwd_pval","bwd_r_squared","bwd_slope","bwd_best_hit","bwd_sig")
cond_pbs <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/conditional_pbs/conditional_psb_top_variants.txt")
colnames(cond_pbs) <- header_cond
cond_fnf <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/conditional_fnf/conditional_fnf_top_variants.txt")
colnames(cond_fnf) <- header_cond


# sample_meta_data

intron_test <- c(cond_pbs$phe_id, cond_fnf$phe_id)
psi_matrix <- ratios_fnf[rownames(ratios_fnf) %in% intron_test,] |> as.data.frame()
intronID <- rownames(psi_matrix)
psi_ratio <- apply(psi_matrix, 2, as.numeric)
rownames(psi_ratio) <- intronID

#interaction testing-------------------------------------------------

analyze_response_qtl <- function(conditional_result, psi, geno.df, meta_df, covariates_df) {
  results_list <- list()
  
  for (i in 1:nrow(conditional_result)) {
    print(i)
    intronID <- conditional_result$phe_id[i]
    
    psi_matrix <- psi[rownames(psi) %in% intronID, ]
    psi_data <- data.frame(sampleID = names(psi_matrix),
                           psi = as.double(psi_matrix))
    
    variantID <- conditional_result$var_id[i]
    geno_data <- geno.df[rownames(geno.df) %in% variantID, ] %>%
      as.data.frame() %>%
      rownames_to_column(var = "sampleID")
    colnames(geno_data)[2] <- "genotype"
    
    meta_data_fi <- geno_data %>%
      left_join(meta_df,by = c("sampleID"="ID")) %>%
      left_join(psi_data,by = c("sampleID"="sampleID")) %>%
      left_join(covariates_df, by=c("sampleID"= "id"))
    
    meta_combined_all <- meta_data_fi %>%
       mutate(
        Donor = as.factor(Donor),
        Condition = ifelse(Condition == "CTL", 0, ifelse(Condition == "FNF", 1, NA))
       )
    
    # Define Formulas
    pc_columns <- paste(colnames(meta_combined_all)[6:length(colnames(meta_combined_all))], collapse=" + ")
    reduced_equation <- as.formula(paste("psi ~ genotype + Condition + (1|Donor)", pc_columns, sep = "+"))
    full_equation <- as.formula(paste("psi ~ genotype + Condition + Condition:genotype + (1|Donor)", pc_columns, sep = "+"))

    # Compare two models
    reduced_model <- lme4::lmer(reduced_equation, meta_combined_all)
    full_model <- lme4::lmer(full_equation, meta_combined_all)

    result <- anova(reduced_model, full_model)
     result$`Pr(>Chisq)`[2]
    # Save results
    results_list[[i]] <- list(
      result = result,
      qtl_model = reduced_model,
      interaction = full_model,
      pval = result[2, "Pr(>Chisq)"]
    )
  }
  return(results_list)
}


final_pbs_sig_qtl_cond <- readRDS("./01.qtltools_re/pbs_cond_sig_beta_maf_added.rds")
final_fnf_sig_qtl_cond <- readRDS("./01.qtltools_re/fnf_cond_sig_beta_maf_added.rds")

response_pbs_results <- analyze_response_qtl(final_pbs_sig_qtl_cond, psi_ratio, all_geno_transpose_df, meta_catl_simple, pc5_cov_all)
response_fnf_results <- analyze_response_qtl(final_fnf_sig_qtl_cond, psi_ratio, all_geno_transpose_df, meta_catl_simple, pc4_cov_all)

save(response_pbs_results, file="01.qtltools_re/conditional_pbs/response_pbs_results.Rdata")
save(response_fnf_results, file="01.qtltools_re/conditional_fnf/response_fnf_results.Rdata")




