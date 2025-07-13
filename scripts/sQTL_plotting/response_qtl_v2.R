# Response sQTL
## Author: Seyoun Byun
## Date: 07.05.2024
## Edited: 07.08.2024
#-------------------------------------------------------------------------------
setwd("/work/users/s/e/seyoun/CQTL_sQTL")
library(tidyverse)
library(tibble)
library(janitor)
library(vcfR)
library(car)
library(lme4)
library(data.table)
library(org.Hs.eg.db)
library(qvalue)
source("./scripts/sQTL_rscripts/utils.R")


pbs_QTL_cond_final <- readRDS("output/01.qtltools_re/pbs_cond_sig_beta_maf_added.rds")
fnf_QTL_cond_final <- readRDS("output/01.qtltools_re/fnf_cond_sig_beta_maf_added.rds")

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
ctl_fnf_ratio <- fread("output/clu_fnf/ratio_fnf.txt") |> as.data.frame()
rownames(ctl_fnf_ratio) <- ctl_fnf_ratio$Junction
ratios_fnf <- ctl_fnf_ratio[,-1]
intron_test <- c(pbs_QTL_cond_final$phe_id, fnf_QTL_cond_final$phe_id) |> unique()
psi_matrix <- ratios_fnf[rownames(ratios_fnf) %in% intron_test,] |> as.data.frame()
intronID <- rownames(psi_matrix)
psi_ratio <- apply(psi_matrix, 2, as.numeric)
rownames(psi_ratio) <- intronID
psi_ratio.df <- t(psi_ratio) |> as.data.frame()

#covariates
pc5_cov_all <- fread("output/gtex_cluster/qtltools_prep/covariates_PC5") %>%
  t() %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  `colnames<-`(.[1, ]) %>%
  .[-1, ] %>%
  cbind(id = rownames(.), .) %>%
  mutate(across(-id, as.numeric))

pc4_cov_all <- fread("output/gtex_cluster/qtltools_prep/covariates_PC4") %>%
  t() %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  `colnames<-`(.[1, ]) %>%
  .[-1, ] %>%
  cbind(id = rownames(.), .) %>%
  mutate(across(-id, as.numeric))


#meta #exclude samples for  ['AM7352', 'AM7244']
load("output/combined_meta_data.RData")
selected_columns <- c("Donor","ID","Condition","Sex", "Age","FragmentBatch","RIN","RNAextractionKitBatch","RNAshippedDate") #select column needed it based
meta_data <- combined_data %>% dplyr::select(all_of(selected_columns))
#Ancestry 
ancestry_df <- fread("/proj/phanstiel_lab/Data/processed/CQTL/geno/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_COA6_COA7/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_COA6_COA7_predictedAncestry.csv")
ancestry_df$Donor <- sub("^.*_(AM[0-9]+)_.*$", "\\1", ancestry_df$Donor)

ancestry_OA_df <- fread("/proj/phanstiel_lab/Data/processed/CQTL/geno/COA8_OA/ancestry/CQTL_COA8_predictedAncestry.csv")
ancestry_OA_df_filtered <- ancestry_OA_df %>%
  dplyr::filter(grepl("^OA\\d+", Donor))
ancestry_OA_df_filtered$Donor <- gsub("_r2","",ancestry_OA_df_filtered$Donor)

ancestry_cqtl <-rbind(ancestry_df,ancestry_OA_df_filtered)
meta_cqtl <- merge(meta_data,ancestry_cqtl,by="Donor",all.x=TRUE)
#write.table(meta_cqtl, file = "output/clu_fnf/meta_cqtl",sep='\t',quote=F,row.names=F,col.names=T)
meta_cqtl <- fread("output/clu_fnf/meta_cqtl")
meta_catl_simple <- meta_cqtl %>%
  dplyr::filter(Condition %in% c("CTL","FNF")) %>%
  dplyr::select(-Sex,-Age,-FragmentBatch,-RIN,-RNAextractionKitBatch,-RNAshippedDate,-Predicted_Ancestry) %>%
  dplyr::filter(!str_detect(Donor, 'AM7352|AM7244'))






#interaction testing-------------------------------------------------

analyze_response_qtl <- function(conditional_result, psi, geno.df, meta_df, covariates_df) {
  
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
    reduced_model <- lme4::lmer(reduced_equation, meta_combined_all, REML=FALSE)
    full_model <- lme4::lmer(full_equation, meta_combined_all, REML=FALSE)
    
    result <- anova(reduced_model, full_model)
    conditional_result$interaction_pval[i] <-  result$`Pr(>Chisq)`[2]
    conditional_result$lrt_chisq <- result$Chisq[2]  # Chi-square statistic
    conditional_result$lrt_df <- result$Df[2]  # Degrees of freedom
  }
  return(conditional_result)
}

response_pbs_results <- analyze_response_qtl(conditional_result=pbs_QTL_cond_final,
                                             psi=psi_ratio, 
                                             geno.df=all_geno_transpose_df,
                                             meta_df=meta_catl_simple, 
                                             covariates_df=pc5_cov_all)
response_fnf_results <- analyze_response_qtl(conditional_result=fnf_QTL_cond_final, 
                                             psi=psi_ratio, 
                                             geno.df=all_geno_transpose_df,
                                             meta_df=meta_catl_simple,
                                             covariates_df=pc4_cov_all)





saveRDS(response_pbs_results, file="output/01.qtltools_re/conditional_pbs/response_pbs_results.rds")
saveRDS(response_fnf_results, file="output/01.qtltools_re/conditional_pbs/response_fnf_results.rds")





