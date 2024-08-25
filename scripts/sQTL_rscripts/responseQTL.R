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
intron_test <- c(pbsQTL_bindConditional$phe_id, fnfQTL_bindConditional$phe_id) |> unique()
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


#------------------------------------------------------------------------------
#qtl data
header_cond <- c("phe_id","phe_chr","phe_from",
                 "phe_to","phe_strd","n_var_in_cis","dist_phe_var","var_id",
                 "var_chr","var_from","var_to","rank",
                 "fwd_pval","fwd_r_squared","fwd_slope","fwd_best_hit","fwd_sig",
                 "bwd_pval","bwd_r_squared","bwd_slope","bwd_best_hit","bwd_sig")
cond_pbs <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/conditional_pbs/conditional_psb_top_variants.txt")
colnames(cond_pbs) <- header_cond
cond_fnf <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/conditional_fnf/conditional_fnf_top_variants.txt")
colnames(cond_fnf) <- header_cond

#QTL
header <- read.table("scripts/sQTL_rscripts/header.txt")
new_header <- cbind(header, "qval","therehod")
sigQTL_fnf <- read.table("output/01.qtltools_re/01.significant/fnf_0.05_pc4.significant.txt") |> as.data.frame()
sigQTL_pbs <- read.table("output/01.qtltools_re/01.significant/pbs_0.05_pc5.significant.txt") |> as.data.frame()
colnames(sigQTL_pbs) <- new_header
colnames(sigQTL_fnf) <- new_header

# threadhold and qval subset



# Lets combine all the qtls and conditional analysis before running the response QTL. 
#First I am going to check that what phed_id has the other ranks
cond_pbs_ID <- cond_pbs %>%
  dplyr::mutate(id = paste(phe_id, var_id, sep = "_")) 
cond_fnf_ID <- cond_fnf %>%
  dplyr::mutate(id = paste(phe_id, var_id, sep = "_")) 

sigQTL_pbs_ID <- sigQTL_pbs %>%
  dplyr::mutate(id = paste(phe_id, var_id, sep = "_")) 
sigQTL_fnf_ID <- sigQTL_fnf %>%
  dplyr::mutate(id = paste(phe_id, var_id, sep = "_")) 


# Filter rows where both phe_id and var_id don't match sigQTL_pbs
pbs_conditionalOnly <- cond_pbs_ID %>%
  dplyr::filter(!(id %in% sigQTL_pbs_ID$id)) |> 
  dplyr::filter(!(rank %in% "0"))

fnf_conditionalOnly <- cond_fnf_ID %>%
  dplyr::filter(!(id %in% sigQTL_fnf_ID$id)) |> 
  dplyr::filter(!(rank %in% "0"))

pbs_conditional_filter <- cond_pbs_ID %>%
  dplyr::filter(!(id %in% sigQTL_pbs_ID$id)) |> 
  dplyr::filter(rank %in% "0") |> dplyr::pull(var_id)

fnf_conditional_filter <- cond_fnf_ID %>%
  dplyr::filter(!(id %in% sigQTL_fnf_ID$id)) |> 
  dplyr::filter(rank %in% "0")|> dplyr::pull(var_id)


conditional_filterout_snps <- unique(c(pbs_conditional_filter,fnf_conditional_filter))
save(conditional_filterout_snps, file="output/01.qtltools_re/conditional_filterout_snps.txt")

# pbs_sqtl_notmathced_conditional_rank0 <- sigQTL_pbs_ID %>%
#   dplyr::filter(!(id %in% cond_pbs_ID$id))
# 
# fnf_sqtl_notmathced_conditional_rank0 <- sigQTL_fnf_ID %>%
#   dplyr::filter(!(id %in% cond_fnf_ID$id))

# Now, let's adjust the ranks
#pbs_conditional_adjustedRank <- adjust_ranks_for_zero(pbs_conditionalOnly)
#fnf_conditional_adjustedRank <- adjust_ranks_for_zero(fnf_conditionalOnly)

#Subset the columns

qtl_pbs_thres_qval <- sigQTL_pbs_ID |> dplyr::select(c("phe_id", "qval","therehod")) 
qtl_fnf_thres_qval <- sigQTL_fnf_ID |> dplyr::select(c("phe_id", "qval","therehod"))

sigQTL_pbs_subset <- sigQTL_pbs_ID |> dplyr::select(c( "phe_id", "phe_chr","phe_from","phe_to","phe_strd",
                                                       "n_var_in_cis","dist_phe_var","var_id","var_chr", "var_from",
                                   "qval","therehod","id")) |>
  dplyr::mutate(rank = 0)

sigQTL_fnf_subset <- sigQTL_fnf_ID |> dplyr::select(c( "phe_id", "phe_chr","phe_from","phe_to","phe_strd",
                                                        "n_var_in_cis","dist_phe_var","var_id","var_chr", "var_from",
                                                        "qval","therehod","id")) |>
  dplyr::mutate(rank = 0)

pbs_cond_temp <- pbs_conditionalOnly |> dplyr::select(c( "phe_id", "phe_chr","phe_from","phe_to","phe_strd",
                                         "n_var_in_cis","dist_phe_var","var_id","var_chr", "var_from" ,"rank",
                                         "id")) |>
  left_join(qtl_pbs_thres_qval, by = "phe_id")
fnf_cond_temp <- fnf_conditionalOnly |> dplyr::select(c( "phe_id", "phe_chr","phe_from","phe_to","phe_strd",
                                                         "n_var_in_cis","dist_phe_var","var_id","var_chr", "var_from" ,"rank",
                                                         "id")) |>
  left_join(qtl_fnf_thres_qval, by = "phe_id")

# Merge cond and qtl

pbsQTL_bindConditional <- bind_rows(sigQTL_pbs_subset,pbs_cond_temp) 
fnfQTL_bindConditional <- bind_rows(sigQTL_fnf_subset,fnf_cond_temp)



# sample_meta_data

intron_test <- c(pbsQTL_bindConditional$phe_id, fnfQTL_bindConditional$phe_id) |> unique()
psi_matrix <- ratios_fnf[rownames(ratios_fnf) %in% intron_test,] |> as.data.frame()
intronID <- rownames(psi_matrix)
psi_ratio <- apply(psi_matrix, 2, as.numeric)
rownames(psi_ratio) <- intronID
psi_ratio.df <- t(psi_ratio) |> as.data.frame()

pbsQTL_annotated_ensg <- annotate_genes(pbsQTL_bindConditional, hg38_intron_sub_select_first, leafcutter_pheno_subset, txdb_genes)
fnfQTL_annotated_ensg <- annotate_genes(fnfQTL_bindConditional, hg38_intron_sub_select_first, leafcutter_pheno_subset, txdb_genes)

#adding beta and MAF and anova p-value
#read MAF  and beta
#read beta both pbs and fnf
 
header_nominal <-  read.table("scripts/sQTL_rscripts/nominal_header.txt",sep = " ") 
                     
beta_pbs <- fread("output/01.qtltools_re/nominal_pbs/pc5_allchr.pbs.cis")
colnames(beta_pbs) <- paste(header_nominal)

beta_pbs_subset <- beta_pbs %>% dplyr::select(c("nom_pval","phe_id","var_id","slope","slope_se")) %>% 
  dplyr::rename("PBS_p" = "nom_pval","PBS_beta" = slope,"PBS_beta_se"=slope_se)
  
beta_fnf <- fread("output/01.qtltools_re/nominal_fnf/pc4_allchr.fnf.cis")
colnames(beta_fnf) <- paste(header_nominal)

beta_fnf_subset <-beta_fnf %>% dplyr::select(c("nom_pval","phe_id","var_id","slope","slope_se")) %>% 
  dplyr::rename("FNF_p" = "nom_pval", "FNF_beta"=slope,"FNF_beta_se" =slope_se)

#saveRDS(beta_pbs_subset,file="./01.qtltools_re/beta_pbs_pc5.rda")
#saveRDS(beta_fnf_subset,file= "./01.qtltools_re/beta_fnf_pc4.rda")


beta_pbs_subset <- readRDS("output/01.qtltools_re/beta_pbs_pc5.rda")
beta_fnf_subset <- readRDS("output/01.qtltools_re/beta_fnf_pc4.rda")

#MAF
maf_all <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/freq/cqtl.frq")
maf_subset <-  maf_all %>% dplyr::select(c("SNP","MAF","A1")) %>% 
  dplyr::rename("var_id" =SNP, "minor_allele" =A1) 

# Function to read and process BED files
read_bed_file <- function(file_path) {
  fread(file_path, col.names = c("chr", "start", "var_id", "ref", "alt", "rsID"))
}

# Get list of BED files
bed_files <- list.files("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/rsid",pattern = "_rsid.bed",full.names = TRUE)
all_bed_data <- rbindlist(lapply(bed_files, read_bed_file))
all_bed_dataID <- all_bed_data |> dplyr::select(c("var_id","rsID"))

#merge rsIDs
maf_id_add <- left_join(maf_subset , all_bed_dataID,by="var_id") 

#save(maf_id_add, file="/work/users/s/e/seyoun/CQTL_sQTL/output/geno/freq/maf_id_add.rds")
load("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/freq/maf_id_add.rds")




#PBS_adding beta and maf
pbs_QTL_cond_beta_maf <- pbsQTL_annotated_ensg %>% 
  left_join(maf_id_add, by="var_id") %>%
  left_join(beta_pbs_subset,by=c("phe_id" = "phe_id", "var_id" = "var_id")) %>%
  left_join(beta_fnf_subset,by=c("phe_id" = "phe_id", "var_id" = "var_id")) %>%
  mutate(delta_beta = FNF_beta - PBS_beta)|>
  dplyr::mutate(ensg = gsub("\\..*$", "", ensg))


fnf_QTL_cond_beta_maf <- fnfQTL_annotated_ensg %>% 
  left_join(maf_id_add, by="var_id") %>%
  left_join(beta_pbs_subset,by=c("phe_id" = "phe_id", "var_id" = "var_id")) %>%
  left_join(beta_fnf_subset,by=c("phe_id" = "phe_id", "var_id" = "var_id")) %>%
  mutate(delta_beta = FNF_beta - PBS_beta)|>
  dplyr::mutate(ensg = gsub("\\..*$", "", ensg))


ensg_id <- unique(c(pbs_QTL_cond_beta_maf$ensg, fnf_QTL_cond_beta_maf$ensg))

annots_unique <- AnnotationDbi::select(org.Hs.eg.db, ensg_id,
                                columns=c("SYMBOL"), keytype="ENSEMBL") |>
  distinct(ENSEMBL, .keep_all = TRUE) |>
  dplyr::rename("ensg" = "ENSEMBL")

pbs_QTL_cond_final <- left_join(pbs_QTL_cond_beta_maf, annots_unique, by = "ensg")
fnf_QTL_cond_final <- left_join(fnf_QTL_cond_beta_maf, annots_unique, by = "ensg")


saveRDS(pbs_QTL_cond_final, file="output/01.qtltools_re/pbs_cond_sig_beta_maf_added.rds")
saveRDS(fnf_QTL_cond_final, file="output/01.qtltools_re/fnf_cond_sig_beta_maf_added.rds")
pbs_QTL_cond_final <- readRDS("output/01.qtltools_re/pbs_cond_sig_beta_maf_added.rds")
fnf_QTL_cond_final <- readRDS("output/01.qtltools_re/fnf_cond_sig_beta_maf_added.rds")



degenes_fnf_shrink <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/quant/de_genes_results.csv") |>
  dplyr::mutate(gene_id = gsub("\\..*$", "", gene_id))
annots <- AnnotationDbi::select(org.Hs.eg.db, degenes_fnf_shrink$gene_id,
                                columns=c("SYMBOL"), keytype="ENSEMBL")
annots <- dplyr::rename(annots, gene_id = ENSEMBL)
annots_unique <- annots %>% distinct(gene_id, .keep_all = TRUE)
degenes_fnf_shrink_symbol <- left_join(degenes_fnf_shrink, annots_unique, by = "gene_id")
degenes_fnf_shrink_symbol <- degenes_fnf_shrink_symbol %>% dplyr::select("gene_id","SYMBOL","log2FoldChange","padj")


  
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
    
    summarize_count <- meta_data_fi |> 
      group_by(genotype) |> 
      summarise(n = n(), .groups = "drop")
    
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
    conditional_result$lrt_chisq[i] <- result$Chisq[2]  # Chi-square statistic
    conditional_result$lrt_df[i] <- result$Df[2]  # Degrees of freedom
    conditional_result$minor_alle_count[i] <-   summarize_count$n[3]/2
    # # Save results
    # results_list[[i]] <- list(
    #   result = result,
    #   qtl_model = reduced_model,
    #   interaction = full_model,
    #   pval = result[2, "Pr(>Chisq)"]
    #)
  }
  return(conditional_result)
}

response_pbs_results <- analyze_response_qtl(conditional_result=pbs_QTL_cond_final,
                                             psi=all_psi_inv, 
                                             geno.df=all_geno_transpose_df,
                                             meta_df=meta_catl_simple, 
                                             covariates_df=pc5_cov_all)



response_fnf_results <- analyze_response_qtl(conditional_result=fnf_QTL_cond_final, 
                                             psi=all_psi_inv, 
                                             geno.df=all_geno_transpose_df,
                                             meta_df=meta_catl_simple,
                                             covariates_df=pc4_cov_all)




#saveRDS(response_pbs_results, file="output/01.qtltools_re/conditional_pbs/response_pbs_results.rds")
#saveRDS(response_fnf_results, file="output/01.qtltools_re/conditional_pbs/response_pbs_results.rds")

saveRDS(response_pbs_results, file="output/01.qtltools_re/conditional_pbs/response_pbs_resInv.rds")
saveRDS(response_fnf_results, file="output/01.qtltools_re/conditional_pbs/response_fnf_resInv.rds")

response_pbs_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_pbs_resInv.rds")
response_fnf_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_fnf_resInv.rds")

##Filter the significant and adding beta value and MAF
threshold <- as.numeric(0.05)

sig_responsePBS <- response_pbs_results %>% 
  dplyr::filter(interaction_pval  < 0.05) %>%
  dplyr::filter(abs(delta_beta ) > 0.2) |>
  dplyr::filter(minor_alle_count >=5) |>
  arrange(interaction_pval)



sig_responseFNF <- response_fnf_results %>% 
  dplyr::filter(interaction_pval  < 0.05) %>%
  dplyr::filter(abs(delta_beta ) > 0.2) |>
  dplyr::filter(minor_alle_count >=5) |>
  arrange(interaction_pval)
 

fnf_pvals <- sapply(response_fnf_results, function(x) x$pval)
final_annotated_fnf$annova_pval <- fnf_pvals

response_conf_fnf_interection <- final_annotated_fnf %>% 
  dplyr::filter(annova_pval < 0.05) %>%
  dplyr::filter(abs(Beta_diff) >=0.2)
  pull(ensg) %>%
  unique() 



pbs_fnf_response_sGenes <- tibble(values = unique(c(response_pbs_ensg, response_fnf_ensg))) %>%
  mutate(PBS = values %in% response_pbs_ensg,
         FNF = values %in% response_fnf_ensg)

ggplot(pbs_fnf_response_sGenes, aes(A = PBS, B = FNF)) +
  geom_venn(set_names = c("PBS", "FN-f"), 
            fill_color = c("#CBD5E8","#FDCDAC"), 
            stroke_color = NA, auto_scale = TRUE, show_percentage = FALSE,
            text_size = 3, set_name_size = 0) +
  coord_fixed()  +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"))
ggsave(filename = "output/results_plots/sqtl_plots/venn_sGenes.pdf", 
       width = 4, height = 4, units = "in")
save(venn_diagram, file = "output/results_plots/sqtl_plots/venn_sGenes.rda")
venn_font(venn_diagram, font = "Helvetica")




