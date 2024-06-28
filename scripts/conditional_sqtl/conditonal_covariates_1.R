# Making the covariates to run all the condtional so it can have p-value
## Author: Seyoun Byun
## Date: 04.18.2024
## Edited:
#-------------------------------------------------------------------------------
setwd("/work/users/s/e/seyoun/CQTL_sQTL")
library(data.table)
library(dplyr)
library(yaml)
library(readr)

#Making multiple different covariates files for both pbs and fnf 

pc5_pbs <- read.table("output/gtex_cluster/qtltools_prep/covariates_PC5") |> t() |> 
  as.data.frame(as.factor=FALSE)
colnames(pc5_pbs) <- pc5_pbs[1,]
pc5_pbs_cov <- pc5_pbs[-1,]

pc4_fnf <- read.table("output/gtex_cluster/qtltools_prep/covariates_PC4")|> t() |> 
  as.data.frame(as.factor=FALSE)
colnames(pc4_fnf) <- pc4_fnf[1,]
pc4_fnf_cov <- pc4_fnf[-1,]

#-------------------------------------------------------------------------------
#pbs- 

header_cond <- c("phe_id","phe_chr","phe_from",
                 "phe_to","phe_strd","n_var_in_cis","dist_phe_var","var_id",
                 "var_chr","var_from","var_to","rank",
                 "fwd_pval","fwd_r_squared","fwd_slope","fwd_best_hit","fwd_sig",
                 "bwd_pval","bwd_r_squared","bwd_slope","bwd_best_hit","bwd_sig")
cond_pbs <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/conditional_pbs/conditional_psb_top_variants.txt")
colnames(cond_pbs) <- header_cond


# First, check if each phe_id includes ranks other than 0
phe_id_with_multiple_ranks <- cond_pbs %>%
  group_by(phe_id) %>%
  summarize(has_multiple_ranks = any(rank != 0)) %>%
  dplyr::filter(has_multiple_ranks)

# Filter original data to include only those IDs
cond_pbs_filtered <- cond_pbs %>%
  dplyr::filter(phe_id %in% phe_id_with_multiple_ranks$phe_id ) %>%
  dplyr::select(phe_id) %>%
  unique()




####Clenup snps
pbsGeno_raw <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/pbs_geno/06.subset_sigSNps/recodeA_pbs.raw")
colnames(pbsGeno_raw) <- sub("_.*", "", colnames(pbsGeno_raw))

# for (i in nrow(cond_pbs_filtered) ) {
#   print(i)
#   temp_table <-cond_pbs %>%
#     filter(phe_id %in% cond_pbs_filtered[i])
#   bed_dir <- paste0("output/gtex_cluster/qtltools_prep/ctrvsfnf_qqnorm_",  temp_table$phe_chr[1],".bed.gz")
#   
#   if(nrow(temp_table) == 2){
#     rank0_var_ids <- temp_table %>%
#       filter(rank == 0) %>%
#       pull(var_id)
#     phe_id <- temp_table$phe_id[1]
#     bed_txt <- paste0("output/01.qtltools_re/conditional_covariates/pbs_cov","/",phe_id,".txt")
#     write.table(phe_id,file=bed_txt,
#                 sep='\t',quote=F,row.names=F,col.names=F)
#     pbsGeno_subset <- pbsGeno_raw %>%
#       select(c("IID", rank0_var_ids[rank0_var_ids %in% names(pbsGeno_raw)]))
#        
#     pbsGeno_subset <- pbsGeno_subset %>%
#       rename(id=IID)
#     
#     pbsGeno_joined <- pbsGeno_subset %>%
#       left_join(pc5_pbs_cov, by = "id")
#     pbs_geno_joined_cov <- t(pbsGeno_joined)
#     colnames(pbs_geno_joined_cov) <- pbs_geno_joined_cov[1,]
#     pbs_geno_joined_cov <- pbs_geno_joined_cov[-1,]
#     pbs_geno_joined_cov <- cbind(rownames(pbs_geno_joined_cov),pbs_geno_joined_cov)
#     colnames(pbs_geno_joined_cov)[1] <- c("id")
#     cov_dir <- paste0("output/01.qtltools_re/conditional_covariates/pbs_cov","/","covariates_PC5_",phe_id,"_rank1")
#     write.table(pbs_geno_joined_cov,
#                 file=cov_dir ,
#                 sep='\t',quote=F,row.names=F,col.names=T)
#     norm_dir <- paste0("output/01.qtltools_re/cond_norm_pbs/",phe_id,"_rank1")
#     perm_dir <- paste0("output/01.qtltools_re/cond_perm_pbs/",phe_id,"_rank1")
#     system(paste0("scripts/conditional_sqtl/qtltools_rerun.sh ",
#                   "output/geno/pbs_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz ",
#                   bed_dir, " ",
#                   norm_dir, " ", 
#                   perm_dir, " ",
#                   bed_txt, " ",
#                   cov_dir, " ",
#                   phe_id))
#   }
#   
#   if(nrow(temp_table) == 3){
#     #rank1
#     rank0_1_var_ids <- temp_table %>%
#       filter(rank != 1) %>%
#       pull(var_id)
#     phe_id <- temp_table$phe_id[1]
#     bed_txt <- paste0("output/01.qtltools_re/conditional_covariates/pbs_cov","/",phe_id,".txt")
#     write.table(phe_id,file=bed_txt,
#                 sep='\t',quote=F,row.names=F,col.names=F)
#     pbsGeno_subset <- pbsGeno_raw %>%
#       select(c("IID", rank0_1_var_ids[rank0_1_var_ids %in% names(pbsGeno_raw)]))
#     
#     pbsGeno_subset <- pbsGeno_subset %>%
#       rename(id=IID)
#     
#     pbsGeno_joined <- pbsGeno_subset %>%
#       left_join(pc5_pbs_cov, by = "id")
#     pbs_geno_joined_cov <- t(pbsGeno_joined)
#     colnames(pbs_geno_joined_cov) <- pbs_geno_joined_cov[1,]
#     pbs_geno_joined_cov <- pbs_geno_joined_cov[-1,]
#     pbs_geno_joined_cov <- cbind(rownames(pbs_geno_joined_cov),pbs_geno_joined_cov)
#     colnames(pbs_geno_joined_cov)[1] <- c("id")
#     cov_dir <- paste0("output/01.qtltools_re/conditional_covariates/pbs_cov","/","covariates_PC5_",phe_id,"_rank1")
#     write.table(pbs_geno_joined_cov,
#                 file=cov_dir ,
#                 sep='\t',quote=F,row.names=F,col.names=T)
#     norm_dir <- paste0("output/01.qtltools_re/cond_norm_pbs/",phe_id,"_rank1")
#     perm_dir <- paste0("output/01.qtltools_re/cond_perm_pbs/",phe_id,"_rank1")
#     system(paste0("scripts/conditional_sqtl/qtltools_rerun.sh ",
#                   "output/geno/pbs_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz ",
#                   bed_dir, " ",
#                   norm_dir, " ", 
#                   perm_dir, " ",
#                   bed_txt, " ",
#                   cov_dir, " ",
#                   phe_id))
#     #rank2
#     rank0_2_var_ids <- temp_table %>%
#       filter(rank != 2) %>%
#       pull(var_id)
#     phe_id <- temp_table$phe_id[1]
#     pbsGeno_subset <- pbsGeno_raw %>%
#       select(c("IID", rank0_2_var_ids[rank0_2_var_ids %in% names(pbsGeno_raw)]))
#     
#     pbsGeno_subset <- pbsGeno_subset %>%
#       rename(id=IID)
#     
#     pbsGeno_joined <- pbsGeno_subset %>%
#       left_join(pc5_pbs_cov, by = "id")
#     pbs_geno_joined_cov <- t(pbsGeno_joined)
#     colnames(pbs_geno_joined_cov) <- pbs_geno_joined_cov[1,]
#     pbs_geno_joined_cov <- pbs_geno_joined_cov[-1,]
#     pbs_geno_joined_cov <- cbind(rownames(pbs_geno_joined_cov),pbs_geno_joined_cov)
#     colnames(pbs_geno_joined_cov)[1] <- c("id")
#     cov_dir <- paste0("output/01.qtltools_re/conditional_covariates/pbs_cov","/","covariates_PC5_",phe_id,"_rank2")
#     write.table(pbs_geno_joined_cov,
#                 file=cov_dir ,
#                 sep='\t',quote=F,row.names=F,col.names=T)
#     norm_dir <- paste0("output/01.qtltools_re/cond_norm_pbs/",phe_id,"_rank2")
#     perm_dir <- paste0("output/01.qtltools_re/cond_perm_pbs/",phe_id,"_rank2")
#     system(paste0("scripts/conditional_sqtl/qtltools_rerun.sh ",
#                   "output/geno/pbs_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz ",
#                   bed_dir, " ",
#                   norm_dir, " ", 
#                   perm_dir, " ",
#                   bed_txt, " ",
#                   cov_dir, " ",
#                   phe_id))
#   }
#   else(nrow(temp_table) == 4){
#     #rank1
#     rank0_1_var_ids <- temp_table %>%
#       filter(rank != 1) %>%
#       pull(var_id)
#     phe_id <- temp_table$phe_id[1]
#     bed_txt <- paste0("output/01.qtltools_re/conditional_covariates/pbs_cov","/",phe_id,".txt")
#     write.table(phe_id,file=bed_txt,
#                 sep='\t',quote=F,row.names=F,col.names=F)
#     pbsGeno_subset <- pbsGeno_raw %>%
#       select(c("IID", rank0_1_var_ids[rank0_1_var_ids %in% names(pbsGeno_raw)]))
#     
#     pbsGeno_subset <- pbsGeno_subset %>%
#       rename(id=IID)
#     
#     pbsGeno_joined <- pbsGeno_subset %>%
#       left_join(pc5_pbs_cov, by = "id")
#     pbs_geno_joined_cov <- t(pbsGeno_joined)
#     colnames(pbs_geno_joined_cov) <- pbs_geno_joined_cov[1,]
#     pbs_geno_joined_cov <- pbs_geno_joined_cov[-1,]
#     pbs_geno_joined_cov <- cbind(rownames(pbs_geno_joined_cov),pbs_geno_joined_cov)
#     colnames(pbs_geno_joined_cov)[1] <- c("id")
#     cov_dir <- paste0("output/01.qtltools_re/conditional_covariates/pbs_cov","/","covariates_PC5_",phe_id,"_rank1")
#     write.table(pbs_geno_joined_cov,
#                 file=cov_dir ,
#                 sep='\t',quote=F,row.names=F,col.names=T)
#     norm_dir <- paste0("output/01.qtltools_re/cond_norm_pbs/",phe_id,"_rank1")
#     perm_dir <- paste0("output/01.qtltools_re/cond_perm_pbs/",phe_id,"_rank1")
#     system(paste0("scripts/conditional_sqtl/qtltools_rerun.sh ",
#                   "output/geno/pbs_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz ",
#                   bed_dir, " ",
#                   norm_dir, " ", 
#                   perm_dir, " ",
#                   bed_txt, " ",
#                   cov_dir, " ",
#                   phe_id))
#     #rank2
#     rank0_2_var_ids <- temp_table %>%
#       filter(rank != 2) %>%
#       pull(var_id)
#     phe_id <- temp_table$phe_id[1]
#     pbsGeno_subset <- pbsGeno_raw %>%
#       select(c("IID", rank0_2_var_ids[rank0_2_var_ids %in% names(pbsGeno_raw)]))
#     
#     pbsGeno_subset <- pbsGeno_subset %>%
#       rename(id=IID)
#     
#     pbsGeno_joined <- pbsGeno_subset %>%
#       left_join(pc5_pbs_cov, by = "id")
#     pbs_geno_joined_cov <- t(pbsGeno_joined)
#     colnames(pbs_geno_joined_cov) <- pbs_geno_joined_cov[1,]
#     pbs_geno_joined_cov <- pbs_geno_joined_cov[-1,]
#     pbs_geno_joined_cov <- cbind(rownames(pbs_geno_joined_cov),pbs_geno_joined_cov)
#     colnames(pbs_geno_joined_cov)[1] <- c("id")
#     cov_dir <- paste0("output/01.qtltools_re/conditional_covariates/pbs_cov","/","covariates_PC5_",phe_id,"_rank2")
#     write.table(pbs_geno_joined_cov,
#                 file=cov_dir ,
#                 sep='\t',quote=F,row.names=F,col.names=T)
#     norm_dir <- paste0("output/01.qtltools_re/cond_norm_pbs/",phe_id,"_rank2")
#     perm_dir <- paste0("output/01.qtltools_re/cond_perm_pbs/",phe_id,"_rank2")
#     system(paste0("scripts/conditional_sqtl/qtltools_rerun.sh ",
#                   "output/geno/pbs_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz ",
#                   bed_dir, " ",
#                   norm_dir, " ", 
#                   perm_dir, " ",
#                   bed_txt, " ",
#                   cov_dir, " ",
#                   phe_id))
#     
#     #rank3
#     rank0_2_var_ids <- temp_table %>%
#       filter(rank != 3) %>%
#       pull(var_id)
#     phe_id <- temp_table$phe_id[1]
#     pbsGeno_subset <- pbsGeno_raw %>%
#       select(c("IID", rank0_2_var_ids[rank0_2_var_ids %in% names(pbsGeno_raw)]))
#     
#     pbsGeno_subset <- pbsGeno_subset %>%
#       rename(id=IID)
#     
#     pbsGeno_joined <- pbsGeno_subset %>%
#       left_join(pc5_pbs_cov, by = "id")
#     pbs_geno_joined_cov <- t(pbsGeno_joined)
#     colnames(pbs_geno_joined_cov) <- pbs_geno_joined_cov[1,]
#     pbs_geno_joined_cov <- pbs_geno_joined_cov[-1,]
#     pbs_geno_joined_cov <- cbind(rownames(pbs_geno_joined_cov),pbs_geno_joined_cov)
#     colnames(pbs_geno_joined_cov)[1] <- c("id")
#     cov_dir <- paste0("output/01.qtltools_re/conditional_covariates/pbs_cov","/","covariates_PC5_",phe_id,"_rank3")
#     write.table(pbs_geno_joined_cov,
#                 file=cov_dir ,
#                 sep='\t',quote=F,row.names=F,col.names=T)
#     norm_dir <- paste0("output/01.qtltools_re/cond_norm_pbs/",phe_id,"_rank3")
#     perm_dir <- paste0("output/01.qtltools_re/cond_perm_pbs/",phe_id,"_rank3")
#     system(paste0("scripts/conditional_sqtl/qtltools_rerun.sh ",
#                   "output/geno/pbs_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz ",
#                   bed_dir, " ",
#                   norm_dir, " ", 
#                   perm_dir, " ",
#                   bed_txt, " ",
#                   cov_dir, " ",
#                   phe_id))
#   }
# }




process_phe <- function(phe_id, rank_ids, suffix) {
  bed_txt <- paste0("output/01.qtltools_re/conditional_covariates/pbs_cov/", phe_id, ".txt")
  write.table(phe_id, file=bed_txt, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  pbsGeno_subset <- pbsGeno_raw %>%
    dplyr::select(c("IID", rank_ids[rank_ids %in% names(pbsGeno_raw)])) %>%
    dplyr::rename(id = IID) %>%
    left_join(pc5_pbs_cov, by = "id") %>%
    { t(.) } %>%
    `colnames<-`(.[1, ]) %>%
    .[-1, ] %>%
    cbind(rownames(.), .) %>%
    { `colnames<-`(. , c("id", colnames(.)[-1])) }
  
  cov_dir <- paste0("output/01.qtltools_re/conditional_covariates/pbs_cov/covariates_PC5_", phe_id, "_", suffix)
  write.table(pbsGeno_subset, file=cov_dir, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
  
  norm_dir <- paste0("output/01.qtltools_re/cond_norm_pbs/", phe_id, "_", suffix)
  perm_dir <- paste0("output/01.qtltools_re/cond_perm_pbs/", phe_id, "_", suffix)
  system(paste0("scripts/conditional_sqtl/qtltools_rerun.sh ", 
                "output/geno/pbs_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz ", 
                bed_dir, " ", norm_dir, " ", perm_dir, " ", bed_txt, " ", cov_dir, " ", phe_id))
}


for (i in seq_len(nrow(cond_pbs_filtered))) {
  print(paste("Processing row:", i))
  
  # Filter the main dataframe based on current `phe_id`
  temp_table <- cond_pbs %>%
    dplyr::filter(phe_id %in% cond_pbs_filtered[i, "phe_id"])  # Assuming 'phe_id' is the column of interest
  
  # Define 'bed_dir' if it varies per iteration or ensure it's already defined
  bed_dir <- paste0("output/gtex_cluster/qtltools_prep/ctrvsfnf_qqnorm_", temp_table$phe_chr[1], ".bed.gz")
  
  # Proceed with the function call if there are at least two rows to process
  if (nrow(temp_table) >= 2) {
    for (current_rank in 1:(nrow(temp_table) - 1)) {
      rank_var_ids <- temp_table %>%
        dplyr::filter(rank != current_rank) %>%
        pull(var_id)
      process_phe(temp_table$phe_id[1], rank_var_ids, paste0("rank", current_rank))
    }
  }
}


for (i in seq_len(nrow(cond_pbs_filtered))) {
  print(paste("Processing row:", i))
  
  # Filter the main dataframe based on current `phe_id`
  temp_table <- cond_pbs %>%
    dplyr::filter(phe_id %in% cond_pbs_filtered[i, "phe_id"])  # Assuming 'phe_id' is the column of interest
  
  # Define 'bed_dir' if it varies per iteration or ensure it's already defined
  bed_dir <- paste0("output/gtex_cluster/qtltools_prep/ctrvsfnf_qqnorm_", temp_table$phe_chr[1], ".bed.gz")
  
  # Proceed with the function call if there are at least two rows to process
  if (nrow(temp_table) >= 2) {
    current_rank  <- 0
      rank_var_ids <- temp_table %>%
        dplyr::filter(rank != current_rank) %>%
        pull(var_id)
      process_phe(temp_table$phe_id[1], rank_var_ids, paste0("rank", current_rank))
    
  }
}

#-------------------------------------------------------------------------------
#fnf
cond_fnf <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/conditional_fnf/conditional_fnf_top_variants.txt")
colnames(cond_fnf) <- header_cond
fnfGeno_raw <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/fnf_geno/06.subset_sigSNps/recodeA_fnf.raw")
colnames(fnfGeno_raw) <- sub("_.*", "", colnames(fnfGeno_raw))

# First, check if each phe_id includes ranks other than 0
phe_id_with_multiple_ranks <- cond_fnf %>%
  group_by(phe_id) %>%
  summarize(has_multiple_ranks = any(rank != 0)) %>%
  dplyr::filter(has_multiple_ranks)

# Filter original data to include only those IDs
cond_fnf_filtered <- cond_fnf %>%
  dplyr::filter(phe_id %in% phe_id_with_multiple_ranks$phe_id ) %>%
  dplyr::select(phe_id) %>%
  unique()



process_phe <- function(phe_id, rank_ids, suffix) {
  bed_txt <- paste0("output/01.qtltools_re/conditional_covariates/fnf_cov/", phe_id, ".txt")
  write.table(phe_id, file=bed_txt, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  fnfGeno_subset <- fnfGeno_raw %>%
    dplyr::select(c("IID", rank_ids[rank_ids %in% names(fnfGeno_raw)])) %>%
    dplyr::rename(id = IID) %>%
    left_join(pc4_fnf_cov, by = "id") %>%
    { t(.) } %>%
    `colnames<-`(.[1, ]) %>%
    .[-1, ] %>%
    cbind(rownames(.), .) %>%
    { `colnames<-`(. , c("id", colnames(.)[-1])) }
  
  cov_dir <- paste0("output/01.qtltools_re/conditional_covariates/fnf_cov/covariates_PC4_", phe_id, "_", suffix)
  write.table(fnfGeno_subset, file=cov_dir, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
  
  norm_dir <- paste0("output/01.qtltools_re/cond_norm_fnf/", phe_id, "_", suffix)
  perm_dir <- paste0("output/01.qtltools_re/cond_perm_fnf/", phe_id, "_", suffix)
  system(paste0("scripts/conditional_sqtl/qtltools_rerun.sh ", 
                "output/geno/fnf_geno/04.final/snpfiltered_wCHR_no_indels.final.recode.vcf.gz ", 
                bed_dir, " ", norm_dir, " ", perm_dir, " ", bed_txt, " ", cov_dir, " ", phe_id))
}


for (i in seq_len(nrow(cond_fnf_filtered))) {
  print(paste("Processing row:", i))
  
  # Filter the main dataframe based on current `phe_id`
  temp_table <- cond_fnf %>%
    dplyr::filter(phe_id %in% cond_fnf_filtered[i, "phe_id"])  # Assuming 'phe_id' is the column of interest
  
  # Define 'bed_dir' if it varies per iteration or ensure it's already defined
  bed_dir <- paste0("output/gtex_cluster/qtltools_prep/ctrvsfnf_qqnorm_", temp_table$phe_chr[1], ".bed.gz")
  
  # Proceed with the function call if there are at least two rows to process
  if (nrow(temp_table) >= 2) {
    for (current_rank in 1:(nrow(temp_table) - 1)) {
      rank_var_ids <- temp_table %>%
        dplyr::filter(rank != current_rank) %>%
        pull(var_id)
      process_phe(temp_table$phe_id[1], rank_var_ids, paste0("rank", current_rank))
    }
  }
}


for (i in seq_len(nrow(cond_fnf_filtered))) {
  print(paste("Processing row:", i))
  
  # Filter the main dataframe based on current `phe_id`
  temp_table <- cond_fnf %>%
    dplyr::filter(phe_id %in% cond_fnf_filtered[i, "phe_id"])  # Assuming 'phe_id' is the column of interest
  
  # Define 'bed_dir' if it varies per iteration or ensure it's already defined
  bed_dir <- paste0("output/gtex_cluster/qtltools_prep/ctrvsfnf_qqnorm_", temp_table$phe_chr[1], ".bed.gz")
  
  # Proceed with the function call if there are at least two rows to process
  if (nrow(temp_table) >= 2) {
    current_rank  <- 0
    rank_var_ids <- temp_table %>%
      dplyr::filter(rank != current_rank) %>%
      pull(var_id)
    process_phe(temp_table$phe_id[1], rank_var_ids, paste0("rank", current_rank))
    
  }
}

