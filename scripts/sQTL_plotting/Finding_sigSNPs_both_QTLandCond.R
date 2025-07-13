#subset all the sigSNPs_all_list 
setwd("/work/users/s/e/seyoun/CQTL_sQTL/")
library(data.table)

# This is for the only main QTL that have tested it
header <- read.table("scripts/sQTL_rscripts/header.txt")
new_header <- cbind(header, "qval","therehod")
sigQTL_fnf <- read.table("output/01.qtltools_re/01.significant/fnf_0.05_pc4.significant.txt") |> as.data.frame()
sigQTL_pbs <- read.table("output/01.qtltools_re/01.significant/pbs_0.05_pc5.significant.txt") |> as.data.frame()
colnames(sigQTL_pbs) <- new_header
colnames(sigQTL_fnf) <- new_header

#This is for the conditional QTL that have tested it. 
pbs_sig_qtl_cond_annot <- readRDS("output/01.qtltools_re/pbs_cond_sig_beta_maf_added.rds")
fnf_sig_qtl_cond_annot <- readRDS("output/01.qtltools_re/fnf_cond_sig_beta_maf_added.rds")

# Now, I need to subset only for the not the first signals. 
# Make new ID column 
pbs_sig_qtl_conditional_ID <- pbs_sig_qtl_cond_annot %>%
  dplyr::mutate(id = paste(phe_id, var_id, sep = "_")) 
fnf_sig_qtl_conditional_ID <- fnf_sig_qtl_cond_annot %>%
  dplyr::mutate(id = paste(phe_id, var_id, sep = "_")) 

sigQTL_pbs_ID <- sigQTL_pbs %>%
  dplyr::mutate(id = paste(phe_id, var_id, sep = "_")) 
sigQTL_fnf_ID <- sigQTL_fnf %>%
  dplyr::mutate(id = paste(phe_id, var_id, sep = "_")) 

# Filter rows where both phe_id and var_id don't match sigQTL_pbs
pbs_conditionalOnly <- pbs_sig_qtl_conditional_ID %>%
  dplyr::filter(!(id %in% sigQTL_pbs_ID$id))

fnf_conditionalOnly <- fnf_sig_qtl_conditional_ID %>%
  dplyr::filter(!(id %in% sigQTL_fnf_ID$id))


pbs_sqtl_notmathced_conditional_rank0 <- sigQTL_pbs_ID %>%
  dplyr::filter(!(id %in% pbs_sig_qtl_conditional_ID$id))

fnf_sqtl_notmathced_conditional_rank0 <- sigQTL_fnf_ID %>%
  dplyr::filter(!(id %in% fnf_sig_qtl_conditional_ID$id))

# I have to make sig SNPs list for the multiple cases.
# First sQTL from QTLtools permuatationw

#To run the GATK, the text has to be ending with .list.
sigSNPs_from_QTL_only <- unique(c(sigQTL_pbs$var_id, sigQTL_fnf$var_id))
write.table(sigSNPs_from_QTL_only,file="/work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/sigSNPs_permQTL.list",
            col.names = F,
            row.names=F,
            quote=F)


# Second conditonal anlaysis ranks not 0 

pbs_ranks_not0_var_id <- pbs_conditionalOnly |> dplyr::filter(rank !=0 ) |> dplyr::pull(var_id)
fnf_ranks_not0_var_id <- fnf_conditionalOnly |> dplyr::filter(rank !=0 ) |> dplyr::pull(var_id)

ranks_not0 <- unique(pbs_ranks_not0_var_id,fnf_ranks_not0_var_id)
write.table(ranks_not0,file="/work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/sigSNPs_cond_ranksnot0Only.list",
            col.names = F,
            row.names=F,
            quote=F)

# all the list including the sQTLs and conditonal all ranks. 

all_snps_lists <- unique(c(pbs_sig_qtl_cond_annot$var_id, fnf_sig_qtl_cond_annot$var_id, ranks_not0, sigSNPs_from_QTL_only))
write.table(all_snps_lists,file="/work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/sigSNPs_all_includeCond.list",
            col.names = F,
            row.names=F,
            quote=F)

# This is only Conditional all cases. 
sigSNPs_all.list <- unique(pbs_sig_qtl_cond_annot$var_id,fnf_sig_qtl_cond_annot$var_id)
#sGene_all.list <- unique(pbs_sig_qtl_cond_annot$ensg,fnf_sig_qtl_cond_annot$ensg)
write.table(sigSNPs_all.list,file="/work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/sigSNPs_all_CondOnly.list",
            col.names = F,
            row.names=F,
            quote=F)



