#subset all the sigSNPs_all_list 
setwd("/work/users/s/e/seyoun/CQTL_sQTL/")
library(data.table)

pbs_sig_qtl_cond_annot <- readRDS("output/01.qtltools_re/pbs_cond_sig_beta_maf_added.rds")
fnf_sig_qtl_cond_annot <- readRDS("output/01.qtltools_re/fnf_cond_sig_beta_maf_added.rds")

sigSNPs_all.list <- unique(pbs_sig_qtl_cond_annot$var_id,fnf_sig_qtl_cond_annot$var_id)
write.table(sigSNPs_all.list,file="/work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/sigSNPs_all.list",
            col.names = F,
            row.names=F,
            quote=F)



