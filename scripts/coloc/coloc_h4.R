setwd("/work/users/s/e/seyoun/CQTL_sQTL/")


coloc_pbs_results <- list.files("output/coloc",pattern="*_.rda",full.names = TRUE)

coloc_re <- data.frame()
for(i in c(1:length(coloc_pbs_results))){
  coloc_h4_df <- data.frame()
 
  load(coloc_pbs_results[i])
  coloc_name_df <- do.call(rbind,str_split(names(coloc_results),"_",5))
  if(is.null(coloc_name_df)){next}
  coloc_re_length <- length(names(coloc_results))
  colnames(coloc_name_df) <- c("subtype","gwas_sig","gwasrsid","qtl_sig","phe_id")
  for(j in c(1:coloc_re_length)){
    coloc_h4 <- coloc_results[[j]]$summary[["PP.H4.abf"]] 
   coloc_h4_df <- rbind(coloc_h4_df,coloc_h4)
  }
  colnames(coloc_h4_df) <- c('H4')
  coloc_df <- cbind(coloc_name_df,coloc_h4_df)
  coloc_re <- rbind(coloc_re,coloc_df)
  
}


coloc_re_sig <- coloc_re |> dplyr::filter(H4 > 0.7)
