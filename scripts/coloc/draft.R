#chr6:18399699:18401957:clu_31772_+
#chr6:33073669:33080680:clu_30945_-
#chr16:69925484:69929448:clu_14156_+
#chr12:123226351:123227463:clu_8761_-
#chr12:122981266:122981595:clu_9868_+
#chr3:52658342:52662133:clu_25608_-
#chr22:37808163:37810027:clu_24955_+

keep_PP4 <- function(x){
  keep(x, ~.$summary["PP.H4.abf"] >= 0.7)
}
signal = x$signal,
PP4 = x$summary[["PP.H4.abf"]]))

load("output/coloc/coloc_results_pbs_THR_.rda")
thr_names <- c(str_split(names(coloc_results),"_",4, simplify = TRUE)[,4])
response_pbs_results |> dplyr::filter(phe_id %in% c("chr6:18399699:18401957:clu_31772_+", "chr6:33073669:33080680:clu_30945_-", "chr16:69925484:69929448:clu_14156_+",
                            "chr12:123226351:123227463:clu_8761_-","chr12:122981266:122981595:clu_9868_+",thr_names,"chr21:39347292:39347379:clu_24187_-"))


response_pbs_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_pbs_re_lmer.rds")
response_fnf_results <- readRDS("output/01.qtltools_re/conditional_fnf/response_fnf_re_lmer.rds")
##Filter the significant and adding beta value and MAF
threshold <- as.numeric(0.05)


sig_PBS_interactionResults <- response_pbs_results %>% 
  dplyr::filter(interaction_pval <= threshold)

pbs_pvals <- sapply(response_pbs_results, function(x) x$pval)
final_annotated_pbs$annova_pval <- pbs_pvals

response_pbs_ensg <- response_pbs_results %>% 
  dplyr::filter(interaction_pval < 0.05) %>%
  dplyr::filter(abs(delta_beta) >=0.2)
pull(ensg) %>%
  unique()

fnf_pvals <- sapply(response_fnf_results, function(x) x$pval)
final_annotated_fnf$annova_pval <- fnf_pvals

response_conf_fnf_interection <- final_annotated_fnf %>% 
  dplyr::filter(annova_pval < 0.05) %>%
  dplyr::filter(abs(Beta_diff) >=0.2)
pull(ensg) %>%
  unique() 

