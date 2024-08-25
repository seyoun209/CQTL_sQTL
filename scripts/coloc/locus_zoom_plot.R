#Making plot for the locus zoom
setwd("/work/users/s/e/seyoun/CQTL_sQTL")
library(plotgardener)
library(ggplot2)
library(grid)
library(data.table)
library(dplyr)
library(stringr)
source("scripts/utils/utils.R")
#library(plotgardenerData)
#data("hg19_insulin_GWAS")

process_var_id <- function(var_ids) {
  var_ids_no_chr <- str_remove(var_ids, "^chr")
  split_elements <- str_split_fixed(var_ids_no_chr, ":", 3)
  processed_ids <- paste(split_elements[, 1], split_elements[, 2], sep = ":")
  return(processed_ids)
}


#preprocessing Input for the GWAS and sQTL
#load(file="output/coloc/gwas_summary_stat_list.rda") #name of the file gwas_summary_stats_list
OAsubtypes <- c("AllOA", "FingerOA", "HandOA", "HipOA", "KneeHipOA", "KneeOA", "THR", "ThumbOA", "TJR", "TKR")
# For The LD < 0.5 there are only 

norm_header <- c("phe_id", "phe_chr", "phe_from", "phe_to", "phe_strd", 
                 "n_var_in_cis", "dist_phe_var", "var_id", "var_chr", "var_from",
                 "var_to", "nom_pval", "r_squared", "slope", "slope_se", "best_hit") 


# load qtl_ld0 and sqtl_pvalue only for the right chromosome
load(file="output/coloc/ld0_leadsnps/ld0_combined_leadsnps.rda")
ld0_qtl <- combined_ld_data

pbs_sig_qtl_cond_annot <- readRDS("output/01.qtltools_re/pbs_cond_sig_beta_maf_added.rds")
fnf_sig_qtl_cond_annot <- readRDS("output/01.qtltools_re/fnf_cond_sig_beta_maf_added.rds")

pbs_sig_qtl_cond_annot[, leadsnp := process_var_id(var_id)]
fnf_sig_qtl_cond_annot[, leadsnp := process_var_id(var_id)]



# load qtl_ld0 and sqtl_pvalue
i <- OAsubtypes[4]

output_dir <- "output/results_plots/coloc"
if (!dir.exists(output_dir)) { 
  dir.create(output_dir)
}


for (i in OAsubtypes) {
  print(i)
  subtype_fi <- list.files("./output/coloc/gwas_qtl_ld05_subset",pattern=i,full.names = T)
  if(length(subtype_fi) == 0 ){
    print(paste0("no-subtype: ",i))
    next
  }
  gwas_df <- fread(subtype_fi[1])
  gwas_leadSig <- gwas_df$`CHR:hg38POS` |> unique() |> tibble(value =  `unique(gwas_df$\`CHR:hg38POS\`)`) |> 
    dplyr::select(-c(`unique(gwas_df$\`CHR:hg38POS\`)`)) %>%
    dplyr::mutate(
      chr = str_split(value, ":", simplify = TRUE)[,1],
      pos = str_split(value, ":", simplify = TRUE)[,2]
    )
  #---------------------------------------------------------------------------
  # qtl
  qtl_df <- fread(subtype_fi[2])
  leadsnps <- unique(qtl_df$leadsnp) |> tibble(value =`unique(qtl_df$leadsnp)`) |> 
    dplyr::select(-c(`unique(qtl_df$leadsnp)`)) %>%
    dplyr::mutate(
      chr = str_split(value, ":", simplify = TRUE)[,1],
      pos = str_split(value, ":", simplify = TRUE)[,2]
    )
  
  pdf(file = paste0("output/results_plots/coloc/ld05_distance_1mb_",i,".pdf"), width = 5, height = 6.3)
  
  for(j in gwas_leadSig$value){
    print(j)
    gwas_leadsig_pos <- gwas_df |> dplyr::filter(`CHR:hg38POS` == j ) |> dplyr::select("chrom","hg38pos") |>  slice_head(n = 1)
    gwas_df_subset <- gwas_df |>
      dplyr::filter(`CHR:hg38POS` %in% j )|>
      dplyr::filter(!grepl("NA", `ldbuddy_CHR:hg38POS`)) |>
      dplyr::select("ldbuddy_CHR:hg38POS", "ldbuddy_R2", "CHR:hg38POS", "rsID", "p") |>
      dplyr::mutate(
        chr_ldbuddy = paste0("chr", str_split(`ldbuddy_CHR:hg38POS`, ":", simplify = TRUE)[, 1]),
        pos_ldbuddy = as.numeric(str_split(`ldbuddy_CHR:hg38POS`, ":", simplify = TRUE)[, 2])
      )
    
    # Get signal region of + or - 500kb
    min_region <- gwas_leadsig_pos$hg38pos - 500000
    max_region <- gwas_leadsig_pos$hg38pos + 500000
    matching_leadsnps <- leadsnps |> dplyr::filter(chr == gwas_leadSig |> dplyr::filter(value == j) |> pull(chr))
    for (k in matching_leadsnps$value) {
      print(k)
      #-------------------------------------------------------------------------
      #QTL
      chr_leadsnps <- paste0("chr",str_split_fixed(k, ":",n=2)[1])
      
      
      ld0_qtl_chr <- ld0_qtl[[chr_leadsnps]] %>% dplyr::filter(leadsnp == k) 
      ld0_qtl_chr_range <- ld0_qtl_chr %>% dplyr::mutate( chr_ldbuddy = paste0("chr",str_split(ld0_qtl_chr$ldbuddy,":",simplify=TRUE)[,1]),
                                                          pos_ldbuddy = as.double(str_split(ld0_qtl_chr$ldbuddy,":",simplify=TRUE)[,2]))
      pbs_check_leadsnp_exist <- pbs_sig_qtl_cond_annot %>% dplyr::filter(leadsnp %in% k)
      fnf_check_leadsnp_exist <- fnf_sig_qtl_cond_annot %>% dplyr::filter(leadsnp %in% k)
      
      pheid_fnf <- NULL
      pheid_pbs <- NULL
      fnf_qtl_region <- NULL
      pbs_qtl_region <- NULL
      if(nrow(pbs_check_leadsnp_exist ) == 0){
        cat("fnf_check_leadsnp_exist contains the matching leadsnp values.\n")
        pheid_fnf <- pbs_check_leadsnp_exist[which.max(abs(fnf_check_leadsnp_exist$delta_beta)),] |> dplyr::select("phe_id","SYMBOL")
        
        fnf_norm_qtl_chr <- fread(paste0("output/nominal_1mb/nominal_fnf/pc4/", chr_leadsnps, ".fnf.cis"), 
                                  select = c(1, 8, 9, 10, 12))
        colnames(fnf_norm_qtl_chr) <- norm_header[c(1, 8, 9, 10, 12)]
        fnf_qtl_region <- fnf_norm_qtl_chr |> dplyr::filter(phe_id %in% pheid_fnf$phe_id)
        fnf_qtl_region <- copy(fnf_qtl_region)
        fnf_qtl_region[, snpchr_loc := process_var_id(fnf_qtl_region$var_id)]
        
        pbs_norm_qtl_chr <- fread(paste0("output/nominal_1mb/nominal_pbs/pc5/", chr_leadsnps, ".pbs.cis"), 
                                  select = c(1, 8, 9, 10, 12)) 
        colnames(pbs_norm_qtl_chr) <- norm_header[c(1, 8, 9, 10, 12)]
        
        pbs_qtl_region <- pbs_norm_qtl_chr %>% dplyr::filter(phe_id %in% pheid_fnf$phe_id) 
        pbs_qtl_region <- copy(pbs_qtl_region)
        pbs_qtl_region[, snpchr_loc := process_var_id(pbs_qtl_region$var_id)]
        
        pbs_qtl_region <- pbs_qtl_region %>%
          group_by(snpchr_loc) %>%
          dplyr::slice(which.min(nom_pval)) %>%
          ungroup()
        
        fnf_qtl_region <- fnf_qtl_region %>%
          group_by(snpchr_loc) %>%
          dplyr::slice(which.min(nom_pval)) %>%
          ungroup()
        
        leftjoin_pbs<- inner_join(pbs_qtl_region,ld0_qtl_chr_range, by=c("snpchr_loc"="ldbuddy"))
        leftjoin_fnf <- inner_join(fnf_qtl_region,ld0_qtl_chr_range, by=c("snpchr_loc"="ldbuddy"))
        
      } else if (nrow(pbs_check_leadsnp_exist ) >  0 & nrow(fnf_check_leadsnp_exist ) >  0) {
        # Identify relevant phe_id values
        pheid_pbs <- pbs_check_leadsnp_exist[which.max(abs(pbs_check_leadsnp_exist$delta_beta)), ] |> dplyr::select("phe_id", "SYMBOL")
        pheid_fnf <- fnf_check_leadsnp_exist[which.max(abs(fnf_check_leadsnp_exist$delta_beta)), ] |> dplyr::select("phe_id", "SYMBOL")
        
        # Load only necessary columns and filter by relevant phe_id values
        pbs_norm_qtl_chr <- fread(paste0("output/nominal_1mb/nominal_pbs/pc5/", chr_leadsnps, ".pbs.cis"), 
                                  select = c(1, 8, 9, 10, 12)) 
        colnames(pbs_norm_qtl_chr) <- norm_header[c(1, 8, 9, 10, 12)]
        pbs_qtl_region <- pbs_norm_qtl_chr |> dplyr::filter(phe_id %in% pheid_pbs$phe_id)
        pbs_qtl_region <- copy(pbs_qtl_region)
        pbs_qtl_region[, snpchr_loc := process_var_id(pbs_qtl_region$var_id)]
        
        fnf_norm_qtl_chr <- fread(paste0("output/nominal_1mb/nominal_fnf/pc4/", chr_leadsnps, ".fnf.cis"), 
                                  select = c(1, 8, 9, 10, 12))
        colnames(fnf_norm_qtl_chr) <- norm_header[c(1, 8, 9, 10, 12)]
        fnf_qtl_region <- fnf_norm_qtl_chr |> dplyr::filter(phe_id %in% pheid_fnf$phe_id)
        fnf_qtl_region <- copy(fnf_qtl_region)
        fnf_qtl_region[, snpchr_loc := process_var_id(fnf_qtl_region$var_id)]
        
        pbs_qtl_region <- pbs_qtl_region %>%
          group_by(snpchr_loc) %>%
          dplyr::slice(which.min(nom_pval)) %>%
          ungroup()
        
        fnf_qtl_region <- fnf_qtl_region %>%
          group_by(snpchr_loc) %>%
          dplyr::slice(which.min(nom_pval)) %>%
          ungroup()
        
        
        leftjoin_pbs <- inner_join(pbs_qtl_region, ld0_qtl_chr_range, by = c("snpchr_loc" = "ldbuddy"))
        leftjoin_fnf <- inner_join(fnf_qtl_region, ld0_qtl_chr_range, by = c("snpchr_loc" = "ldbuddy"))

        
      } else {
        cat("pbs_check_leadsnp_exist contains the matching leadsnp values.\n")
          pheid_pbs <- pbs_check_leadsnp_exist[which.max(abs(pbs_check_leadsnp_exist$delta_beta)),] |> dplyr::select("phe_id","SYMBOL")
          pbs_norm_qtl_chr <- fread(paste0("output/nominal_1mb/nominal_pbs/pc5/", chr_leadsnps, ".pbs.cis"), 
                                    select = c(1, 8, 9, 10, 12)) 
          colnames(pbs_norm_qtl_chr) <- norm_header[c(1, 8, 9, 10, 12)]
          pbs_qtl_region <- pbs_norm_qtl_chr |> dplyr::filter(phe_id %in% pheid_pbs$phe_id)
          pbs_qtl_region <- copy(pbs_qtl_region)
          pbs_qtl_region[, snpchr_loc := process_var_id(pbs_qtl_region$var_id)]
          
          
          fnf_norm_qtl_chr <- fread(paste0("output/nominal_1mb/nominal_fnf/pc4/", chr_leadsnps, ".fnf.cis"), 
                                    select = c(1, 8, 9, 10, 12))
          colnames(fnf_norm_qtl_chr) <- norm_header[c(1, 8, 9, 10, 12)]
          
          fnf_qtl_region <- fnf_norm_qtl_chr |> dplyr::filter(phe_id %in% pheid_pbs$phe_id)
          fnf_qtl_region <- copy(fnf_qtl_region)
          fnf_qtl_region[, snpchr_loc := process_var_id(fnf_qtl_region$var_id)]
          
          pbs_qtl_region <- pbs_qtl_region %>%
            group_by(snpchr_loc) %>%
            dplyr::slice(which.min(nom_pval)) %>%
            ungroup()
          
          fnf_qtl_region <- fnf_qtl_region %>%
            group_by(snpchr_loc) %>%
            dplyr::slice(which.min(nom_pval)) %>%
            ungroup()
        
        leftjoin_pbs<- inner_join(pbs_qtl_region,ld0_qtl_chr_range, by=c("snpchr_loc"="ldbuddy"))
        leftjoin_fnf <- inner_join(fnf_qtl_region,ld0_qtl_chr_range, by=c("snpchr_loc"="ldbuddy"))
      
        }
      
      
      #-------------------------------------------------------------------------
      #Plotting data
      gwas_df_group <- gwas_df_subset |>  mutate(LDgrp = cut(ldbuddy_R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))
      gwas_df_group$LDgrp <- addNA(gwas_df_group$LDgrp)
      gwas_fi_plot <- gwas_df_group |> dplyr::rename(chrom ="chr_ldbuddy",
                                                     pos = "pos_ldbuddy",
                                                     p="p",
                                                     LD ="ldbuddy_R2",
                                                     snp ="ldbuddy_CHR:hg38POS")  |> data.frame()
      gwas_fi_plot <- gwas_fi_plot %>%
        dplyr::select(chrom, pos, p, snp, LD,"CHR.hg38POS","rsID","LDgrp")
      


        qtl_fi_group <- leftjoin_pbs |>  mutate(LDgrp = cut(ldbuddy_R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))
        qtl_fi_group$LDgrp <- addNA(qtl_fi_group$LDgrp)
        
        qtl_fi_plot <- qtl_fi_group |> dplyr::rename(chrom ="chr_ldbuddy",
                                                     pos = "pos_ldbuddy",
                                                     p="nom_pval",
                                                     LD ="ldbuddy_R2",
                                                     snp ="snpchr_loc") |> data.frame()
        
        qtl_fi_plot <- qtl_fi_plot %>%
          dplyr::select(chrom, pos, p, snp, LD,"leadsnp","var_chr", "var_from","var_id","LDgrp")
        
        
        fnf_fi_group <- leftjoin_fnf |>  mutate(LDgrp = cut(ldbuddy_R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))
        fnf_fi_group$LDgrp <- addNA(fnf_fi_group$LDgrp)
        
        fnf_fi_plot <- fnf_fi_group |> dplyr::rename(chrom ="chr_ldbuddy",
                                                     pos = "pos_ldbuddy",
                                                     p="nom_pval",
                                                     LD ="ldbuddy_R2",
                                                     snp ="snpchr_loc") |> data.frame()
        
        fnf_fi_plot <- fnf_fi_plot %>%
          dplyr::select(chrom, pos, p, snp, LD,"leadsnp","var_chr", "var_from","var_id","LDgrp")
        
        gwas_lead <- gwas_fi_plot$CHR.hg38POS |> unique() 
        chrom <- paste0("chr",gwas_leadsig_pos$chrom)
        
        pbs_lead <- qtl_fi_plot$leadsnp  |> unique()
        fnf_lead <- fnf_fi_plot$leadsnp  |> unique()
        
        if(isEmpty(pheid_fnf$phe_id)){
          pheid_fnf$phe_id <- pheid_pbs$phe_id
        }
        
        
        if(isEmpty(pheid_pbs$phe_id)){
          pheid_pbs$phe_id <- pheid_fnf$phe_id
        }
        
        if(pbs_lead == fnf_lead){
          
          #Visualize locus plot
          pageCreate(width = 5, height = 6.3, showGuides = F)
          
          ## Plot all GWAS data
          
          gwas_region <- pgParams(assembly = "hg38",chrom = chrom,
                                  chromstart = min_region,
                                  chromend = max_region,
                                  x = 0.55, width = 3.75)
          gwas_ylim <- ceiling(max(log10(gwas_fi_plot$p)*-1)) + 2
          
        
          
          gwas_locus <- plotManhattan(data = gwas_fi_plot ,
                                      params = gwas_region,
                                      range = c(0, gwas_ylim),
                                      fill = colorby("LDgrp",
                                                     palette = colorRampPalette(c("#262C74",
                                                                                  "#98CDED",
                                                                                  "#499A53",
                                                                                  "#EEA741",
                                                                                  "#DD3931",
                                                                                  "grey"))),
                                      y = 0.5, height = 1.5,
                                      snpHighlights = data.frame(snp = c(gwas_lead, pbs_lead),
                                                                 pch = c(24, 23),
                                                                 cex = c(0.75, 0.75),
                                                                 col = c("black", "black")))
          annoYaxis(plot = gwas_locus, at = seq(0, gwas_ylim, 2),
                    axisLine = TRUE, fontsize = 8)
          plotText(
            label = "-log10(p-value)", x = 0.24, y = 1.25, rot = 90,
            fontsize = 8, just = "center",
            default.units = "inches", fontfamily = "Helvetica"
          )
          plotText(label = "GWAS", x = 0.65, y = 0.5, just = c("left", "top"),
                   fontfamily = "Helvetica", fontsize = 11)
          subtype <- i
          plotText(label = subtype, x = 0.5, y = 0.15, just = c("left", "top"),
                   fontfamily = "Helvetica", fontsize = 11, fontface = "bold")

          
          plotLegend(legend = c("0.8 - 1.0",
                                "0.6 - 0.8",
                                "0.4 - 0.6",
                                "0.2 - 0.4",
                                "0.0 - 0.2"),
                     fill = c("#DD3931", "#EEA741", "#499A53","#98CDED","#262C74"),
                     x = 3.5, y = 0.5, width = 0.1, height = 0.4, border = FALSE,
                     fontsize = 6)
          
          
          ## PBS QTL data
          
         pbs_ylim <- ceiling(max(log10(qtl_fi_plot$p)*-1)) + 2
         fnf_ylim <- ceiling(max(log10(fnf_fi_plot$p)*-1)) + 2
        
         
         y_lim <- max(pbs_ylim, fnf_ylim)
          
          pbs_locus <- plotManhattan(data = qtl_fi_plot ,
                                     params = gwas_region,
                                     range = c(0, y_lim),
                                     fill = colorby("LDgrp",
                                                    palette = colorRampPalette(c("#262C74",
                                                                                 "#98CDED",
                                                                                 "#499A53",
                                                                                 "#EEA741",
                                                                                 "#DD3931",
                                                                                 "grey"))),
                                     y = 2.2, height = 1.5,
                                     snpHighlights = data.frame(snp = c(gwas_lead, pbs_lead),
                                                                pch = c(24, 23),
                                                                cex = c(0.75, 0.75),
                                                                col = c("black", "black")))
          
          
          annoYaxis(plot = pbs_locus, at = seq(0, y_lim, 2),
                    axisLine = TRUE, fontsize = 8)
          plotText(
            label = "-log10(p-value)", x = 0.24, y = 2.95, rot = 90,
            fontsize = 8, just = "center",
            default.units = "inches", fontfamily = "Helvetica"
          )
          plotText(label = "PBS", x = 0.6, y = 2.1, just = c("left", "top"),
                   fontfamily = "Helvetica", fontsize = 11)
          
          plotText(label = pheid_pbs$phe_id, x =1 , y =2.1, just = c("left", "top"),
                   fontfamily = "Helvetica", fontsize = 8)
          
          
          grid.points(x = 3.5, y = 2.4, default.units = "native", pch = 24,
                      size = unit(0.5, "char"))
          plotText(label = paste0(gwas_lead),
                   fontsize = 7, fontfamily = "Helvetica",
                   just = "left", x = 3.6, y = 2.4)
          
          grid.points(x = 3.5, y = 2.5, default.units = "native", pch = 23,
                      size = unit(0.5, "char"))
          plotText(label =pbs_lead,
                   fontsize = 7, fontfamily = "Helvetica",
                   just = "left", x = 3.6, y = 2.5)
          
          
          #FN-f locus zoom

          
          fnf_locus <- plotManhattan(data = fnf_fi_plot ,
                                     params = gwas_region,
                                     range = c(0, y_lim),
                                     fill = colorby("LDgrp",
                                                    palette = colorRampPalette(c("#262C74",
                                                                                 "#98CDED",
                                                                                 "#499A53",
                                                                                 "#EEA741",
                                                                                 "#DD3931",
                                                                                 "grey"))),
                                     y = 3.9, height = 1.5,
                                     snpHighlights = data.frame(snp = c(gwas_lead, pbs_lead),
                                                                pch = c(24, 23),
                                                                cex = c(0.75, 0.75),
                                                                col = c("black", "black")))
          
          
          annoYaxis(plot = fnf_locus, at = seq(0, y_lim, 2),
                    axisLine = TRUE, fontsize = 8)
          plotText(
            label = "-log10(p-value)", x = 0.24, y = 2.95, rot = 90,
            fontsize = 8, just = "center",
            default.units = "inches", fontfamily = "Helvetica"
          )
          plotText(label = "Fn-f", x = 0.6, y = 3.8, just = c("left", "top"),
                   fontfamily = "Helvetica", fontsize = 11)

          plotText(label = pheid_fnf$phe_id, x =1 , y =3.8, just = c("left", "top"),
                   fontfamily = "Helvetica", fontsize = 8)
          
  
          pbs_gene <- pheid_pbs$SYMBOL
          plotgenes <- plotGenes(params = gwas_region, y = 5.5,
                                 height = 0.5,
                                 geneHighlights = data.frame("gene" = pbs_gene,
                                                             "color" = "#37a7db"))
          annoGenomeLabel(plot = plotgenes, params = gwas_region, fontsize = 8,y=6.05)
          
          
        } else {
          #Visualize locus plot
          pageCreate(width = 5, height = 6.3, showGuides = F)
          
          ## Plot all GWAS data
          
          gwas_region <- pgParams(assembly = "hg38",chrom = chrom,
                                  chromstart = min_region,
                                  chromend = max_region,
                                  x = 0.55, width = 3.75)
          gwas_ylim <- ceiling(max(log10(gwas_fi_plot$p)*-1)) + 2
          
          gwas_locus <- plotManhattan(data = gwas_fi_plot ,
                                      params = gwas_region,
                                      range = c(0, gwas_ylim),
                                      fill = colorby("LDgrp",
                                                     palette = colorRampPalette(c("#262C74",
                                                                                  "#98CDED",
                                                                                  "#499A53",
                                                                                  "#EEA741",
                                                                                  "#DD3931",
                                                                                  "grey"))),
                                      y = 0.5, height = 1.5,
                                      snpHighlights = data.frame(snp = c(gwas_lead, pbs_lead,fnf_lead),
                                                                 pch = c(24, 23,23),
                                                                 cex = c(0.75, 0.75,0.75),
                                                                 col = c("black", "black","black")))
          annoYaxis(plot = gwas_locus, at = seq(0, gwas_ylim, 2),
                    axisLine = TRUE, fontsize = 8)
          plotText(
            label = "-log10(p-value)", x = 0.24, y = 1.25, rot = 90,
            fontsize = 8, just = "center",
            default.units = "inches", fontfamily = "Helvetica"
          )
          plotText(label = "GWAS", x = 0.65, y = 0.5, just = c("left", "top"),
                   fontfamily = "Helvetica", fontsize = 11)
          subtype <- i
          plotText(label = subtype, x = 0.5, y = 0.15, just = c("left", "top"),
                   fontfamily = "Helvetica", fontsize = 11, fontface = "bold")
          
          
          plotLegend(legend = c("0.8 - 1.0",
                                "0.6 - 0.8",
                                "0.4 - 0.6",
                                "0.2 - 0.4",
                                "0.0 - 0.2"),
                     fill = c("#DD3931", "#EEA741", "#499A53","#98CDED","#262C74"),
                     x = 3.5, y = 0.5, width = 0.1, height = 0.4, border = FALSE,
                     fontsize = 6)
          
          
          ## PBS QTL data
          
          pbs_ylim <- ceiling(max(log10(qtl_fi_plot$p)*-1)) + 2
          fnf_ylim <- ceiling(max(log10(fnf_fi_plot$p)*-1)) + 2
          
          
          y_lim <- max(pbs_ylim, fnf_ylim)
          
          pbs_locus <- plotManhattan(data = qtl_fi_plot ,
                                     params = gwas_region,
                                     range = c(0, y_lim),
                                     fill = colorby("LDgrp",
                                                    palette = colorRampPalette(c("#262C74",
                                                                                 "#98CDED",
                                                                                 "#499A53",
                                                                                 "#EEA741",
                                                                                 "#DD3931",
                                                                                 "grey"))),
                                     y = 2.2, height = 1.5,
                                     snpHighlights = data.frame(snp = c(gwas_lead, pbs_lead,fnf_lead),
                                                                pch = c(24, 23,23),
                                                                cex = c(0.75, 0.75,0.75),
                                                                col = c("black", "black","black")))
          
          
          annoYaxis(plot = pbs_locus, at = seq(0, y_lim, 2),
                    axisLine = TRUE, fontsize = 8)
          plotText(
            label = "-log10(p-value)", x = 0.24, y = 2.95, rot = 90,
            fontsize = 8, just = "center",
            default.units = "inches", fontfamily = "Helvetica"
          )
          plotText(label = "PBS", x = 0.6, y = 2.1, just = c("left", "top"),
                   fontfamily = "Helvetica", fontsize = 11)
          
          plotText(label = pheid_pbs$phe_id, x =1 , y =2.1, just = c("left", "top"),
                   fontfamily = "Helvetica", fontsize = 8)
          
          
          grid.points(x = 3.5, y = 2.4, default.units = "native", pch = 24,
                      size = unit(0.5, "char"))
          plotText(label = paste0(gwas_lead),
                   fontsize = 7, fontfamily = "Helvetica",
                   just = "left", x = 3.6, y = 2.4)
          
          grid.points(x = 3.5, y = 2.5, default.units = "native", pch = 23,
                      size = unit(0.5, "char"))
          plotText(label =pbs_lead,
                   fontsize = 7, fontfamily = "Helvetica",
                   just = "left", x = 3.6, y = 2.5)
          
          
          #FN-f locus zoom

          
          fnf_locus <- plotManhattan(data = fnf_fi_plot ,
                                     params = gwas_region,
                                     range = c(0, y_lim),
                                     fill = colorby("LDgrp",
                                                    palette = colorRampPalette(c("#262C74",
                                                                                 "#98CDED",
                                                                                 "#499A53",
                                                                                 "#EEA741",
                                                                                 "#DD3931",
                                                                                 "grey"))),
                                     y = 3.9, height = 1.5,
                                     snpHighlights = data.frame(snp = c(gwas_lead, pbs_lead,fnf_lead),
                                                                pch = c(24, 23,23),
                                                                cex = c(0.75, 0.75,0.75),
                                                                col = c("black", "black","black")))
          
          
          annoYaxis(plot = fnf_locus, at = seq(0, y_lim, 2),
                    axisLine = TRUE, fontsize = 8)
          plotText(
            label = "-log10(p-value)", x = 0.24, y = 2.95, rot = 90,
            fontsize = 8, just = "center",
            default.units = "inches", fontfamily = "Helvetica"
          )
          plotText(label = "Fn-f", x = 0.6, y = 3.8, just = c("left", "top"),
                   fontfamily = "Helvetica", fontsize = 11)
          
          plotText(label = pheid_fnf$phe_id, x =1 , y =3.8, just = c("left", "top"),
                   fontfamily = "Helvetica", fontsize = 8)
          
          grid.points(x = 3.5, y = 3.9, default.units = "native", pch = 24,
                      size = unit(0.5, "char"))
          plotText(label = paste0(gwas_lead),
                   fontsize = 7, fontfamily = "Helvetica",
                   just = "left", x = 3.6, y = 3.9)
          
          grid.points(x = 3.5, y = 4.0, default.units = "native", pch = 23,
                      size = unit(0.5, "char"))
          plotText(label =fnf_lead,
                   fontsize = 7, fontfamily = "Helvetica",
                   just = "left", x = 3.6, y = 4.0)
          
          pbs_gene <- pheid_pbs$SYMBOL
          fnf_gene <- pheid_fnf$SYMBOL
          plotgenes <- plotGenes(params = gwas_region, y = 5.5,
                                 height = 0.5,
                                 geneHighlights = data.frame("gene" = c(pbs_gene,fnf_gene),
                                                             "color" = "#37a7db"))
          annoGenomeLabel(plot = plotgenes, params = gwas_region, fontsize = 8,y=6.05)
          
          
        }
    }
  }
  dev.off()
}







