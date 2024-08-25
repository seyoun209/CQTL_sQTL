#Finding the sQTL within 50kb of GWAS lead signals
setwd("/work/users/s/e/seyoun/CQTL_sQTL")
library(plotgardener)
library(ggplot2)
library(grid)
library(data.table)
library(dplyr)
library(stringr)
library(GenomicRanges)
source("scripts/utils/utils.R")


# Finding the overlaps for the coloc  SQTL and GWAS


# First, Finding the overlaps between GWAS lead signals within sQTL.

pbs_sig_qtl_cond_annot <- readRDS("output/01.qtltools_re/pbs_cond_sig_beta_maf_added.rds")
fnf_sig_qtl_cond_annot <- readRDS("output/01.qtltools_re/fnf_cond_sig_beta_maf_added.rds")

pbs_sig_qtl_cond_annot[, leadsnp := process_var_id(var_id)]
fnf_sig_qtl_cond_annot[, leadsnp := process_var_id(var_id)]

load(file="output/coloc/ld0_leadsnps/ld0_combined_leadsnps.rda")
ld0_qtl <- combined_ld_data
OAsubtypes <- c("AllOA", "FingerOA", "HandOA", "HipOA", "KneeHipOA", "KneeOA", "THR", "ThumbOA", "TJR", "TKR")


norm_header <- c("phe_id", "phe_chr", "phe_from", "phe_to", "phe_strd", 
                 "n_var_in_cis", "dist_phe_var", "var_id", "var_chr", "var_from",
                 "var_to", "nom_pval", "r_squared", "slope", "slope_se", "best_hit") 


output_dir <- "output/results_plots/coloc/50kb"
if (!dir.exists(output_dir)) { 
  dir.create(output_dir)
}  


for (subtype in OAsubtypes) {
  print(subtype)
  gwas_data <- fread(paste0("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg38_LD/", subtype, "/leads/EUR_", subtype, "_leads_ld_final.csv"))
  #gwas_data <- fread(paste0("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg38_LD/", subtype, "/leads/ALL_", subtype, "_leads_ld_final.csv"))
  subtype_QTL_ld05sharedGWAS_subset <- data.table()
  subtype_GWAS_sharedQTL_subset <- data.table()
  pdf(file = paste0("output/results_plots/coloc/50kb/distance_1mb",subtype,".pdf"), width = 5, height = 6.3)
  for (chr in c(1:22)) {
    print(chr)
    gwas_chr_data <- gwas_data %>% filter(chrom == chr)
    
    if(isEmpty(gwas_chr_data)) {
      next
    }
    combined_ld_qtl <- combined_ld_data[[paste0("chr", chr)]]
    
    gwas_leadSig <- gwas_chr_data$`CHR:hg38POS` |> unique() |> tibble(value =  `unique(gwas_chr_data$\`CHR:hg38POS\`)`) |> 
      dplyr::select(-c(`unique(gwas_chr_data$\`CHR:hg38POS\`)`)) %>%
      dplyr::mutate(
        chr = paste0("chr",str_split(value, ":", simplify = TRUE)[,1]),
        pos = as.numeric(str_split(value, ":", simplify = TRUE)[,2])
      ) |>
      dplyr::mutate(min_pos = pos -500000,
                      max_pos = pos + 500000) # This is for the 50k

    
    gwas_gr <- GRanges(
      seqnames =gwas_leadSig$chr,
      ranges = IRanges(start = gwas_leadSig$min_pos, end = gwas_leadSig$max_pos),
      value = gwas_leadSig$value
    )
    
    fnf_gr <- GRanges(
      seqnames = fnf_sig_qtl_cond_annot$var_chr,
      ranges = IRanges(start = fnf_sig_qtl_cond_annot$var_from, end = fnf_sig_qtl_cond_annot$var_to),
      var_id = fnf_sig_qtl_cond_annot$var_id,
      snp_pos = fnf_sig_qtl_cond_annot$leadsnp,
      phe_id = fnf_sig_qtl_cond_annot$phe_id
    )
    pbs_gr <- GRanges(
      seqnames = pbs_sig_qtl_cond_annot$var_chr,
      ranges = IRanges(start = pbs_sig_qtl_cond_annot$var_from, end = pbs_sig_qtl_cond_annot$var_to),
      var_id = pbs_sig_qtl_cond_annot$var_id,
      snp_pos = pbs_sig_qtl_cond_annot$leadsnp,
      phe_id = pbs_sig_qtl_cond_annot$phe_id
    )
    
    #Find overlap PBS
    pbs_overlapPairs <- findOverlapPairs(pbs_gr,gwas_gr)
    pbs_overlap_snp <- first(pbs_overlapPairs)$snp_pos |> unique()
    gwas_related_to_pbs <- second(pbs_overlapPairs)[which(first(pbs_overlapPairs)$snp_pos %in% pbs_overlap_snp)]
    gwas_related_pbs_snp <- second(pbs_overlapPairs)$value |> unique()
    
    #Find overlap FNF
    fnf_overlapPairs <- findOverlapPairs(fnf_gr,gwas_gr)
    fnf_overlap_snp <- first(fnf_overlapPairs)$snp_pos |> unique()
    gwas_related_fnf_snp <- second(pbs_overlapPairs)$value |> unique()
    gwas_related_to_fnf <- second(fnf_overlapPairs)[which(first(fnf_overlapPairs)$snp_pos %in% fnf_overlap_snp)]
    
    gwas_related_qtl <- unique(gwas_related_pbs_snp, gwas_related_fnf_snp)
    
    for (i in gwas_related_qtl) {
      print(i)
      qtls_snpID <- unique(first(pbs_overlapPairs)[which(second(pbs_overlapPairs)$value %in% i)]$snp_pos,
                           first(fnf_overlapPairs)[which(second(fnf_overlapPairs)$value %in% i)]$snp_pos)
      
      chr_gwas_leadsig <- first(pbs_overlapPairs)[which(second(pbs_overlapPairs)$value %in% i)][1] |> seqnames() |>as.factor()
      
      gwas_subset_leadsig_df <- gwas_chr_data |> dplyr::filter(`CHR:hg38POS` == i ) 
      if(nrow(gwas_subset_leadsig_df) == 1 ){
        next
      }
      gwas_df_subset <- gwas_subset_leadsig_df |>
        dplyr::filter(`CHR:hg38POS` %in% i )|>
        dplyr::filter(!grepl("NA", `ldbuddy_CHR:hg38POS`)) |>
        dplyr::select("ldbuddy_CHR:hg38POS", "ldbuddy_R2", "CHR:hg38POS", "rsID", "p") |>
        dplyr::mutate(
          chr_ldbuddy = paste0("chr", str_split(`ldbuddy_CHR:hg38POS`, ":", simplify = TRUE)[, 1]),
          pos_ldbuddy = as.numeric(str_split(`ldbuddy_CHR:hg38POS`, ":", simplify = TRUE)[, 2]))
      min_region <- gwas_leadSig |> dplyr::filter(value==i) |> dplyr::pull(min_pos)
      max_region <- gwas_leadSig |> dplyr::filter(value==i) |> dplyr::pull(max_pos)
      

            for (k in qtls_snpID) {
            
            pheid_pbs <- pbs_sig_qtl_cond_annot |> dplyr::filter(leadsnp %in% k) |> dplyr::select("phe_id", "SYMBOL")
            pheid_fnf <- fnf_sig_qtl_cond_annot |> dplyr::filter(leadsnp %in% k) |> dplyr::select("phe_id", "SYMBOL")
            ld0_qtl_chr <- ld0_qtl[[c(matrix(chr_gwas_leadsig))]] %>% dplyr::filter(leadsnp == k) 
            ld0_qtl_chr_range <- ld0_qtl_chr %>% dplyr::mutate( chr_ldbuddy = paste0("chr",str_split(ld0_qtl_chr$ldbuddy,":",simplify=TRUE)[,1]),
                                                                pos_ldbuddy = as.double(str_split(ld0_qtl_chr$ldbuddy,":",simplify=TRUE)[,2]))
            
            if (isEmpty(pheid_fnf)){
              cat("PBS_only \n")
              pbs_norm_qtl_chr <- fread(paste0("output/nominal_1mb/nominal_pbs/pc5/", c(matrix(chr_gwas_leadsig)), ".pbs.cis"), 
                                        select = c(1, 8, 9, 10, 12)) 
              colnames(pbs_norm_qtl_chr) <- norm_header[c(1, 8, 9, 10, 12)]
              pbs_qtl_region <- pbs_norm_qtl_chr |> dplyr::filter(phe_id %in% pheid_pbs$phe_id)
              pbs_qtl_region <- copy(pbs_qtl_region)
              pbs_qtl_region[, snpchr_loc := process_var_id(pbs_qtl_region$var_id)]
              
              fnf_norm_qtl_chr <- fread(paste0("output/nominal_1mb/nominal_fnf/pc4/", c(matrix(chr_gwas_leadsig)), ".fnf.cis"), 
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
              leftjoin_pbs <- inner_join(pbs_qtl_region, ld0_qtl_chr_range, by = c("snpchr_loc" = "ldbuddy"))
              leftjoin_fnf <- inner_join(fnf_qtl_region, ld0_qtl_chr_range, by = c("snpchr_loc" = "ldbuddy"))
              
              
            } else if(nrow(pheid_pbs) > 0 & nrow(pheid_fnf) > 0) {
              
    
        # Load only necessary columns and filter by relevant phe_id values
        pbs_norm_qtl_chr <- fread(paste0("output/nominal_1mb/nominal_pbs/pc5/", chr_gwas_leadsig, ".pbs.cis"), 
                                  select = c(1, 8, 9, 10, 12)) 
        colnames(pbs_norm_qtl_chr) <- norm_header[c(1, 8, 9, 10, 12)]
        pbs_qtl_region <- pbs_norm_qtl_chr |> dplyr::filter(phe_id %in% pheid_pbs$phe_id)
        pbs_qtl_region <- copy(pbs_qtl_region)
        pbs_qtl_region[, snpchr_loc := process_var_id(pbs_qtl_region$var_id)]
        
        fnf_norm_qtl_chr <- fread(paste0("output/nominal_1mb/nominal_fnf/pc4/", chr_gwas_leadsig, ".fnf.cis"), 
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
              cat("FNF only \n")
              fnf_norm_qtl_chr <- fread(paste0("output/nominal_1mb/nominal_fnf/pc4/", chr_gwas_leadsig, ".fnf.cis"), 
                                        select = c(1, 8, 9, 10, 12))
              colnames(fnf_norm_qtl_chr) <- norm_header[c(1, 8, 9, 10, 12)]
              fnf_qtl_region <- fnf_norm_qtl_chr |> dplyr::filter(phe_id %in% pheid_fnf$phe_id)
              fnf_qtl_region <- copy(fnf_qtl_region)
              fnf_qtl_region[, snpchr_loc := process_var_id(fnf_qtl_region$var_id)]
              
              pbs_norm_qtl_chr <- fread(paste0("output/nominal_1mb/nominal_pbs/pc5/", chr_gwas_leadsig, ".pbs.cis"), 
                                        select = c(1, 8, 9, 10, 12)) 
              colnames(pbs_norm_qtl_chr) <- norm_header[c(1, 8, 9, 10, 12)]
              pbs_qtl_region <- pbs_norm_qtl_chr |> dplyr::filter(phe_id %in% pheid_pbs$phe_id)
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
              
              leftjoin_pbs <- inner_join(pbs_qtl_region, ld0_qtl_chr_range, by = c("snpchr_loc" = "ldbuddy"))
              leftjoin_fnf <- inner_join(fnf_qtl_region, ld0_qtl_chr_range, by = c("snpchr_loc" = "ldbuddy"))
              
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
            chrom <- c(matrix(chr_gwas_leadsig))
            
            pbs_lead <- qtl_fi_plot$leadsnp  |> unique()
            fnf_lead <- fnf_fi_plot$leadsnp  |> unique()
            
            if(isEmpty(pheid_fnf$phe_id)){
              pheid_fnf$phe_id <- 'PBS_specific sQTL'
            }
            
            
            if(isEmpty(pheid_pbs$phe_id)){
              pheid_pbs$phe_id <- 'FNF_specific sQTL'
            }
            
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
            subtype <- subtype
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
              label = "-log10(p-value)", x = 0.24, y = 4.7, rot = 90,
              fontsize = 8, just = "center",
              default.units = "inches", fontfamily = "Helvetica"
            )
            plotText(label = "Fn-f", x = 0.6, y = 3.8, just = c("left", "top"),
                     fontfamily = "Helvetica", fontsize = 11)
            
            plotText(label = pheid_fnf$phe_id, x =1 , y =3.8, just = c("left", "top"),
                     fontfamily = "Helvetica", fontsize = 8)
            
            
            pbs_gene <- c(pheid_pbs$SYMBOL, pheid_fnf$SYMBOL) |> unique()
            plotgenes <- plotGenes(params = gwas_region, y = 5.5,
                                   height = 0.5,
                                   geneHighlights = data.frame("gene" = pbs_gene,
                                                               "color" = "#37a7db"))
            annoGenomeLabel(plot = plotgenes, params = gwas_region, fontsize = 8,y=6.05)
            
      
            }
      
    }
  }
  dev.off()
}

    
    
    



