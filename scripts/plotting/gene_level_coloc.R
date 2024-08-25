# Gene level
library(dplyr)
library(tidyverse)
library(plotgardener)
library(data.table)
library(grid)
library(RColorBrewer)
library(ggtext)
library(GenomicFeatures)
library(ggpubr)
library(gt)
library(bstfun)
setwd("/work/users/s/e/seyoun/CQTL_sQTL")
yl_gn_bu <- brewer.pal(n = 9, name = "YlGnBu")
# Signal data 
dir_merged <- "/work/users/s/e/seyoun/CQTL_sQTL/output/signals/merged_norm/"
ctl_signal <- paste0(dir_merged,"CTL_norm.bw")
fnf_signal <- paste0(dir_merged,"FNF_norm.bw")
oa_signal <- paste0(dir_merged,"OA_norm.bw")


load("/work/users/s/e/seyoun/Ref/genome/gencode.v45.annotation.TxDb")

txdb <- loadDb("../crispr/02.test_seq/gencode.v45.annotation.TxDb")
txdb_genes <- genes(txdb)

remove_version <- function(ensg) {
  gsub("\\.\\d+$", "", ensg)
}

# Prepare txdb_genes data
txdb_genes_df <- as.data.frame(txdb_genes) %>%
  mutate(gene_id = remove_version(gene_id)) %>%
  dplyr::select(gene_id, seqnames, start, end)

norm_header <- c("phe_id", "phe_chr", "phe_from", "phe_to", "phe_strd", 
                 "n_var_in_cis", "dist_phe_var", "var_id", "var_chr", "var_from",
                 "var_to", "nom_pval", "r_squared", "slope", "slope_se", "best_hit")


pbs_norm_qtl_pc5_list <- readRDS("output/nominal_1mb/nominal_pbs/pbs_norm_qtl_pc5_list.rds")
fnf_norm_qtl_pc4_list <- readRDS("output/nominal_1mb/nominal_fnf/fnf_norm_qtl_pc4_list.rds")

load("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/freq/maf_id_add.rds") #maf_id_add  is the name
maf_subset_rsID <- maf_id_add |> dplyr::select(-c("minor_allele","MAF"))

#make a plot for the response QTL
response_pbs_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_pbs_results.rds") 
response_fnf_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_fnf_results.rds")
#This is ratios and meta data
ctl_fnf_ratio <- fread("output/clu_fnf/ratio_fnf.txt") |> as.data.frame()
rownames(ctl_fnf_ratio) <- ctl_fnf_ratio$Junction
ratios_fnf <- ctl_fnf_ratio[,-1]
intron_test <- c(response_pbs_results$phe_id, response_fnf_results$phe_id) |> unique()
psi_matrix <- ratios_fnf[rownames(ratios_fnf) %in% intron_test,] |> as.data.frame()
intronID <- rownames(psi_matrix)
psi_ratio <- apply(psi_matrix, 2, as.numeric)
rownames(psi_ratio) <- intronID
psi_ratio.df <- t(psi_ratio) |> as.data.frame()

clustersIDs <- str_split_fixed(rownames(ratios_fnf) ,":",4)[,4]

pbsGeno_raw <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/pbs_geno/06.subset_sigSNps/recodeA_pbs.raw")
colnames(pbsGeno_raw) <- sub("_.*", "", colnames(pbsGeno_raw))

fnfGeno_raw <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/fnf_geno/06.subset_sigSNps/recodeA_fnf.raw")
colnames(fnfGeno_raw) <- sub("_.*", "", colnames(fnfGeno_raw))

#oaGeno_raw <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/oa_geno/06.subset_sigSNps/recodeA_OA.raw")
#colnames(oaGeno_raw) <- sub("_.*", "", colnames(oaGeno_raw))

all_geno_matrix <- bind_rows(pbsGeno_raw,fnfGeno_raw)
all_geno_transpose_df <- all_geno_matrix %>%
  column_to_rownames(var = "IID")  %>%
  dplyr::select(-FID, -MAT,-PAT, -SEX, -PHENOTYPE) %>% 
  as.data.frame() %>%  # Ensure it's still a data frame
t() 
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




#Plotting for per all the coloc FNF genes
coloc_results_fnf <- fread("output/coloc/coloc_result_table_fnf.txt")
coloc_results_pbs <- fread("output/coloc/coloc_result_table.txt")
sig_fnf_coloc <- coloc_results_fnf |> dplyr::filter(coloc_H4 >= 0.7) |> distinct(phe_id,sQTL_rsID, .keep_all = TRUE)
sig_fnf_coloc <- coloc_results_fnf |> dplyr::filter(coloc_H4+coloc_H3 >= 0.7) |> distinct(phe_id,sQTL_rsID, .keep_all = TRUE)
sig_pbs_coloc <- coloc_results_pbs |> dplyr::filter(coloc_H4 >= 0.7) |> distinct(phe_id,sQTL_rsID, .keep_all = TRUE)
sig_pbs_coloc <- coloc_results_pbs |> dplyr::filter(coloc_H4+coloc_H3 >= 0.7) |> distinct(phe_id,sQTL_rsID, .keep_all = TRUE)
# sig_fnf_coloc_df <- cbind(sig_fnf_coloc, 
#       str_split(sig_fnf_coloc$phe_id ,":",simplify = TRUE) %>%
#   data.frame() |>
#   setNames(c("phe_chr","phe_start", "phe_end", "clusterID")))



#-------------------------------------------------------------------------------
pdf(file = paste0("output/results_plots/coloc/coloc_sig/fnf_significant_p4p3_0.7.pdf"), width = 4.5, height = 8.0)
for(i in 1:nrow(sig_fnf_coloc)){
fnf_results <- response_fnf_results |> dplyr::filter(phe_id %in% sig_fnf_coloc$phe_id[i])
#intronID <- response_pbs_results |> dplyr::filter(SYMBOL == gene_name) |> dplyr::pull("phe_id")
minorAllele <- unique(fnf_results$minor_allele)

psi_matrix <- psi_ratio[rownames(psi_ratio) %in% fnf_results$phe_id, ]
psi_data <- data.frame(sampleID = names(psi_matrix), psi = as.double(psi_matrix))

variantID <- unique(fnf_results$var_id)
rsID <- unique(fnf_results$rsID)
alleles <- unlist(strsplit(variantID, ":"))
ref_allele <- alleles[3]
alt_allele <- alleles[4]
protective_allele <- ifelse(minorAllele == ref_allele, alt_allele, ref_allele)

geno_data <- all_geno_transpose_df[rownames(all_geno_transpose_df) %in% variantID, ] %>%
  as.data.frame() %>%
  rownames_to_column(var = "sampleID")
colnames(geno_data)[2] <- "genotype"

meta_data_fi <- geno_data %>%
  left_join(meta_catl_simple,by = c("sampleID"="ID")) %>%
  left_join(psi_data,by = c("sampleID"="sampleID"))

meta_combined_all <-  meta_data_fi %>%
  mutate(
    genotype = factor(genotype, levels = c("0", "1", "2")),
    Condition = ifelse(Condition == "CTL", "PBS", ifelse(Condition == "FNF", "FN-f", NA)),
    Condition = factor(Condition, levels = c("PBS", "FN-f")),
    genotype = case_when(
      genotype == "2" ~ paste(minorAllele, minorAllele, sep = "/"),
      genotype == "1" ~ paste(minorAllele, protective_allele, sep = "/"),
      genotype == "0" ~ paste(protective_allele, protective_allele, sep = "/")
    )
  )  %>%
  mutate(
    genotype = factor(genotype, levels = c(
      paste(protective_allele, protective_allele, sep = "/"),
      paste(minorAllele, protective_allele, sep = "/"),
      paste(minorAllele, minorAllele, sep = "/")
    ))
  )


maxrange <- range(meta_combined_all$psi, na.rm = TRUE)
yaxis_range_min <- min(0,round(maxrange[1],1))
yaxis_range_max <- round(maxrange[2], 1)

if(yaxis_range_min == 0 & yaxis_range_max == 0){
  yaxis_range_min <-  maxrange[1]
  yaxis_range_max <-  maxrange[2]
}

genotype_plot <- ggplot(meta_combined_all, aes(x = genotype, y = psi, fill = Condition)) +
  geom_boxplot(outlier.shape = NA,
               linewidth = 0.25, alpha = 0.7) +
  geom_point(color = "grey40", position = position_jitterdodge(), size = 0.25) +
  geom_smooth(aes(group = Condition, color = Condition), method = "lm", se = FALSE, linetype = "solid") +
  labs(x = "Genotype", y = "Intron usage") +
  scale_fill_manual(values = c('#1e87a5','#FFB81C')) +
  scale_color_manual(values = c('#1e87a5','#FFB81C')) +
  #scale_y_continuous(limits = c(yaxis_range_min, yaxis_range_max), breaks = seq(yaxis_range_min, yaxis_range_max, length.out=5)) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(strip.placement = "outside",
        axis.line.y =  element_line(linewidth = 0.25),
        axis.line.x =  element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 0.25),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(size = 6, family= "Helvetica",
                                        margin = margin(r = -0.5)),
        text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.text.x = element_text(color = "black", size = 6, margin = margin(b = 0)),
        strip.background = element_blank(),
        strip.text = element_blank(),  # This removes the facet titles
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        legend.position = "none",  # This removes the legend
        panel.spacing.x = unit(0.5, "cm")) 
#facet_wrap(~ Condition, ncol = 2)



clustersIDs_fnf <- str_split_fixed(fnf_results$phe_id ,":",4)[,4]
psi_subset_clusters <- ratios_fnf[which(clustersIDs %in% clustersIDs_fnf),] |> t()
if(nrow(psi_subset_clusters)== 1){
  psi_subset_clusters <- ratios_fnf[which(clustersIDs == clustersIDs_fnf),] |> data.frame()
  colnames(psi_subset_clusters) <- fnf_results$phe_id
}
cluster_posInfo <- str_split(colnames(psi_subset_clusters),":",simplify = TRUE) |> 
  data.frame() |> setNames(c("cluster_chr","cluster_start","cluster_end","clusterID"))
psi_subset_clusters <- tibble::rownames_to_column(as.data.frame(psi_subset_clusters), var = "sampleID")
meta_combined_all_clusters <- left_join(meta_combined_all,psi_subset_clusters,by="sampleID")

table_expression <- meta_combined_all_clusters %>%
  group_by(genotype, Condition) %>%
  summarise(across(starts_with("chr"), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  ungroup() |> arrange(Condition, genotype)
table_expression_rounded <- table_expression %>%
  mutate(across(where(is.numeric), ~round(., 4)))


gg_table <- ggtexttable(table_expression_rounded, 
                        rows = NULL, 
                        theme = ttheme(
                          colnames.style = colnames_style(color = "black", fill = "grey", size = 7),
                          tbody.style = tbody_style(color = "black", fill = "white", size = 7)
                        )) +     theme_minimal() +
  theme(panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        
        panel.spacing.x = unit(0.5, "cm"))


ggsave(filename =   paste0("output/results_plots/coloc/coloc_sig/fnf_",fnf_results$phe_id,"_",fnf_results$SYMBOL,"_" ,
                           fnf_results$var_id,".pdf"), 
       plot =gg_table,width = 15, height = 5, units = "in")


chrom <- fnf_results$phe_chr
ld_file <- paste0("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld/",chrom,"/",variantID,".ld")
ld_calc <- fread(ld_file) |> dplyr::select(c("SNP_B","R2")) |>
  dplyr::rename("var_id" = "SNP_B")


pbs_norminal_qtl <- pbs_norm_qtl_pc5_list[[chrom]]
pbs_qtl_region <- pbs_norminal_qtl |> dplyr::filter(phe_id %in% fnf_results$phe_id)
pbs_qtl_pval_rsid <- left_join(pbs_qtl_region,maf_subset_rsID  ,by="var_id")

fnf_norminal_qtl <- fnf_norm_qtl_pc4_list[[chrom]]
fnf_qtl_region <- fnf_norminal_qtl |> dplyr::filter(phe_id %in% fnf_results$phe_id)
fnf_qtl_pval_rsid <- left_join(fnf_qtl_region,maf_subset_rsID  ,by="var_id")


leftjoin_pbs <- inner_join(pbs_qtl_pval_rsid, ld_calc, by = c("var_id" = "var_id"))
leftjoin_fnf <- inner_join(fnf_qtl_pval_rsid, ld_calc, by = c("var_id" = "var_id"))

#GWAS data subsetting ####
gwas_data <- fread(paste0("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg38_LD/", sig_fnf_coloc$subtype[i], "/leads/EUR_", sig_fnf_coloc$subtype[i], "_leads_ld_final.csv"))

gwas_data_subset <- gwas_data |> dplyr::filter(rsID  == sig_fnf_coloc$GWAS_rsID[i]) |>
  dplyr::filter(!grepl("NA", `ldbuddy_CHR:hg38POS`)) |>
  dplyr::select("ldbuddy_CHR:hg38POS", "ldbuddy_R2", "CHR:hg38POS", "ldbuddy_rsID","rsID","CHR:hg38POS", "p")|>
  dplyr::mutate(
    chr_ldbuddy = paste0("chr", str_split(`ldbuddy_CHR:hg38POS`, ":", simplify = TRUE)[, 1]),
    pos_ldbuddy = as.numeric(str_split(`ldbuddy_CHR:hg38POS`, ":", simplify = TRUE)[, 2])
  )


#-------------------------------------------------------------------------
#Plotting data

gwas_df_group <- gwas_data_subset |>  mutate(LDgrp = cut(ldbuddy_R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))
gwas_df_group$LDgrp <- addNA(gwas_df_group$LDgrp)
gwas_df_LDinclude <- gwas_df_group |> dplyr::rename(chrom ="chr_ldbuddy",
                                               pos = "pos_ldbuddy",
                                               p="p",
                                               LD ="ldbuddy_R2",
                                               snp ="ldbuddy_rsID")  |> data.frame()
gwas_fi_plot <- gwas_df_LDinclude %>%
  dplyr::select(chrom, pos, p, snp, LD,"ldbuddy_CHR.hg38POS","CHR.hg38POS","rsID","LDgrp")
gwas_rsID <- sig_fnf_coloc$GWAS_rsID[i]
gwas_leadsig_pos <- str_split(gwas_fi_plot$CHR.hg38POS,":",simplify = TRUE)[,2] |> unique() |> as.double()
gwas_min <-  gwas_leadsig_pos- 100000
gwas_max <- gwas_leadsig_pos + 100000
qtl_min <- fnf_results$var_from   -100000
qtl_max  <-  fnf_results$var_from +100000
minregion <- min(gwas_min, qtl_min)
maxregion <- max(gwas_max, qtl_max)

pbs_ld <- leftjoin_pbs |>  mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))
pbs_ld$LDgrp <- addNA(pbs_ld$LDgrp)

pbs_locus_plot <- pbs_ld |> dplyr::rename(chrom ="var_chr",
                                          pos = "var_from",
                                          p="nom_pval",
                                          LD ="R2",
                                          snp ="rsID") |> data.frame()

pbs_locus_plot <- pbs_locus_plot %>%
  dplyr::select("chrom", "pos", "p", "snp", "LD","LDgrp","phe_id","var_id")


fnf_ld <- leftjoin_fnf |>   mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))
fnf_ld$LDgrp <- addNA(fnf_ld$LDgrp)

fnf_locus_plot <- fnf_ld |>  dplyr::rename(chrom ="var_chr",
                                           pos = "var_from",
                                           p="nom_pval",
                                           LD ="R2",
                                           snp ="rsID") |> data.frame()


fnf_locus_plot <- fnf_locus_plot %>%
  dplyr::select("chrom", "pos", "p", "snp", "LD","LDgrp","phe_id","var_id")

#------------------------------------------------------------------------
#plotting

pageCreate(width = 4.5, height =7.5, showGuides = FALSE)

region_pg <- pgParams(assembly = "hg38",chrom = chrom,
                      chromstart = minregion,
                      chromend = maxregion,
                      x = 0.55, width = 3.5)
pbs_ylim <- ceiling(max(log10(pbs_locus_plot$p)*-1)) + 2
fnf_ylim <- ceiling(max(log10(fnf_locus_plot$p)*-1)) + 2
gwas_ylim <- ceiling(max(log10(gwas_fi_plot$p)*-1)) + 2
ylim_pg <- max(pbs_ylim,fnf_ylim)



gwas_locus <- plotManhattan(data = gwas_fi_plot ,
                            params = region_pg,
                            range = c(0, gwas_ylim),
                            fill = colorby("LDgrp",
                                           palette = colorRampPalette(c("#262C74",
                                                                        "#98CDED",
                                                                        "#499A53",
                                                                        "#EEA741",
                                                                        "#DD3931",
                                                                        "grey"))),
                            y = 0.5, height = 1,
                            snpHighlights = data.frame(snp = c(gwas_rsID, rsID),
                                                       pch = c(24, 23),
                                                       cex = c(0.75, 0.75),
                                                       col = c("black", "black")))
annoYaxis(plot = gwas_locus, at = seq(0, gwas_ylim, 2),
          axisLine = TRUE, fontsize = 8)
plotText(
  label = "-log10(p-value)", x = 0.24, y = 1, rot = 90,
  fontsize = 8, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = "GWAS", x = 0.65, y = 0.5, just = c("left", "top"),
         fontfamily = "Helvetica", fontsize = 8)


pbs_locus <- plotManhattan(data = pbs_locus_plot ,
                           params = region_pg,
                           range = c(0, ylim_pg),
                           fill = colorby("LDgrp",
                                          palette = colorRampPalette(c("#262C74",
                                                                       "#98CDED",
                                                                       "#499A53",
                                                                       "#EEA741",
                                                                       "#DD3931",
                                                                       "grey"))),
                           y = 1.5+0.1, height = 1,
                           snpHighlights = data.frame(snp = c(gwas_rsID, rsID),
                                                      pch = c(24, 23),
                                                      cex = c(0.75, 0.75),
                                                      col = c("black", "black")))
annoYaxis(plot = pbs_locus, at = seq(0, ylim_pg, 2),
          axisLine = TRUE, fontsize = 8)
plotText(
  label = "-log10(p-value)", x = 0.24, y = 2.2, rot = 90,
  fontsize = 8, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = "PBS sQTL", x = 0.65, y = 1.6, just = c("left", "top"),
         fontfamily = "Helvetica", fontsize = 8)


#plotLegend(legend = c("0.8 - 1.0",
#                      "0.6 - 0.8",
#                      "0.4 - 0.6",
#                      "0.2 - 0.4",
#                      "0.0 - 0.2"),
#           fill = c("#DD3931", "#EEA741", "#499A53","#98CDED","#262C74"),
#           x = 3.9, y = 2, width = 0.075, height = 0.3, border = FALSE,
#           fontsize = 6)

#fnf plot
fnf_locus <- plotManhattan(data = fnf_locus_plot ,
                           params = region_pg,
                           range = c(0, ylim_pg),
                           fill = colorby("LDgrp",
                                          palette = colorRampPalette(c("#262C74",
                                                                       "#98CDED",
                                                                       "#499A53",
                                                                       "#EEA741",
                                                                       "#DD3931",
                                                                       "grey"))),
                           y = 2.7, height = 1,
                           snpHighlights = data.frame(snp = c(gwas_rsID, rsID),
                                                      pch = c(24, 23),
                                                      cex = c(0.75, 0.75),
                                                      col = c("black", "black")))


annoYaxis(plot = fnf_locus, at = seq(0, ylim_pg, 2),
          axisLine = TRUE, fontsize = 8)
plotText(
  label = "-log10(p-value)", x = 0.24, y = 3.2, rot = 90,
  fontsize = 8, just = "center",
  default.units = "inches", fontfamily = "Helvetica"
)
plotText(label = "FN-f sQTL", x = 0.6, y = 2.7, just = c("left", "top"),
         fontfamily = "Helvetica", fontsize = 8)


grid.points(x = 2.7, y = 0.56, default.units = "native", pch = 24,
            size = unit(0.5, "char"))
grid.points(x = 2.7, y = 1.76, default.units = "native", pch = 23,
            size = unit(0.5, "char"))
grid.points(x = 2.7, y = 2.76, default.units = "native", pch = 23,
            size = unit(0.5, "char"))

plotText(label =  paste0(gwas_rsID, "-(" ,sig_fnf_coloc$subtype[i]," index)"),
         fontsize = 8, fontfamily = "Helvetica",
         just = "left", x = 2.8, y =0.55)
plotText(label = paste0(rsID, "-(sQTL index)"),
         fontsize = 8, fontfamily = "Helvetica",
         just = "left", x = 2.8, y = 1.75)
plotText(label = paste0(rsID, "-(sQTL index)"),
         fontsize = 8, fontfamily = "Helvetica",
         just = "left", x = 2.8, y = 2.75)


#RNA_signals <- plotMultiSignal(data = list(ctl_signal,
#                                           fnf_signal),
#                               params = region_pg,
#                               y = 3.85, height = 0.74, linecolor = c(yl_gn_bu[7],yl_gn_bu[7]), 
#                               fill = c(yl_gn_bu[7],yl_gn_bu[7]),
#                               default.units = "inches",
#                               gapdistance = 0.02)
#plotText(label = "PBS",
#         fontsize = 7, x = 0.5, y =3.9 , just = "left", fontfamily = "Helvetica")
#plotText(label = "FN-f",
#         fontsize = 7, x = 0.5, y =4.3, just = "left", fontfamily = "Helvetica")
#plotText(label = "OA",
#         fontsize = 7, x = 0.5, y =4.65, just = "left", fontfamily = "Helvetica")
#plotText(
#  label = "RNA", x = 0.35, y = 4.2, rot = 90,
#  fontsize = 8, just = "center",
#  default.units = "inches", fontfamily = "Helvetica"
#)

gene_hl <- fnf_results$SYMBOL
plotgenes <- plotGenes(params = region_pg, y = 3.8,
                       height = 0.5,
                       geneHighlights = data.frame("gene" = gene_hl,
                                                   "color" = "#37a7db"))

annoGenomeLabel(plot = plotgenes, params = region_pg, fontsize = 8,y=4.35)




# Process pbs_results
gene_info <- fnf_results %>%
  left_join(txdb_genes_df, by = c("ensg" = "gene_id")) %>%
  dplyr::rename(gene_chr = seqnames, gene_start = start, gene_end = end) |> 
  dplyr::select("gene_chr","gene_start","gene_end")


plot_selectedGene <- plotGenes( chrom =gene_info$gene_chr,
                                chromstart = gene_info$gene_start,
                                chromend = gene_info$gene_end,
                                assembly = "hg38",
                                y = 4.5,
                                x=0.55,
                                width=3.5,
                                height = 0.75, strandLabels = FALSE,fontsize = 0,
                                geneHighlights = data.frame("gene" = gene_hl,
                                                            "color" = "#37a7db"))



hl_start <-  min(cluster_posInfo$cluster_start) |>  as.numeric()
hl_end <-  max(cluster_posInfo$cluster_end) |> as.numeric()
anno_highlight_plot <- annoHighlight(
  plot = plot_selectedGene,
  chrom = gene_info$gene_chr, 
  chromstart = hl_start,
  chromend = hl_end,
  y = 4.7, height = 0.3, just = c("left", "top"),
  default.units = "inches"
)


plotGG(plot =  genotype_plot, x = 0.2, y = 5.25, height = 2.35, width = 2.35)

plotText(label =  fnf_results$phe_id, 
         x = 0.6,
         y=5.25,
         fontsize = 8, fontfamily = "Helvetica",
         just="left",
         fontcolor = "black")

plotText(label =  fnf_results$rsID, 
         x = 1.25,
         y=7.6,
         fontsize = 8, fontfamily = "Helvetica",
         just="left",
         fontcolor = "black")
}

dev.off()


#-------------------------------------------------------------------------------

#pbs

pdf(file = paste0("output/results_plots/coloc/coloc_sig/pbs_significant_p4p3_0.7.pdf"), width = 4.5, height = 8.0)
for(i in 1:nrow(sig_pbs_coloc)){
  pbs_results <- response_pbs_results |> dplyr::filter(phe_id %in% sig_pbs_coloc$phe_id[i])
  #intronID <- response_pbs_results |> dplyr::filter(SYMBOL == gene_name) |> dplyr::pull("phe_id")
  minorAllele <- unique(pbs_results$minor_allele)
  
  psi_matrix <- psi_ratio[rownames(psi_ratio) %in% pbs_results$phe_id, ]
  psi_data <- data.frame(sampleID = names(psi_matrix), psi = as.double(psi_matrix))
  
  variantID <- unique(pbs_results$var_id)
  rsID <- unique(pbs_results$rsID)
  alleles <- unlist(strsplit(variantID, ":"))
  ref_allele <- alleles[3]
  alt_allele <- alleles[4]
  protective_allele <- ifelse(minorAllele == ref_allele, alt_allele, ref_allele)
  
  geno_data <- all_geno_transpose_df[rownames(all_geno_transpose_df) %in% variantID, ] %>%
    as.data.frame() %>%
    rownames_to_column(var = "sampleID")
  colnames(geno_data)[2] <- "genotype"
  
  meta_data_fi <- geno_data %>%
    left_join(meta_catl_simple,by = c("sampleID"="ID")) %>%
    left_join(psi_data,by = c("sampleID"="sampleID"))
  
  meta_combined_all <-  meta_data_fi %>%
    mutate(
      genotype = factor(genotype, levels = c("0", "1", "2")),
      Condition = ifelse(Condition == "CTL", "PBS", ifelse(Condition == "FNF", "FN-f", NA)),
      Condition = factor(Condition, levels = c("PBS", "FN-f")),
      genotype = case_when(
        genotype == "2" ~ paste(minorAllele, minorAllele, sep = "/"),
        genotype == "1" ~ paste(protective_allele, minorAllele, sep = "/"),
        genotype == "0" ~ paste(protective_allele, protective_allele, sep = "/")
      )
    )  %>%
    mutate(
      genotype = factor(genotype, levels = c(
        paste(protective_allele, protective_allele, sep = "/"),
        paste(protective_allele, minorAllele, sep = "/"),
        paste(minorAllele, minorAllele, sep = "/")
      ))
    )
  
  
  maxrange <- range(meta_combined_all$psi, na.rm = TRUE)
  yaxis_range_min <- min(0,round(maxrange[1],1))
  yaxis_range_max <- round(maxrange[2], 1)
  
  if(yaxis_range_min == 0 & yaxis_range_max == 0){
    yaxis_range_min <-  maxrange[1]
    yaxis_range_max <-  maxrange[2]
  }
  
  genotype_plot <- ggplot(meta_combined_all, aes(x = genotype, y = psi, fill = Condition)) +
    geom_boxplot(outlier.shape = NA,
                 linewidth = 0.25, alpha = 0.7) +
    geom_point(color = "grey40", position = position_jitterdodge(), size = 0.25) +
    geom_smooth(aes(group = Condition, color = Condition), method = "lm", se = FALSE, linetype = "solid") +
    labs(x = "Genotype", y = "Intron usage") +
    scale_fill_manual(values = c('#1e87a5','#FFB81C')) +
    scale_color_manual(values = c('#1e87a5','#FFB81C')) +
    #scale_y_continuous(limits = c(yaxis_range_min, yaxis_range_max), breaks = seq(yaxis_range_min, yaxis_range_max, length.out=5)) +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(strip.placement = "outside",
          axis.line.y =  element_line(linewidth = 0.25),
          axis.line.x =  element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(color = "black", linewidth = 0.25),
          axis.ticks.length.y = unit(-0.1, "cm"),
          axis.title.x = element_blank(),
          axis.title.y = element_markdown(size = 6, family= "Helvetica",
                                          margin = margin(r = -0.5)),
          text = element_text(family = "Helvetica"),
          axis.text.y = element_text(color = "black", size = 6),
          axis.text.x = element_text(color = "black", size = 6, margin = margin(b = 0)),
          strip.background = element_blank(),
          strip.text = element_blank(),  # This removes the facet titles
          panel.background = element_rect(fill = "transparent", color = "transparent"),
          plot.background = element_rect(fill = "transparent", color = "transparent"),
          panel.grid = element_blank(),
          legend.position = "none",  # This removes the legend
          panel.spacing.x = unit(0.5, "cm")) 
    #facet_wrap(~ Condition, ncol = 2)
  
  
  
  clustersIDs_pbs <- str_split_fixed(pbs_results$phe_id ,":",4)[,4]
  psi_subset_clusters <- ratios_fnf[which(clustersIDs %in% clustersIDs_pbs),] |> t()
  if(nrow(psi_subset_clusters)== 1){
    psi_subset_clusters <- ratios_fnf[which(clustersIDs == clustersIDs_pbs),] |> data.frame()
    colnames(psi_subset_clusters) <- pbs_results$phe_id
  }
  cluster_posInfo <- str_split(colnames(psi_subset_clusters),":",simplify = TRUE) |> 
    data.frame() |> setNames(c("cluster_chr","cluster_start","cluster_end","clusterID"))
  psi_subset_clusters <- tibble::rownames_to_column(as.data.frame(psi_subset_clusters), var = "sampleID")
  meta_combined_all_clusters <- left_join(meta_combined_all,psi_subset_clusters,by="sampleID")
  
  table_expression <- meta_combined_all_clusters %>%
    group_by(genotype, Condition) %>%
    summarise(across(starts_with("chr"), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
    ungroup() |> arrange(Condition, genotype)
  table_expression_rounded <- table_expression %>%
    mutate(across(where(is.numeric), ~round(., 4)))
  

  gg_table <- ggtexttable(table_expression_rounded, 
                          rows = NULL, 
                          theme = ttheme(
                            colnames.style = colnames_style(color = "black", fill = "grey", size = 7),
                            tbody.style = tbody_style(color = "black", fill = "white", size = 7)
                          )) +     theme_minimal() +
    theme(panel.background = element_rect(fill = "transparent", color = "transparent"),
          plot.background = element_rect(fill = "transparent", color = "transparent"),
          panel.grid = element_blank(),

          panel.spacing.x = unit(0.5, "cm"))
  

  ggsave(filename =   paste0("output/results_plots/coloc/coloc_sig/pbs_",pbs_results$phe_id,"_",pbs_results$SYMBOL,"_" ,pbs_results$var_id,".pdf"), 
         plot =gg_table,width = 15, height = 5, units = "in")
  

  chrom <- pbs_results$phe_chr
  ld_file <- paste0("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/ld/",chrom,"/",variantID,".ld")
  ld_calc <- fread(ld_file) |> dplyr::select(c("SNP_B","R2")) |>
    dplyr::rename("var_id" = "SNP_B")
  
  
  pbs_norminal_qtl <- pbs_norm_qtl_pc5_list[[chrom]]
  pbs_qtl_region <- pbs_norminal_qtl |> dplyr::filter(phe_id %in% pbs_results$phe_id)
  pbs_qtl_pval_rsid <- left_join(pbs_qtl_region,maf_subset_rsID  ,by="var_id")
  
  fnf_norminal_qtl <- fnf_norm_qtl_pc4_list[[chrom]]
  fnf_qtl_region <- fnf_norminal_qtl |> dplyr::filter(phe_id %in% pbs_results$phe_id)
  fnf_qtl_pval_rsid <- left_join(fnf_qtl_region,maf_subset_rsID  ,by="var_id")
  
  
  leftjoin_pbs <- inner_join(pbs_qtl_pval_rsid, ld_calc, by = c("var_id" = "var_id"))
  leftjoin_fnf <- inner_join(fnf_qtl_pval_rsid, ld_calc, by = c("var_id" = "var_id"))
  
  #GWAS data subsetting ####
  gwas_data <- fread(paste0("/proj/phanstiel_lab/External/gwas/OA/Boer_2021_hg38_LD/", sig_pbs_coloc$subtype[i], "/leads/EUR_", sig_pbs_coloc$subtype[i], "_leads_ld_final.csv"))
  
  gwas_data_subset <- gwas_data |> dplyr::filter(rsID  == sig_pbs_coloc$GWAS_rsID[i]) |>
    dplyr::filter(!grepl("NA", `ldbuddy_CHR:hg38POS`)) |>
    dplyr::select("ldbuddy_CHR:hg38POS", "ldbuddy_R2", "CHR:hg38POS", "ldbuddy_rsID","rsID","CHR:hg38POS", "p")|>
    dplyr::mutate(
      chr_ldbuddy = paste0("chr", str_split(`ldbuddy_CHR:hg38POS`, ":", simplify = TRUE)[, 1]),
      pos_ldbuddy = as.numeric(str_split(`ldbuddy_CHR:hg38POS`, ":", simplify = TRUE)[, 2])
    )
  
  
  #-------------------------------------------------------------------------
  #Plotting data
  
  gwas_df_group <- gwas_data_subset |>  mutate(LDgrp = cut(ldbuddy_R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))
  gwas_df_group$LDgrp <- addNA(gwas_df_group$LDgrp)
  gwas_df_LDinclude <- gwas_df_group |> dplyr::rename(chrom ="chr_ldbuddy",
                                                      pos = "pos_ldbuddy",
                                                      p="p",
                                                      LD ="ldbuddy_R2",
                                                      snp ="ldbuddy_rsID")  |> data.frame()
  gwas_fi_plot <- gwas_df_LDinclude %>%
    dplyr::select(chrom, pos, p, snp, LD,"ldbuddy_CHR.hg38POS","CHR.hg38POS","rsID","LDgrp")
  gwas_rsID <- sig_pbs_coloc$GWAS_rsID[i]
  gwas_leadsig_pos <- str_split(gwas_fi_plot$CHR.hg38POS,":",simplify = TRUE)[,2] |> unique() |> as.double()
  gwas_min <-  gwas_leadsig_pos- 100000
  gwas_max <- gwas_leadsig_pos + 100000
  qtl_min <- pbs_results$var_from   -100000
  qtl_max  <-  pbs_results$var_from +100000
  minregion <- min(gwas_min, qtl_min)
  maxregion <- max(gwas_max, qtl_max)
  
  pbs_ld <- leftjoin_pbs |>  mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))
  pbs_ld$LDgrp <- addNA(pbs_ld$LDgrp)
  
  pbs_locus_plot <- pbs_ld |> dplyr::rename(chrom ="var_chr",
                                            pos = "var_from",
                                            p="nom_pval",
                                            LD ="R2",
                                            snp ="rsID") |> data.frame()
  
  pbs_locus_plot <- pbs_locus_plot %>%
    dplyr::select("chrom", "pos", "p", "snp", "LD","LDgrp","phe_id","var_id")
  
  
  fnf_ld <- leftjoin_fnf |>   mutate(LDgrp = cut(R2, c(0, 0.2, 0.4, 0.6, 0.8, 1)))
  fnf_ld$LDgrp <- addNA(fnf_ld$LDgrp)
  
  fnf_locus_plot <- fnf_ld |>  dplyr::rename(chrom ="var_chr",
                                             pos = "var_from",
                                             p="nom_pval",
                                             LD ="R2",
                                             snp ="rsID") |> data.frame()
  
  
  fnf_locus_plot <- fnf_locus_plot %>%
    dplyr::select("chrom", "pos", "p", "snp", "LD","LDgrp","phe_id","var_id")
  
  #------------------------------------------------------------------------
  #plotting
  
  pageCreate(width = 4.5, height =8, showGuides = FALSE)
  
  region_pg <- pgParams(assembly = "hg38",chrom = chrom,
                        chromstart = minregion,
                        chromend = maxregion,
                        x = 0.55, width = 3.5)
  pbs_ylim <- ceiling(max(log10(pbs_locus_plot$p)*-1)) + 2
  fnf_ylim <- ceiling(max(log10(fnf_locus_plot$p)*-1)) + 2
  gwas_ylim <- ceiling(max(log10(gwas_fi_plot$p)*-1)) + 2
  ylim_pg <- max(pbs_ylim,fnf_ylim)
  
  
  
  gwas_locus <- plotManhattan(data = gwas_fi_plot ,
                              params = region_pg,
                              range = c(0, gwas_ylim),
                              fill = colorby("LDgrp",
                                             palette = colorRampPalette(c("#262C74",
                                                                          "#98CDED",
                                                                          "#499A53",
                                                                          "#EEA741",
                                                                          "#DD3931",
                                                                          "grey"))),
                              y = 0.5, height = 1,
                              snpHighlights = data.frame(snp = c(gwas_rsID, rsID),
                                                         pch = c(24, 23),
                                                         cex = c(0.75, 0.75),
                                                         col = c("black", "black")))
  annoYaxis(plot = gwas_locus, at = seq(0, gwas_ylim, 2),
            axisLine = TRUE, fontsize = 8)
  plotText(
    label = "-log10(p-value)", x = 0.24, y = 1, rot = 90,
    fontsize = 8, just = "center",
    default.units = "inches", fontfamily = "Helvetica"
  )
  plotText(label = "GWAS", x = 0.65, y = 0.5, just = c("left", "top"),
           fontfamily = "Helvetica", fontsize = 8)
  
  
  pbs_locus <- plotManhattan(data = pbs_locus_plot ,
                             params = region_pg,
                             range = c(0, ylim_pg),
                             fill = colorby("LDgrp",
                                            palette = colorRampPalette(c("#262C74",
                                                                         "#98CDED",
                                                                         "#499A53",
                                                                         "#EEA741",
                                                                         "#DD3931",
                                                                         "grey"))),
                             y = 1.5+0.1, height = 1,
                             snpHighlights = data.frame(snp = c(gwas_rsID, rsID),
                                                        pch = c(24, 23),
                                                        cex = c(0.75, 0.75),
                                                        col = c("black", "black")))
  annoYaxis(plot = pbs_locus, at = seq(0, ylim_pg, 2),
            axisLine = TRUE, fontsize = 8)
  plotText(
    label = "-log10(p-value)", x = 0.24, y = 2.2, rot = 90,
    fontsize = 8, just = "center",
    default.units = "inches", fontfamily = "Helvetica"
  )
  plotText(label = "PBS sQTL", x = 0.65, y = 1.6, just = c("left", "top"),
           fontfamily = "Helvetica", fontsize = 8)
  
  
  #plotLegend(legend = c("0.8 - 1.0",
  #                      "0.6 - 0.8",
  #                      "0.4 - 0.6",
  #                      "0.2 - 0.4",
  #                      "0.0 - 0.2"),
  #           fill = c("#DD3931", "#EEA741", "#499A53","#98CDED","#262C74"),
  #           x = 3.9, y = 2, width = 0.075, height = 0.3, border = FALSE,
  #           fontsize = 6)
  
  #fnf plot
  fnf_locus <- plotManhattan(data = fnf_locus_plot ,
                             params = region_pg,
                             range = c(0, ylim_pg),
                             fill = colorby("LDgrp",
                                            palette = colorRampPalette(c("#262C74",
                                                                         "#98CDED",
                                                                         "#499A53",
                                                                         "#EEA741",
                                                                         "#DD3931",
                                                                         "grey"))),
                             y = 2.7, height = 1,
                             snpHighlights = data.frame(snp = c(gwas_rsID, rsID),
                                                        pch = c(24, 23),
                                                        cex = c(0.75, 0.75),
                                                        col = c("black", "black")))
  
  
  annoYaxis(plot = fnf_locus, at = seq(0, ylim_pg, 2),
            axisLine = TRUE, fontsize = 8)
  plotText(
    label = "-log10(p-value)", x = 0.24, y = 3.2, rot = 90,
    fontsize = 8, just = "center",
    default.units = "inches", fontfamily = "Helvetica"
  )
  plotText(label = "FN-f sQTL", x = 0.6, y = 2.7, just = c("left", "top"),
           fontfamily = "Helvetica", fontsize = 8)
  
  
  grid.points(x = 2.7, y = 0.56, default.units = "native", pch = 24,
              size = unit(0.5, "char"))
  grid.points(x = 2.7, y = 1.76, default.units = "native", pch = 23,
              size = unit(0.5, "char"))
  grid.points(x = 2.7, y = 2.76, default.units = "native", pch = 23,
              size = unit(0.5, "char"))
  
  plotText(label =  paste0(gwas_rsID, "-(" ,sig_pbs_coloc$subtype[i]," index)"),
           fontsize = 8, fontfamily = "Helvetica",
           just = "left", x = 2.8, y =0.55)
  plotText(label = paste0(rsID, "-(sQTL index)"),
           fontsize = 8, fontfamily = "Helvetica",
           just = "left", x = 2.8, y = 1.75)
  plotText(label = paste0(rsID, "-(sQTL index)"),
           fontsize = 8, fontfamily = "Helvetica",
           just = "left", x = 2.8, y = 2.75)
  
  
  #RNA_signals <- plotMultiSignal(data = list(ctl_signal,
  #                                           fnf_signal),
  #                               params = region_pg,
  #                               y = 3.85, height = 0.74, linecolor = c(yl_gn_bu[7],yl_gn_bu[7]), 
  #                               fill = c(yl_gn_bu[7],yl_gn_bu[7]),
  #                               default.units = "inches",
  #                               gapdistance = 0.02)
  #plotText(label = "PBS",
  #         fontsize = 7, x = 0.5, y =3.9 , just = "left", fontfamily = "Helvetica")
  #plotText(label = "FN-f",
  #         fontsize = 7, x = 0.5, y =4.3, just = "left", fontfamily = "Helvetica")
  #plotText(label = "OA",
  #         fontsize = 7, x = 0.5, y =4.65, just = "left", fontfamily = "Helvetica")
  #plotText(
  #  label = "RNA", x = 0.35, y = 4.2, rot = 90,
  #  fontsize = 8, just = "center",
  #  default.units = "inches", fontfamily = "Helvetica"
  #)
  
  gene_hl <- pbs_results$SYMBOL

  plotgenes <- plotGenes(params = region_pg, y = 3.8,
                         height = 0.5,
                        geneHighlights = data.frame("gene" = gene_hl,
                                                 "color" = "#37a7db"))
  
  annoGenomeLabel(plot = plotgenes, params = region_pg, fontsize = 8,y=4.35)
  

  
  
  # Process pbs_results
  gene_info <- pbs_results %>%
    left_join(txdb_genes_df, by = c("ensg" = "gene_id")) %>%
    dplyr::rename(gene_chr = seqnames, gene_start = start, gene_end = end) |> 
    dplyr::select("gene_chr","gene_start","gene_end")
  
  
plot_selectedGene <- plotGenes( chrom =gene_info$gene_chr,
                          chromstart = gene_info$gene_start,
                          chromend = gene_info$gene_end,
                          assembly = "hg38",
                          y = 4.5,
                          x=0.55,
                          width=3.5,
                         height = 0.75, strandLabels = FALSE,fontsize = 0,
                         geneHighlights = data.frame("gene" = gene_hl,
                         "color" = "#37a7db"))
#   
  # plot_selectedtranscript <- plotTranscripts(
  #   chrom =gene_info$gene_chr,
  #   chromstart = gene_info$gene_start,
  #   chromend = gene_info$gene_end,
  #   assembly = "hg38",
  #   y = 4.35,
  #   x=0.55,
  #   width=3.5,
  #  height = 1.5,
  #   just = c("left", "top"), default.units = "inches",
  #   #transcriptHighlights=trans_highlights,
  #   labels = "both",spaceHeight =0,strandSplit = FALSE,
  #  fontsize = 5,spaceWidth = 0,boxHeight = unit(2,"mm"),colorbyStrand = "#37a7db"
  # )


hl_start <-  min(cluster_posInfo$cluster_start) |>  as.numeric()
hl_end <-  max(cluster_posInfo$cluster_end) |> as.numeric()
  anno_highlight_plot <- annoHighlight(
    plot = plot_selectedGene,
    chrom = gene_info$gene_chr, 
    chromstart = hl_start,
    chromend = hl_end,
    y = 4.7, height = 0.3, just = c("left", "top"),
    default.units = "inches"
  )
  
  
  plotGG(plot =  genotype_plot, x = 0.2, y = 5.25, height = 2.35, width = 2.35)
  
  plotText(label =  pbs_results$phe_id, 
           x = 0.6,
           y=5.25,
           fontsize = 8, fontfamily = "Helvetica",
           just="left",
           fontcolor = "black")
  
  plotText(label =  pbs_results$rsID, 
           x = 1.25,
           y=7.6,
           fontsize = 8, fontfamily = "Helvetica",
           just="left",
           fontcolor = "black")
  
  #plotGG(plot =  gg_table, x = 1.5, y = 7.5, height = 1.5, width = 1)
  
  # plotText(label =  "PBS", 
  #          x = 0.65,
  #          y=4.7,
  #          fontsize = 8, fontfamily = "Helvetica",
  #          just="left",
  #          fontcolor = "#1e87a5")
  # 
  # plotText(label =  "FN-f", 
  #          x = 1.5,
  #          y=4.7,
  #          fontsize = 8, fontfamily = "Helvetica",
  #          just="left",
  #          fontcolor = "#FFB81C")
  # 
  # plotText(label =  paste0("\u03B2 = ", round(pbs_results$PBS_beta, 3)), 
  #          x = 0.65,
  #          y=4.85,
  #          fontsize = 7, fontfamily = "Helvetica",
  #          just="left",
  #          fontcolor = "black")
  # 
  # plotText(label =   paste0("\u03B2 = ", round(pbs_results$FNF_beta, 3)), 
  #          x = 1.5,
  #          y=4.85,
  #          fontsize = 7, fontfamily = "Helvetica",
  #          just="left",
  #          fontcolor = "black")
  # 
  # 
  # plotText(label =  paste0("pvalue = ", round(pbs_results$PBS_p, 3)), 
  #          x = 0.65,
  #          y=4.95,
  #          fontsize = 7, fontfamily = "Helvetica",
  #          just="left",
  #          fontcolor = "black")
  # 
  # plotText(label =   paste0("pvalue = ", round(pbs_results$FNF_p, 3)), 
  #          x = 1.5,
  #          y=4.95,
  #          fontsize = 7, fontfamily = "Helvetica",
  #          just="left",
  #          fontcolor = "black")
  # 
  # plotText(label =   paste0("pvalue = ", round(pbs_results$FNF_p, 3)), 
  #          x = 1.5,
  #          y=4.95,
  #          fontsize = 7, fontfamily = "Helvetica",
  #          just="left",
  #          fontcolor = "black")
  
  
  
 
}

dev.off()



