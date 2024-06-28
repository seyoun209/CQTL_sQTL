# Create sQTL relevent plots

setwd("/work/users/s/e/seyoun/CQTL_sQTL/output")
source("scripts/sQTL_rscripts/utils.R")

#-------------------------------------------------------------------------------
# Venn diagram - 1 
#-------------------------------------------------------------------------------

#This is sGene for the spling QTL including the condtional analysis

saveRDS(final_pbs_sig_qtl_cond, file="./01.qtltools_re/pbs_cond_sig_beta_maf_added.rds")
saveRDS(final_fnf_sig_qtl_cond, file="./01.qtltools_re/fnf_cond_sig_beta_maf_added.rds")
pbs_sig_qtl_cond_annot <- readRDS("./01.qtltools_re/pbs_cond_sig_beta_maf_added.rds")
fnf_sig_qtl_cond_annot <- readRDS("./01.qtltools_re/fnf_cond_sig_beta_maf_added.rds")

# adding Symboles

library(rtracklayer)
gtf_path <- "/work/users/s/e/seyoun/Ref/genome/gencode.v45.annotation.gtf"
gtf_data <- import(gtf_path)
gene_info <- gtf_data[gtf_data$type == "gene"]

genes_df <- data.frame(
  gene_id = mcols(gtf_data)$gene_id,
  gene_name = mcols(gtf_data)$gene_name,
  stringsAsFactors = FALSE
)

genes_df_unique <- genes_df %>%
  dplyr::rename(ensg = gene_id) %>%
  dplyr::distinct(ensg, gene_name, .keep_all = TRUE)  


sig_sGene_pbs_symbol <- pbs_sig_qtl_cond_annot %>%
  left_join(genes_df_unique, by="ensg")

sig_sGene_fnf_symbol <- fnf_sig_qtl_cond_annot %>%
  left_join(genes_df_unique, by="ensg")


pbs_sGene_sig <- sig_sGene_pbs_symbol %>%
  dplyr::filter(!is.na(gene_name)) %>%
  dplyr::select(gene_name) %>%
  unique() %>%
  unlist()

fnf_sGene_sig <- sig_sGene_fnf_symbol %>%
  dplyr::filter(!is.na(gene_name)) %>%
  dplyr::select(gene_name) %>%
  unique() %>%
  unlist()


#sGene -1 
sgene_list <- list("sGene_PBS"=pbs_sGene_sig, "sGene_FN-f"=fnf_sGene_sig)

sig_sGenes_all_tibble <- tibble(values = unique(c(pbs_sGene_sig, fnf_sGene_sig))) %>%
  mutate(PBS = values %in% pbs_sGene_sig,
         FNF = values %in% fnf_sGene_sig)


sGene_venn <- ggplot(sig_sGenes_all_tibble, aes(A = PBS, B = FNF)) +
  geom_venn(set_names = c("PBS", "FN-f"), 
            fill_color = c("PBS"="#BFDDFF","FN-f"="#FFDDA2"),
            set_name_color= c("PBS" = "#2057A7","FN-f" = "#F2BC40"),
            stroke_color = NA, auto_scale = TRUE, show_percentage = FALSE,
            text_size = 3, set_name_size = 4) +
  coord_fixed()  +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background =element_blank(),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        plot.margin = margin(t = 5, r = 10, b = 5, l = 5))
ggsave(filename = "output/results_plots/sqtl_plots/venn_sGenes.pdf", 
       plot =sGene_venn,width = 4, height = 4, units = "in")
save(sGene_venn, file = "output/results_plots/sqtl_plots/venn_sGenes.rda")

#Intron junctions -2

pbs_sGene_introns <- sig_sGene_pbs_symbol %>%
  dplyr::filter(!is.na(genomicLoc)) %>%
  dplyr::select(genomicLoc) %>%
  unique() %>%
  unlist()

fnf_sGene_introns <- sig_sGene_fnf_symbol %>%
  dplyr::filter(!is.na(genomicLoc)) %>%
  dplyr::select(genomicLoc) %>%
  unique() %>%
  unlist()

sIntrons_list <- list("Introns_PBS"=pbs_sGene_introns, "Introns_FN-f"=fnf_sGene_introns)

sig_introns_all_tibble <- tibble(values = unique(c(pbs_sGene_introns, fnf_sGene_introns))) %>%
  mutate(PBS = values %in% pbs_sGene_introns,
         FNF = values %in% fnf_sGene_introns)

sIntrons_venn <- ggplot(sig_introns_all_tibble, aes(A = PBS, B = FNF)) +
  geom_venn(set_names = c("PBS", "FN-f"), 
            fill_color = c("PBS"="#BFDDFF","FN-f"="#FFDDA2"),
            set_name_color= c("PBS" = "#2057A7","FN-f" = "#F2BC40"),
            stroke_color = NA, auto_scale = TRUE, show_percentage = FALSE,
            text_size = 3, set_name_size = 4) +
  coord_fixed()  +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background =element_blank(),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        plot.margin = margin(t = 5, r = 10, b = 5, l = 5))
ggsave(filename = "output/results_plots/sqtl_plots/venn_introns.pdf", 
       plot =sIntrons_venn,width = 4, height = 4, units = "in")
save(sIntrons_venn, file = "output/results_plots/sqtl_plots/venn_introns.rda")




#sSNPs_intronjunction_pairs -3


pbs_sGene_snpIntron_pairs <- sig_sGene_pbs_symbol %>%
  dplyr::mutate(snp_intron_pair = paste(var_id, genomicLoc, sep = "_")) %>%
  dplyr::filter(!is.na(snp_intron_pair)) %>%
  dplyr::select(snp_intron_pair) %>%
  unique() %>%
  unlist()

fnf_sGene_snpIntron_pairs <- sig_sGene_fnf_symbol %>%
  dplyr::mutate(snp_intron_pair = paste(var_id, genomicLoc, sep = "_")) %>%
  dplyr::filter(!is.na(snp_intron_pair)) %>%
  dplyr::select(snp_intron_pair) %>%
  unique() %>%
  unlist()


pairs_list <- list("pairs_PBS"=pbs_sGene_snpIntron_pairs, "paris_FN-f"=fnf_sGene_snpIntron_pairs)

sig_pairs_all_tibble <- tibble(values = unique(c(pbs_sGene_snpIntron_pairs, fnf_sGene_snpIntron_pairs))) %>%
  mutate(PBS = values %in% pbs_sGene_snpIntron_pairs,
         FNF = values %in% fnf_sGene_snpIntron_pairs)

snp_introns_pairs_venn <- ggplot(sig_pairs_all_tibble, aes(A = PBS, B = FNF)) +
  geom_venn(set_names = c("PBS", "FN-f"), 
            fill_color = c("PBS"="#BFDDFF","FN-f"="#FFDDA2"),
            set_name_color= c("PBS" = "#2057A7","FN-f" = "#F2BC40"),
            stroke_color = NA, auto_scale = TRUE, show_percentage = FALSE,
            text_size = 3, set_name_size = 4) +
  coord_fixed()  +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background =element_blank(),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        plot.margin = margin(t = 5, r = 10, b = 5, l = 5))
#p <-venn_font(sIntrons_venn, font = "Helvetica")
ggsave(filename = "output/results_plots/sqtl_plots/venn_snp_introns_pairs.pdf", 
       plot =snp_introns_pairs_venn,width = 4, height = 4, units = "in")
save(snp_introns_pairs_venn, file = "output/results_plots/sqtl_plots/enn_snp_introns_pairs.rda")


#Supp figure 4e venn diagrams
pdf(file = "./results_plots/sqtl_plots/fig4_supp_ven.pdf", width = 3, height = 5.1)
pageCreate(width = 3, height =5.1 , showGuides = FALSE)
plotText("D", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
load("results_plots/sqtl_plots/enn_snp_introns_pairs.rda")
load("results_plots/sqtl_plots/venn_introns.rda")
load("results_plots/sqtl_plots/venn_sGenes.rda")


plotText(label =  "Intron Junctions", 
         x = 1.75,
         y= 0.4,
         fontsize = 10, fontfamily = "Helvetica",
         just="center",
         fontcolor = "black")
plotGG(plot =  sIntrons_venn, x = 0.5, y = 0, height = 2.5, width = 2.5)

plotText(label =  "sQTL-sGene pairs", 
         x =1.75,
         y= 1.9,
         fontsize = 10, fontfamily = "Helvetica",
         just="center",
         fontcolor = "black")
plotGG(plot =  sGene_venn, x = 0.5, y = 1.6, height = 2.5, width = 2.5)

plotText(label = "sSNP-intron junction pairs", 
         x = 1.75,
         y=3.45,
         fontsize = 10, fontfamily = "Helvetica",
         just="center",
         fontcolor = "black")
plotGG(plot =  snp_introns_pairs_venn, x = 0.25, y = 3.1, height = 2.5, width = 2.5)

dev.off()




#-------------------------------------------------------------------------------
#Create venndiagram for the response QTL Figure 3B

# Function to count the rank
rank_counts <- function(counts, dataset_name) {
  df <- counts |> table() |> as_tibble()  |> dplyr::rename(rank="counts") |> dplyr::rename( counts ="n") %>%
    dplyr::mutate(dataset_name=dataset_name) %>% dplyr::rename(group="dataset_name")
  return(df)
}

# Prepare data for both PBS and FNF
all_cond_count_rank.df <- rbind(rank_counts(pbs_sig_qtl_cond_annot$rank, "PBS"), rank_counts(fnf_sig_qtl_cond_annot$rank, "FN-f"))
all_cond_count_rank <- all_cond_count_rank.df %>%
  mutate(rank = paste0(rank, "°"))
all_cond_count_rank <- all_cond_count_rank %>%
  mutate(rank = factor(rank, levels = c("3°", "2°", "1°", "0°"), ordered = TRUE)) %>%
  mutate(group = factor(group, levels = c("PBS","FN-f"), ordered = TRUE))
#dfm_only_cond <- dfm %>% dplyr::filter(rank != "0°")


library(ggpubr)
# Correctly assigning colors to the levels
palette_colr <- c("FN-f" = "#F2BC40", "PBS" = "#2057A7")

# Now call ggdotchartd

introns_dotchart <- ggdotchart(all_cond_count_rank, x = "rank", y = "counts",
           color = "group",                          # Column to define color groups
           palette = palette_colr,                   # Define colors for each group
           add = "segments",                         # Add segments from zero to the dot's value
           position = position_dodge(0.3),           # Dodge positions for clarity
           rotate = TRUE,                            # Horizontal plot
           group = "group",                          # Group data by 'dataset' for plotting
           dot.size = 3,                             # Size of dots
           ggtheme = theme_pubr()) +  
  labs(title = "Number of intron junctions") +       # Use ggpubr theme for aesthetics
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_markdown(size = 6, family = "Helvetica", margin = margin(r = -15)),
        text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.text.x = element_text(color = "black", size = 6),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        legend.title = element_blank(),
        legend.position = c(0.9, 0.15),
        plot.title = element_text(hjust = 0.5, family = "Helvetica", size = 10, margin = margin(b = -3))) +
  scale_x_discrete(limits = c("3°", "2°", "1°", "0°")) 
  
save(introns_dotchart, file = "results_plots/sqtl_plots/introns_dotchart.rda")
ggsave(filename ="results_plots/sqtl_plots/introns_dotchart.pdf",
       plot = introns_dotchart, width = 5, height = 5, units = "in")



#------------------------------------------------------------------------------
#boxplots for log2 Fold change of counts
library(plyr)

all_cond_count_rank_log <- all_cond_count_rank %>% add_row(rank = "3°", counts = NA, group = "FN-f") %>%
  mutate(log2 = log2(counts)) %>% 
  mutate(group = factor(group, levels = c("PBS","FN-f"), ordered = TRUE))

counts_rank_barplot <- ggplot(all_cond_count_rank_log, aes(x = rank, y = sapply(log2, FUN=function(x) ifelse(x==0, 0.1,x) ), fill = group)) +
  geom_bar(stat = "identity", position =position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("PBS"="#BFDDFF","FN-f"="#FFDDA2")) +
  scale_y_continuous( name = "log~2~ (counts) of sIntron junctions",
                      limits = c(0, 13), breaks = seq(0, 13, 1)) +
  coord_cartesian(clip = "off") +
  #geom_text(aes(label = sprintf("%.2f", round(log2, 2))), position = position_dodge(width = 0.9), vjust = -0.5, size = 3)+
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
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
        strip.text.x.top = element_text(size = 8, margin = margin(b = 5)),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.9, 0.9),
        panel.spacing.y = unit(0.5, "cm")) +
  scale_x_discrete(limits = c("0°", "1°", "2°", "3°")) 

save(counts_rank_barplot, file = "output/results_plots/sqtl_plots/counts_rank_vennplot.rda")
ggsave(filename ="output/results_plots/sqtl_plots/counts_rank_barplot.pdf",
       plot = counts_rank_barplot, width = 5, height = 5, units = "in")

#-------------------------------------------------------------------------------
# Venndiagram for the response QTL
#This will be the lmer version
response_pbs_results <- readRDS("01.qtltools_re/conditional_pbs/response_pbs_re_lmer.rds")
response_fnf_results <- readRDS("01.qtltools_re/conditional_fnf/response_fnf_re_lmer.rds")

#add symbols
response_pbs_re_wSymbol <- response_pbs_results %>%
  left_join(genes_df_unique, by="ensg")
response_fnf_re_wSymbol <- response_fnf_results %>%
  left_join(genes_df_unique, by="ensg")

response_pbs_results_nolmer <- readRDS("01.qtltools_re/conditional_pbs/response_pbs_re_no_lmer.rds")
response_fnf_results_nolmer <- readRDS("01.qtltools_re/conditional_fnf/response_fnf_re_no_lmer.rds")

#add symbols
response_pbs_re_wSymbol_nolmer <- response_pbs_results_nolmer %>%
  left_join(genes_df_unique, by="ensg")
response_fnf_re_wSymbol_nolmer <- response_fnf_results_nolmer %>%
  left_join(genes_df_unique, by="ensg")

#geno_data to count each 5 donors with each variant genotype
pbsGeno_raw <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/pbs_geno/06.subset_sigSNps/recodeA_pbs.raw")
colnames(pbsGeno_raw) <- sub("_.*", "", colnames(pbsGeno_raw))

fnfGeno_raw <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/fnf_geno/06.subset_sigSNps/recodeA_fnf.raw")
colnames(fnfGeno_raw) <- sub("_.*", "", colnames(fnfGeno_raw))

pbs_transpose_df <- pbsGeno_raw %>%
  column_to_rownames(var = "IID")  %>%
  dplyr::select(-FID, -MAT,-PAT, -SEX, -PHENOTYPE) %>% 
  as.data.frame() %>%  # Ensure it's still a data frame
  t() |> as.data.frame() |>
  rownames_to_column(var = "var_id")

fnf_transpose_df <- fnfGeno_raw %>%
  column_to_rownames(var = "IID")  %>%
  dplyr::select(-FID, -MAT,-PAT, -SEX, -PHENOTYPE) %>% 
  as.data.frame() %>%  # Ensure it's still a data frame
  t() |> as.data.frame()|>
  rownames_to_column(var = "var_id")

pbs_allele_counts <- pbs_transpose_df %>%
  rowwise() %>%
  dplyr::mutate(count_0 = sum(across(everything()) == 0),
         count_1 = sum(across(everything()) == 1),
         count_2 = sum(across(everything()) == 2)) %>%
  ungroup() %>%
  dplyr::select(var_id,count_0, count_1, count_2)

pbs_allele_counts_subset <- pbs_allele_counts %>%
  dplyr::filter(count_0 >= 5 & count_1 >= 5 & count_2 >= 5)


fnf_allele_counts <- fnf_transpose_df %>%
  rowwise() %>%
  dplyr::mutate(count_0 = sum(across(everything()) == 0),
                count_1 = sum(across(everything()) == 1),
                count_2 = sum(across(everything()) == 2)) %>%
  ungroup() %>%
  dplyr::select(var_id,count_0, count_1, count_2)

fnf_allele_counts_subset <- fnf_allele_counts %>%
  dplyr::filter(count_0 >= 5 & count_1 >= 5 & count_2 >= 5)

# 
# #lmer
# highconf_response_pbs <- response_pbs_re_wSymbol |> dplyr::filter(interaction_pval <= 0.05) |>
#   dplyr::filter(abs(delta_beta) >= 0.2) |>
#   dplyr::filter(var_id  %in% fnf_allele_counts_subset$var_id) # you can use either one. Since it's genotype and technically same.
# 
# highconf_response_fnf <- response_fnf_re_wSymbol |> dplyr::filter(interaction_pval <= 0.05) |>
#   dplyr::filter(abs(delta_beta) >= 0.2) |>
#   dplyr::filter(var_id  %in% fnf_allele_counts_subset$var_id)# you can use either one. Since it's genotype and technically same.
# 
# sig_response_pbs <-  response_pbs_re_wSymbol |> dplyr::filter(interaction_pval <= 0.05)
# sig_response_fnf <-  response_fnf_re_wSymbol |> dplyr::filter(interaction_pval <= 0.05)



highconf_response_pbs_nolmer <- response_pbs_re_wSymbol_nolmer |> dplyr::filter(interaction_pval <= 0.05) |>
  dplyr::filter(abs(delta_beta) >= 0.2) |>
  dplyr::filter(var_id  %in% fnf_allele_counts_subset$var_id) # you can use either one. Since it's genotype and technically same.

highconf_response_fnf_nolmer <- response_fnf_re_wSymbol_nolmer |> dplyr::filter(interaction_pval <= 0.05) |>
  dplyr::filter(abs(delta_beta) >= 0.2) |>
  dplyr::filter(var_id  %in% fnf_allele_counts_subset$var_id)# you can use either one. Since it's genotype and technically same.

sig_response_pbs_nolmer <-  response_pbs_re_wSymbol_nolmer |> dplyr::filter(interaction_pval <= 0.05)
sig_response_fnf_nolmer <-  response_fnf_re_wSymbol_nolmer |> dplyr::filter(interaction_pval <= 0.05)

write.csv(sig_response_pbs_nolmer, "output/01.qtltools_re/response_qtl/reseponsesQTL_PBS_significant.csv")
write.csv(sig_response_fnf_nolmer, "output/01.qtltools_re/response_qtl/reseponsesQTL_FNF_significant.csv")


#venndiagram response sQTL
#ensg_re_sig_pbs <- highconf_response_pbs_nolmer$ensg[!is.na(highconf_response_pbs_nolmer$ensg)]
#ensg_re_sig_fnf <- highconf_response_fnf_nolmer$ensg[!is.na(highconf_response_fnf_nolmer$ensg)]
ensg_re_sig_pbs <- sig_response_pbs_nolmer$ensg[!is.na(sig_response_pbs_nolmer$ensg)]
ensg_re_sig_fnf <- sig_response_fnf_nolmer$ensg[!is.na(sig_response_fnf_nolmer$ensg)]

sgene_list <- list("sGene_PBS"=ensg_re_sig_pbs, "sGene_FN-f"=ensg_re_sig_fnf)

sig_sGenes_list_response <- tibble(values = unique(c(ensg_re_sig_pbs, ensg_re_sig_fnf))) %>%
  dplyr::mutate("PBS" = values %in% ensg_re_sig_pbs,
                "FNF" = values %in% ensg_re_sig_fnf)

sig_response_sGenesVenn <- ggplot(sig_sGenes_list_response, aes(A = PBS, B = FNF)) +
  geom_venn( fill_color =c("PBS_sGenes" = "#A5C5FF","FN-f_sGenes" = "#FCCE52"),
               stroke_color = NA,
               show_percentage = FALSE,
               text_size = 3,
             set_names= c(PBS="PBS_sGenes",FNF="FN-f_sGenes"),
             set_name_color= c("PBS_sGenes" = "#2057A7","FN-f_sGenes" = "#F2BC40"),
               set_name_size = 4) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "transparent", color = "transparent"))
#sGene_venn <- venn_font(sGene_venn, font = "Helvetica")
ggsave(filename = "output/results_plots/sqtl_plots/response_venndiagram.pdf", 
       plot =sig_response_sGenesVenn,width = 4, height = 4, units = "in")
save(sig_response_sGenesVenn, file = "output/results_plots/sqtl_plots/response_venndiagram.rda")



#high confidence response QTL pie chart

ensg_re_highconfid_pbs <- highconf_response_pbs_nolmer$ensg[!is.na(highconf_response_pbs_nolmer$ensg)] |> unique()
ensg_re_highconfid_fnf <- highconf_response_fnf_nolmer$ensg[!is.na(highconf_response_fnf_nolmer$ensg)] |> unique()

ensg_re_highconfid_pbs_df <-data.frame(ensg = ensg_re_highconfid_pbs) %>%
  mutate(group = ifelse(ensg %in% ensg_re_highconfid_fnf,  "shared","high-confidence"))

ensg_re_highconfid_fnf_df <-data.frame(ensg = ensg_re_highconfid_fnf) %>%
  mutate(group = ifelse(ensg %in% ensg_re_highconfid_pbs,  "shared","high-confidence"))

fnf_group_chart <- table(ensg_re_highconfid_fnf_df$group) |> as.data.frame()
colnames(fnf_group_chart)[1] <- "group"

fnf_pieColors <- c("high-confidence" = lighten("#F2BC40",0.2),
                   "shared" = lighten("#FFDDA2",0.4))

fnf_label_text <- paste("n=", fnf_group_chart$Freq[1])
high_conf_FNF_Pie <- create_pie_chart(fnf_group_chart, fnf_pieColors, 'FNF')
save(high_conf_FNF_Pie, file = "output/results_plots/sqtl_plots/high_conf_FNF_Pie.rda")
pbs_group_chart <- table(ensg_re_highconfid_fnf_df$group) |> as.data.frame()
colnames(pbs_group_chart)[1] <- "group"


pbs_pieColors <- c("high-confidence" = darken("#BFDDFF",0.4),
               "shared" = lighten("#BFDDFF",0.4))

pbs_label_text <- paste("n=", fnf_group_chart$Freq[1])
high_conf_PBS_Pie <- create_pie_chart(pbs_group_chart, pbs_pieColors,"PBS")

save(high_conf_FNF_Pie, file = "output/results_plots/sqtl_plots/high_conf_PBS_Pie.rda")



#-------------------------------------------------------------------------------
#Density plot for sQTL and rs-QTL

# sQTl distance calc first -from splice Donor 
library(scales)

pbs_sig_qtl_distance <- pbs_sig_qtl_cond_annot %>%
  mutate(distance_from_sd = case_when(
    phe_strd == "-" ~ phe_to - var_from,
    phe_strd == "+" ~ phe_from - var_from
  )) %>%
  mutate(abs_distance = abs(distance_from_sd),
         group="PBS")

fnf_sig_qtl_distance <- fnf_sig_qtl_cond_annot %>%
  mutate(distance_from_sd = case_when(
    phe_strd == "-" ~ phe_to - var_from,
    phe_strd == "+" ~ phe_from - var_from
  )) %>%
  mutate(abs_distance = abs(distance_from_sd),
         group="FNF")

all_distance_info_sig_qtl <- rbind(pbs_sig_qtl_distance,fnf_sig_qtl_distance)
all_distance_info_sig_qtl$group <- factor(all_distance_info_sig_qtl$group ,levels = c("PBS","FNF"))

binwidth <- 5000


summarize_sig_qtl_dis <- all_distance_info_sig_qtl %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(
    mean = mean(abs_distance), 
    median = median(abs_distance))

pbs_y <- data.frame(group ="PBS",
y_at_median = get_y_for_median(pbs_sig_qtl_distance,summarize_sig_qtl_dis$median[1],binwidth),
y_at_mean =get_y_for_median(pbs_sig_qtl_distance,summarize_sig_qtl_dis$mean[1],binwidth))
fnf_y <- data.frame(group ="FNF",
                    y_at_median = get_y_for_median(fnf_sig_qtl_distance,summarize_sig_qtl_dis$median[2],binwidth),
                    y_at_mean =get_y_for_median(fnf_sig_qtl_distance,summarize_sig_qtl_dis$mean[2],binwidth))
y_all <- rbind(pbs_y,fnf_y)
summarize_dis_qtl_fi <- left_join(summarize_sig_qtl_dis,y_all,by="group")

#line_plot

distQTL_line_plot <-ggplot(all_distance_info_sig_qtl, aes(x = abs_distance, group=group,colour=group)) +
  geom_line(stat = "bin", binwidth = binwidth) +
  ggtitle("Distance to sQTLs to splicing donor site")+
  scale_x_continuous(limits = c(0, 240000), name = "Distance (in kb) to splicing donor site (SD)", expand = c(0,0),
                     labels = label_number(scale = .001, big.mark = "")) +
  scale_y_continuous(name = "Number of sQTLs", expand = c(0,0),
                     limits = c(0,1200)) +
  scale_color_manual(values = c("PBS" = "#0067B9", "FNF" = "#FCCE52"),labels = c("PBS" = "PBS", "FNF" = "FN-f")) +
  geom_point(data = summarize_dis_qtl_fi, aes(x = mean[1]-500, y = y_at_mean[1]), color = "#2057A7", size = 2) +
  geom_text(data = summarize_dis_qtl_fi, aes(x =mean[1]+35000, y =y_at_mean[1]), label = sprintf("PBS-Mean = %.1f kb", summarize_dis_qtl_fi$mean[1]/1000),
            color = darken("#2057A7",0.2),size = 3) +
  geom_point(data = summarize_dis_qtl_fi, aes(x = mean[2]+750, y = y_at_mean[2]), color = "#F2BC40", size = 2) +
  geom_text(data = summarize_dis_qtl_fi, aes(x =mean[2]+40000, y =y_at_mean[2]), label = sprintf("FN-f-Mean = %.1f kb", summarize_dis_qtl_fi$mean[2]/1000),
            color = darken("#F2BC40",0.2),size = 3) +
  theme(text = element_text(family = "Helvetica"),
        panel.background = element_rect(fill = 'transparent', color = "transparent"),
        plot.background = element_rect(fill = 'transparent', color = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.9,0.8),
        legend.background = element_blank(),
        plot.title = element_text(hjust = 0.5, family = "Helvetica",
                                  size = 10, margin = margin(b=-3)),
        axis.title = element_text(size = 9),
        axis.line = element_line(color = "black", 
                                 linewidth = 0.25),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.ticks.y = element_line(color = "black",
                                    linewidth = 0.25),
        axis.text = element_text(color = "black", size = 8),
        plot.margin = margin(t = 5, r = 10, b = 5, l = 5))

save(distQTL_line_plot, file = "output/results_plots/sqtl_plots/distQTL_line_plot.rda")

# rs-qtl plot

pbs_response_qtl_distance <- sig_response_pbs_nolmer %>%
  mutate(distance_from_sd = case_when(
    phe_strd == "-" ~ phe_to - var_from,
    phe_strd == "+" ~ phe_from - var_from
  )) %>%
  mutate(abs_distance = abs(distance_from_sd),
         group="PBS")

fnf_response_qtl_distance <- sig_response_fnf_nolmer %>%
  mutate(distance_from_sd = case_when(
    phe_strd == "-" ~ phe_to - var_from,
    phe_strd == "+" ~ phe_from - var_from
  )) %>%
  mutate(abs_distance = abs(distance_from_sd),
         group="FNF")

all_distance_info_rs_qtl <- rbind(pbs_response_qtl_distance,fnf_response_qtl_distance)
all_distance_info_rs_qtl$group <- factor(all_distance_info_rs_qtl$group ,levels = c("PBS","FNF"))

binwidth <- 5000

summarize_rs_qtl_dis <- all_distance_info_rs_qtl %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(
    mean = mean(abs_distance), 
    median = median(abs_distance))

pbs_rs_y <- data.frame(group ="PBS",
                    y_at_median = get_y_for_median(pbs_response_qtl_distance,summarize_rs_qtl_dis$median[1],binwidth),
                    y_at_mean =get_y_for_median(pbs_response_qtl_distance,summarize_rs_qtl_dis$mean[1],binwidth))
fnf_rs_y <- data.frame(group ="FNF",
                    y_at_median = get_y_for_median(fnf_response_qtl_distance,summarize_rs_qtl_dis$median[2],binwidth),
                    y_at_mean =get_y_for_median(fnf_response_qtl_distance,summarize_rs_qtl_dis$mean[2],binwidth))
y_rs_all <- rbind(pbs_rs_y,fnf_rs_y)
summarize_rs_qtl_fi <- left_join(summarize_rs_qtl_dis,y_rs_all,by="group")

#rs-qtl distance line plot
dist_rsQTL_line_plot <- ggplot(all_distance_info_rs_qtl, aes(x = abs_distance, group=group,colour=group)) +
  geom_line(stat = "bin", binwidth = binwidth) +
  ggtitle("Distance to response-sQTLs to splicing donor site")+
  scale_x_continuous(limits = c(0, 240000), name = "Distance (in kb) to splicing donor site (SD)", expand = c(0,0),
                     labels = label_number(scale = .001, big.mark = "")) +
  scale_y_continuous(name = "Number of response sQTLs", expand = c(0,0),
                     limits = c(0,450)) +
  scale_color_manual(values = c("PBS" = "#0067B9", "FNF" = "#FCCE52"),labels = c("PBS" = "PBS", "FNF" = "FN-f")) +
  geom_point(data = summarize_rs_qtl_fi, aes(x = mean[1]+1500, y = y_at_mean[1]), color = "#2057A7", size = 2) +
  geom_text(data = summarize_rs_qtl_fi, aes(x =mean[1]+45000, y =y_at_mean[1]), label = sprintf("PBS-Mean = %.1f kb", summarize_rs_qtl_fi$mean[1]/1000),
            color = darken("#2057A7",0.2),size = 3) +
  geom_point(data = summarize_rs_qtl_fi, aes(x = mean[2]-1000, y = y_at_mean[2]), color = "#F2BC40", size = 2) +
  geom_text(data = summarize_rs_qtl_fi, aes(x =mean[2]+40000, y =y_at_mean[2]), label = sprintf("FN-f-Mean = %.1f kb", summarize_rs_qtl_fi$mean[2]/1000),
            color = darken("#F2BC40",0.2),size = 3) +
  theme(text = element_text(family = "Helvetica"),
        panel.background = element_rect(fill = 'transparent', color = "transparent"),
        plot.background = element_rect(fill = 'transparent', color = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5, family = "Helvetica",
                                  size = 10, margin = margin(b=-3)),
        legend.position = c(0.9,0.8),
        legend.background = element_blank(),
        axis.title = element_text(size = 9),
        axis.line = element_line(color = "black", 
                                 linewidth = 0.25),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.ticks.y = element_line(color = "black",
                                    linewidth = 0.25),
        axis.text = element_text(color = "black", size = 8),
        plot.margin = margin(t = 5, r = 10, b = 5, l = 5))

save(dist_rsQTL_line_plot, file = "output/results_plots/sqtl_plots/dist_rsQTL_line_plot.rda")



