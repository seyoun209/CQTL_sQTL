# Making OA donors compared plot with the current CTL vs FNF 
## Author: Seyoun Byun
## Date: 03.28.2024
## Edited:06.15.2024 (Tried to see categorical helps with the confounder)

#-------------------------------------------------------------------------------
setwd("/work/users/s/e/seyoun/CQTL_sQTL")
library(limma)
library(magrittr)
library(data.table)
library(dplyr)
library(leafviz)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(RColorBrewer)
library(plotgardener)
library(grid)
library(sva)
library(gridtext)
library(ggtext)
source("scripts/utils/utils.R")

pbs_oa <- load("output/clu_oa/PBSvsOA.Rdata")
for (name in pbs_oa) {
  # Construct the new name by prefixing with "fnf_"
  new_name <- paste0(name,"_oa")
  
  # Assign the object to the new name in the global environment
  assign(new_name, get(name))
  remove(name)
}


ratios_oa <- counts_oa %>%
  mutate(clu = str_split_fixed(rownames(counts_oa), ":", 4)[,4]) %>%
  group_by(clu) %>%
  mutate_all(normalize_column) %>%
  ungroup() %>%
  as.data.frame() %>%
  set_rownames(rownames(counts_oa)) %>%
  dplyr::select(-clu)

ratios_oa_qc = ratios_oa[rowMeans(is.na(ratios_oa)) <= 0.4,,drop=F ] #Try to remove 40% or less NA values in each rows
#From 134376 to 133944 --> it is dropped 432
row_means = rowMeans(ratios_oa_qc, na.rm = T) #calculate the mean without NAs in the rows. 
row_means_outer = outer(row_means, rep(1,ncol(ratios_oa_qc)))# making outlier to the NAs
ratios_oa_qc[is.na(ratios_oa_qc)] = row_means_outer[is.na(ratios_oa_qc)] # instead of Na, add the rowmean
colnames(ratios_oa_qc) <- colnames(ratios_oa_qc)
ratios_oa_qc <- cbind(rownames(ratios_oa_qc), ratios_oa_qc)
colnames(ratios_oa_qc)[1] <- c('Junction')
#write.table(ratios_oa_qc, file = "output/clu_oa/ratio_oa.txt",sep='\t',quote=F,row.names=F,col.names=T)
ratios_oa_qc <- read.table("output/clu_oa/ratio_oa.txt",sep='\t',header=T) |> as.data.frame()


#-------------------------------------------------------------------------------
#Limma for the differential analysis - FNF
#-------------------------------------------------------------------------------
config <- yaml::read_yaml("config/rna_prcoess.yaml")
conditions_to_include <- unlist(strsplit("CTL FNF", " "))

load("output/combined_meta_data.RData") # It is loading name is combined_data

selected_columns <- c("Donor","ID","Condition","Sex", "Age","FragmentBatch","RIN","RNAextractionKitBatch","RNAshippedDate") #select column needed it based
meta_data <- combined_data %>% dplyr::select(all_of(selected_columns))
#Ancestry 
ancestry_df <- fread("/proj/phanstiel_lab/Data/processed/CQTL/geno/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_COA6_COA7/CQTL_COA_01_GDA8_COA2_01_COA3_01_GDA8_COA4_COA5_COA6_COA7_predictedAncestry.csv")
ancestry_df$Donor <- sub("^.*_(AM[0-9]+)_.*$", "\\1", ancestry_df$Donor)

ancestry_OA_df <- fread("/proj/phanstiel_lab/Data/processed/CQTL/geno/COA8_OA/ancestry/CQTL_COA8_predictedAncestry.csv")
ancestry_OA_df_filtered <- ancestry_OA_df %>%
  filter(grepl("^OA\\d+", Donor))
ancestry_OA_df_filtered$Donor <- gsub("_r2","",ancestry_OA_df_filtered$Donor)

ancestry_cqtl <-rbind(ancestry_df,ancestry_OA_df_filtered)
meta_cqtl <- merge(meta_data,ancestry_cqtl,by="Donor",all.x=TRUE)
#write.table(meta_cqtl, file = "output/clu_fnf/meta_cqtl",sep='\t',quote=F,row.names=F,col.names=T)

meta_cqtl <- fread("output/clu_fnf/meta_cqtl")
meta_ctl_fnf <- meta_cqtl[!meta_cqtl$Donor %in% config$samples_to_omit, ] %>%
  dplyr::filter(Condition %in% conditions_to_include)
meta_ctl_oa <- meta_ctl_fnf %>% dplyr::filter(Condition %in% c('CTL','OA')) %>%
  mutate(FragmentBatch = ifelse(Condition == "OA", "batch0", FragmentBatch))
  
#-------------------------------------------------------------------------------
#corrected batch with Limma
#-------------------------------------------------------------------------------
# PBS vs OA won't fix the batch correction because only PBS data have Fragmentbatch information and RNAextractionKitbatch or SequencingBatch have more tha 1 in PBS data. 
# rownames(ratios_oa_qc) <- ratios_oa_qc[,1]
# ratios_oa_qc <- ratios_oa_qc[, -1]
# 
# model_oa <-model.matrix(~Condition,data=meta_ctl_oa)
# 
# limma_oa_batchcorr <- limma::removeBatchEffect(ratios_oa_qc,batch=as.factor(meta_ctl_oa$Donor))
# limma_oa_batchcorr_junc <- cbind(rownames(limma_oa_batchcorr), limma_oa_batchcorr)
# colnames(limma_oa_batchcorr_junc)[1] <- c('Junction')
# write.table(limma_oa_batchcorr_junc, file = "output/clu_oa/psi_oa_limma_batchCorr",sep='\t',quote=F,row.names=F,col.names=T)
# limma_oa_batchcorr_df <- limma_oa_batchcorr |> as.data.frame()
# limma_oa_batchcorr_df <- calculate_delta_psi(limma_oa_batchcorr_df, "CTL", "OA")

OA_deltapsiCalc_df <- calculate_delta_psi(ratios_oa_qc, "CTL", "OA")

#-------------------------------------------------------------------------------
#Draw PCA plot
remove_vars <- c("Row.names","Donor","RIN","ID", "Condition","Sex","Age","Predicted_Ancestry","FragmentBatch","RNAextractionKitBatch","RNAshippedDate")
oa_pca_prep <-pca_prep(ratios_oa_qc,meta_ctl_oa,remove_vars)

before_batch_oa <- ggplot(oa_pca_prep$pca_data, aes(x = PC1, y = PC2, color = Condition ), label = oa_pca_prep$Donor ) +
  geom_point() +
  #geom_text_repel(size = 2, box.padding = unit(0.1, "lines")) +
  scale_color_manual(values = c("CTL" = "#74CCEC", "OA" = "#FAB394")) +
  labs(
    title = "Before batch correction",
    x = paste0("PC1 Variance:", round(oa_pca_prep$variance_explained[1] * 100, 2), "%"),
    y = paste0("PC2 Variance:", round(oa_pca_prep$variance_explained[2] * 100, 2), "%")) +
  theme(plot.background = element_rect(fill='transparent', color=NA))+
  theme_bw(base_size = 10) +
  theme(aspect.ratio = 1,panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"))

ggsave(filename = "output/results_plots/Supplementary_figures/supp_fig1_pca_oa.pdf",
       plot = before_batch_oa, width = 4.5, height = 4.5, units = "in")
save(before_batch_oa, file = "output/results_plots/Supplementary_figures/supp_figure1_pca_oa.rda")


#-------------------------------------------------------------------------------
# Finding the background genes for OA

split_names <- strsplit(rownames(counts_oa), ":")
split_df <- matrix(unlist(split_names), ncol = 4, byrow = TRUE) |> as.data.frame()
split_df <- split_df %>%
  dplyr::mutate(start = as.numeric(V2),
                end = as.numeric(V3)) %>%
  dplyr::select(-V2 ,-V3)
colnames(split_df) <- c("chr",  "clusterID","start", "end")
combined_df <- inner_join(split_df, introns_gencode_renamed, by = c("chr", "start", "end"))
combined_df_nodot <- combined_df %>%
  mutate(V5 = sub("\\..*$", "", V5))

subset_background_genes_ensg <- combined_df_nodot $V5 |> unique() 
#write.table(subset_background_genes_ensg, file = "output/clu_oa/background_OA_genes_set.txt",sep='\t',quote=F,row.names=F,col.names=F)

#-------------------------------------------------------------------------------
introns_oa_pval_include <- join_introns_deltapsi_fdr(OA_deltapsiCalc_df,introns_oa,"./output/clu_oa/ctlvsoa_ds_cluster_significance.txt")
introns_oa_sig <- introns_oa_pval_include %>% dplyr::filter(p.adjust <= 0.05)
oa_maxCluster <- introns_oa_sig %>%
  mutate(abs_deltapsi = abs(deltapsi_batch)) %>%  # Add a new column for the absolute value of deltapsi
  group_by(clusterID) %>%
  # Use slice_max to select the row with the maximum absolute deltapsi value
  dplyr::slice(which.max(abs_deltapsi)) %>%
  # Ensure that ensemblID is not "."
  dplyr::filter(ensemblID != ".") %>%
  # Optionally, remove the abs_deltapsi column if it's no longer needed
  dplyr::select(-abs_deltapsi)

sig_psi_maxCluster_oa_psi2 <- oa_maxCluster[abs(oa_maxCluster$deltapsi_batch) >= 0.20,]
#Gene for the cluster Max
sig_psi_maxCluster_oa_psi2_modified <- sig_psi_maxCluster_oa_psi2 %>%
  mutate(ensemblID = sub("\\..*$", "", ensemblID))
sig_gene_20_percent_diff_oa <- sig_psi_maxCluster_oa_psi2_modified$ensemblID |> unique()  #Save it to ENSG

sig_psi_maxCluster_oa_psi15 <- oa_maxCluster[abs(oa_maxCluster$deltapsi_batch) >= 0.15,]
#Gene for the cluster Max
sig_psi_maxCluster_oa_psi15_modified <- sig_psi_maxCluster_oa_psi15 %>%
  mutate(ensemblID = sub("\\..*$", "", ensemblID))
sig_gene_15_percent_diff_oa <- sig_psi_maxCluster_oa_psi15_modified$ensemblID |> unique()  #Save it to ENSG
write.table(sig_gene_20_percent_diff_oa, file = "output/clu_oa/sig_gene_20_percent_diff.txt",sep='\t',quote=F,row.names=F,col.names=F)
write.table(sig_gene_15_percent_diff_oa, file = "output/clu_oa/sig_gene_15_percent_diff.txt",sep='\t',quote=F,row.names=F,col.names=F)

#-------------------------------------------------------------------------------
#KEGG GO
#-------------------------------------------------------------------------------

## all significant differential splicing genes.
#system("scripts/Differential_splicing/run_homer.sh /work/users/s/e/seyoun/CQTL_sQTL/output/clu_oa/sig_gene_20_percent_diff.txt /work/users/s/e/seyoun/CQTL_sQTL/output/clu_oa/background_OA_genes_set.txt /work/users/s/e/seyoun/CQTL_sQTL/output/clu_oa/homer/homer_sig_diffsplicing_all_fdr05_psi2")
#system("scripts/Differential_splicing/run_homer.sh /work/users/s/e/seyoun/CQTL_sQTL/output/clu_oa/sig_gene_15_percent_diff.txt /work/users/s/e/seyoun/CQTL_sQTL/output/clu_oa/background_OA_genes_set.txt /work/users/s/e/seyoun/CQTL_sQTL/output/clu_oa/homer/homer_sig_diffsplicing_all_fdr05_psi15")

#GO

go_data <- read_delim("output/clu_oa/homer/homer_sig_diffsplicing_all_fdr05_psi15/biological_process.txt") |>
  mutate(pval = exp(1)^logP) |>
  dplyr::filter(pval < 0.01)
sig_go <- reduceGO(go_data,
                   category = "GO")

go_table <- sig_go |>
  dplyr::select(-`Entrez Gene IDs`, -pval, -logP, -size, -termUniqueness, -termUniquenessWithinCluster,
                -termDispensability, -category) |>
  relocate(`-log10pval`, .after = Enrichment) |>
  arrange(desc(`-log10pval`))

write_csv(go_table, file = "output/clu_oa/table/GO_sig.csv")


Siggo_plotting <- sig_go |>
  dplyr::filter(Term == parentTerm) |>
  dplyr::filter(parentTerm %in%  unique(go_table$parentTerm)[1:10]) |>
  arrange(`-log10pval`)

Siggo_plotting$parentTerm <- factor(Siggo_plotting$parentTerm, levels = Siggo_plotting$parentTerm)



OAGO_barplots <- ggplot(Siggo_plotting, aes(x = `-log10pval`, y = parentTerm, fill = category)) +
  geom_vline(xintercept = 2, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 4, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 5, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 8, color = "grey75", alpha = 0.4) +
  #geom_vline(xintercept =15 , color = "grey75", alpha = 0.4) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 8), expand = c(0, 0), name = "-log~10~pval",
                     breaks = seq(0, 8, 2)) +
  scale_fill_manual(values = "#C8F0BF") +
  facet_wrap(~category, ncol = 1, strip.position = "left", scales = "free_y") +
  geom_text(aes(x = 0, label = parentTerm), hjust = 0, family = "Helvetica",
            size = 2.5) +
  theme(panel.background = element_rect(fill = 'transparent', color = "transparent"),
        plot.background = element_rect(fill = 'transparent', color = "transparent"),
        text = element_text(family = "Helvetica"),
        legend.position = "None",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown(size = 6),
        axis.text.x = element_text(color = "black", size = 4),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_line(linewidth = 0.25),
        strip.text = element_blank(),
        panel.spacing = unit(0, "mm"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8)) +
  ggtitle("GO Terms")

ggsave(filename = "output/results_plots/Figure1_differntial_splicing/OAGO_barplots.pdf",
       plot = OAGO_barplots, width = 5, height = 5, units = "in")
save(OAGO_barplots, file = "output/results_plots/Figure1_differntial_splicing/OAGO_barplots.rda")

##pathway Reacome + kegg
# Read in from Homer
reactome_data <- read_delim("output/clu_oa/homer/homer_sig_diffsplicing_all_fdr05_psi15/reactome.txt") |>
  mutate(pval = exp(1)^logP) |>
  dplyr::filter(pval < 0.01) |>
  distinct(Term, .keep_all = TRUE) |>
  mutate(`-log10pval` = -log10(pval)) |>
  mutate(category = "Reactome Pathway")

kegg_data <- read_delim("output/clu_oa/homer/homer_sig_diffsplicing_all_fdr05_psi15/kegg.txt") |>
  mutate(pval = exp(1)^logP) |>
  dplyr::filter(pval < 0.01) |>
  distinct(Term, .keep_all = TRUE) |>
  mutate(`-log10pval` = -log10(pval)) |>
  mutate(category = "KEGG Pathway")

combined_pathway <- bind_rows(reactome_data, kegg_data)

## Format and write to table
pathway_table <- combined_pathway |>
  dplyr::select(-logP, -pval, -`Entrez Gene IDs`) |>
  relocate(`-log10pval`, .after = Enrichment) |>
  arrange(desc(`-log10pval`))
write_csv(pathway_table, file = "output/clu_oa/table/pathway_table_sig.csv")

# Plot top 10 significant for each category
pathway_plotting <- pathway_table |>
  dplyr::filter(Term %in% unique(pathway_table$Term)
  ) |>
  arrange(`-log10pval`)

pathway_plotting$Term <- factor(pathway_plotting$Term, levels = pathway_plotting$Term)

OApathway_barplots <- ggplot(pathway_plotting, aes(x = `-log10pval`, y = Term)) +
  geom_vline(xintercept = 1, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 2, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 3, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 4, color = "grey75", alpha = 0.4) +
  # geom_vline(xintercept = 5, color = "grey75", alpha = 0.4) +
  geom_bar(stat = "identity", fill="#C8F0BF") +
  scale_x_continuous(expand = c(0, 0), name = "-log~10~pval", limits = c(0, 4),
                     breaks = seq(0, 4, 1)) +
  #facet_wrap(~category, ncol = 1, strip.position = "left", scales = "free_y") +
  geom_text(aes(x = 0, label = Term), hjust = 0, family = "Helvetica",
            size = 2.5) +
  theme(panel.background = element_rect(fill = 'transparent', color = "transparent"),
        plot.background = element_rect(fill = 'transparent', color = "transparent"),
        text = element_text(family = "Helvetica"),
        legend.position = "None",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_markdown(size = 6),
        axis.text.x = element_text(color = "black", size = 4),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_line(linewidth = 0.25),
        strip.text = element_blank(),
        panel.spacing = unit(0, "mm"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8)) +
  ggtitle("Pathways")
ggsave(filename = "output/results_plots/Figure1_differntial_splicing/OApathway_barplot.pdf",
       plot = OApathway_barplots, width = 5, height = 5, units = "in")
save(OApathway_barplots, file = "output/results_plots/Figure1_differntial_splicing/OApathway_barplots.rda")
#-------------------------------------------------------------------------------
#boxplot for Figure1_E

introns_oa_subset_sig <- introns_oa_sig[abs(introns_oa_sig$deltapsi_batch) >= 0.15,]
introns_oa_subset_sig$loc <- paste0(introns_oa_subset_sig$chr,":",introns_oa_subset_sig$start,":",introns_oa_subset_sig$end)

introns_oa_subset_sig_up <- introns_oa_subset_sig %>% dplyr::filter(deltapsi > 0) %>% 
  dplyr::select(loc,deltapsi_batch,gene)
introns_oa_subset_down <- introns_oa_subset_sig %>% dplyr::filter(deltapsi < 0) %>% 
  dplyr::select(loc,deltapsi_batch,gene)

fnf_all_OA_subset <- bind_rows(introns_oa_subset_sig_up |> mutate(group = "Up in OA"),
                               introns_oa_subset_down |> mutate(group = "Down in OA")) |>
  mutate(group = factor(group, levels = c("Up in OA", "Down in OA")))

test <-c("NR4A1","PIK3CD","ABLIM1","STS","ANKRD36C","ABLIM1","SNRNP70","NCOR2","PRKX","CAMKK2")


fnf_all_OA_subset <- fnf_all_OA_subset |>
  mutate(highlight = case_when(gene %in% c("STS", "CAMKK2", "ANKRD36C", "NCOR2","PRKX") ~ "up",
                               gene %in% c("SNRNP70", "NR4A1", "PIK3CD", "ABLIM1") ~ "down")) |>
  mutate(highlight = factor(highlight, levels = c("up", "down")))






fnf_psi_batchCorr_junction <- cbind(rownames(limma_psi_batchremove.df), limma_psi_batchremove.df)
colnames(fnf_psi_batchCorr_junction)[1] <- c('Junction')
ratios_fnf <- fnf_psi_batchCorr_junction %>%
  separate(Junction, into = c("part1", "part2", "part3","part4"), sep = ":", fill = "right", extra = "drop") %>%
  mutate(loc = paste(part1, part2, part3, sep = ":")) %>%
  dplyr::select(-part1, -part2, -part3, -part4)

sig_OA_subset_psi15 <- fnf_all_OA_subset %>%
  inner_join(ratios_fnf %>% dplyr::select(loc, deltaPSI), by = "loc")

OA_from_all_fnf_up <- ratios_fnf %>% 
  dplyr::filter(loc %in% introns_oa_subset_sig_up$loc) %>%
  dplyr::select(loc,deltaPSI)
OA_from_all_fnf_down <- ratios_fnf %>% 
  dplyr::filter(loc %in% introns_oa_subset_down$loc) %>%
  dplyr::select(loc,deltaPSI)


introns_ids_common <- c(rownames(OA_from_all_fnf_up), rownames(OA_from_all_fnf_down))

down_fnf_common_subset <- introns_fnf_sig %>% 
  dplyr::filter(phe_id %in% rownames(OA_from_all_fnf_down))  %>% 
  arrange(deltapsi_batch) %>%
  dplyr::filter(deltapsi_batch < 0)
up_fnf_common_subset <- introns_fnf_sig %>% 
  dplyr::filter(phe_id %in% rownames(OA_from_all_fnf_up))  %>% 
  arrange(deltapsi_batch) %>%
  dplyr::filter(deltapsi_batch > 0)

intersect(unique(introns_fnf_all_sig$gene), unique(c(down_fnf_common_subset$genes,up_fnf_common_subset$genes)))

# wilcox test

up_test_OA <- wilcox.test(x = as.numeric(OA_from_all_fnf_up$deltaPSI),
                         mu=0,
                          alternative = "greater")

down_test_OA <- wilcox.test(x = as.numeric(OA_from_all_fnf_down$deltaPSI),
                            mu=0,
                            alternative = "less")

cat("p-value:", format.pval(up_test_OA$p.value, digits = 3), "\n")
cat("p-value:", format.pval(down_test_OA$p.value, digits = 3), "\n")



#oa_boxplot_plot(data = fnf_all_OA_subset, up_wilcox_test = up_test_OA,
#                down_wilcox_test = down_test_OA,
#                x = unit(5.5, "native"),
#                y = unit(5.6, "native"), width = 3.5, height = 3)

oa_boxplots <- ggplot(sig_OA_subset_psi15, aes(x = group, y = deltaPSI, fill = group)) +
  geom_hline(yintercept = 0, lty = 2, color = "grey25", linewidth = 0.25) +
  geom_jitter(width = 0.2, color = "grey40", size = 0.25) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.5, alpha = 0.4, width = 0.5, color = c("#FFB81C", '#005587')) +
  stat_boxplot(geom = "errorbar", width = 0.5, color = c("#FFB81C", '#005587')) +
  scale_color_manual(values = c(darken("#FFB81C", 0.3), darken('#005587', 0.3))) +
  scale_fill_manual(values = c("#FFB81C", '#005587')) +
  scale_y_continuous(name = "delatPSI in response to FN-f",
                     limits = c(-0.6, 0.6), breaks = seq(-0.6, 0.6, 0.3)) +
  coord_cartesian(clip = "off") +
  geom_text(aes(label = "*"), y = 0.55, vjust = -1, size = 4, fontface = "bold") +
  theme(strip.placement = "outside",
        axis.line.y = element_line(linewidth = 0.35),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_markdown(size = 8, family = "Helvetica",
                                        margin = margin(r = -15)),
        text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.text.x = element_text(size = 8, margin = margin(b = -1),colour=c(darken("#FFB81C", 0.3), darken('#005587', 0.3))),
        strip.background = element_blank(),
        strip.text.x.bottom = element_markdown(size = 8, margin = margin(t = 1)),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, family = "Helvetica",
                                  size = 10, margin = margin(b = -3)))

ggsave(filename = "output/results_plots/Figure1_differntial_splicing/Fig1E_OA_boxplot.pdf",
       plot = oa_boxplots, width = 6, height = 4.5, units = "in")
save(oa_boxplots, file = "output/results_plots/Supplementary_figures/Fig1E_OA_boxplot.rda")





