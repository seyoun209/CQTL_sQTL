# Making boxplots for the QTL for shared, PBS or FNF treated
## Author: Seyoun Byun
## Date: 03.08.2024
## Edited:06.13.2024
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
library(ggrepel)
library(leafcutter)
library(gridExtra)
library(grid)
library(rrvgo)
library(gprofiler2)
source("scripts/utils/utils.R")


#leafviz("output/clu_fnf/PBSvsFNF.Rdata")
pbs_fnf <- load("output/clu_fnf/PBSvsFNF.Rdata")
introns_gencode <- fread("tools/leafcutter/gencode_hg38_all_introns.bed.gz")
introns_gencode_renamed <- introns_gencode %>% dplyr::rename(chr = "V1", start = "V2", end = "V3")
#make normalized counts for the PSI

for (name in pbs_fnf) {
  # Construct the new name by prefixing with "fnf_"
  new_name <- paste0(name,"_fnf")
  
  # Assign the object to the new name in the global environment
  assign(new_name, get(name))
  remove(name)
}

ratios_fnf <- counts_fnf %>%
  mutate(clu = str_split_fixed(rownames(counts_fnf), ":", 4)[,4]) %>%
  group_by(clu) %>%
  mutate_all(normalize_column) %>%
  ungroup() %>%
  as.data.frame() %>%
  set_rownames(rownames(counts_fnf)) %>%
  dplyr::select(-clu)

ratios_fnf_qc = ratios_fnf[rowMeans(is.na(ratios_fnf)) <= 0.4,,drop=F ] #Try to remove 40% or less NA values in each rows
#147127 to 144656 --> 2471 rows dropped
row_means_fnf = rowMeans(ratios_fnf_qc, na.rm = T) #calculate the mean without NAs in the rows. 
row_means_outer_fnf = outer(row_means_fnf, rep(1,ncol(ratios_fnf_qc)))# making outlier to the NAs
ratios_fnf_qc[is.na(ratios_fnf_qc)] = row_means_outer_fnf[is.na(ratios_fnf_qc)] # instead of Na, add the rowmean
colnames(ratios_fnf_qc) <- colnames(counts_fnf)
ratios_fnf_qc <- cbind(rownames(ratios_fnf_qc), ratios_fnf_qc)
colnames(ratios_fnf_qc)[1] <- c('Junction')
#write.table(ratios_fnf_qc, file = "output/clu_fnf/ratio_fnf.txt",sep='\t',quote=F,row.names=F,col.names=T)



#-------------------------------------------------------------------------------
#Sample adding ancestry and covariates - FNF
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



#-------------------------------------------------------------------------------
#Limma for the differential analysis - FNF
#-------------------------------------------------------------------------------

# Omit samples and filter rows based on configuration and conditions
meta_ctl_fnf <- meta_cqtl[!meta_cqtl$Donor %in% config$samples_to_omit, ] %>%
  dplyr::filter(Condition %in% conditions_to_include)
#write.table(meta_ctl_fnf, file = "output/clu_fnf/meta_ctl_fnf",sep='\t',quote=F,row.names=F,col.names=T)
meta_ctl_fnf <- fread("output/clu_fnf/meta_ctl_fnf")
psi_fnf <- read.table("output/clu_fnf/ratio_fnf.txt",sep="\t",header=T)
rownames(psi_fnf) <- psi_fnf[,1]
psi_fnf <- psi_fnf[, -1]

model_cov_fnf <-model.matrix(~Condition+as.factor(Donor)+
                               as.factor(RNAshippedDate)+as.factor(RNAextractionKitBatch)+
                             as.factor(FragmentBatch),
                             data=meta_ctl_fnf)
limma_psi_batchremove <- limma::removeBatchEffect(psi_fnf,
                                                  batch= as.factor(meta_ctl_fnf$RNAshippedDate),
                                                  batch1=as.factor(meta_ctl_fnf$RNAextractionKitBatch),
                                                  batch2=as.factor(meta_ctl_fnf$FragmentBatch),
                                                  batch3=as.factor(meta_ctl_fnf$Donor))
limma_psi_batchremove_junction <- cbind(rownames(limma_psi_batchremove), limma_psi_batchremove)
colnames(limma_psi_batchremove_junction)[1] <- c('Junction')
#write.table(limma_psi_batchremove_junction, file = "output/clu_fnf/psi_fnf_limma_batch_corrected",sep='\t',quote=F,row.names=F,col.names=T)

limma_psi_batchremove_junction <- fread("output/clu_fnf/psi_fnf_limma_batch_corrected")
rownames(limma_psi_batchremove_junction) <- c(limma_psi_batchremove_junction$Junction)
limma_psi_batchremove.df <- limma_psi_batchremove_junction %>% dplyr::select(-"Junction")



#------------------------------------------------------------------------------
#corrected PCA
remove_vars <- c("Row.names","Donor","RIN","ID", "Condition","Sex","Age","Predicted_Ancestry","FragmentBatch","RNAextractionKitBatch","RNAshippedDate")
fnf_pca_batchCorr <-pca_prep(limma_psi_batchremove.df,meta_ctl_fnf,remove_vars)
ratios_fnf_qc_num <- ratios_fnf_qc[,-1]
fnf_pca_prep <-pca_prep(ratios_fnf_qc_num,meta_ctl_fnf,remove_vars)

after_batch_fnf <- ggplot(fnf_pca_batchCorr$pca_data, aes(x = PC1, y = PC2, color = Condition )) + #label = Donor, )) +
  geom_point() +
  #geom_text_repel(size = 2, box.padding = unit(0.1, "lines")) +
  scale_color_manual(values = c("CTL" = "#74CCEC", "FNF" = "#FAB394")) +
  labs(
    title = "After batch correction",
    x = paste0("PC1 Variance:", round(fnf_pca_batchCorr$variance_explained[1] * 100, 2), "%"),
    y = paste0("PC2 Variance:", round(fnf_pca_batchCorr$variance_explained[2] * 100, 2), "%")) +
  theme(plot.background = element_rect(fill='transparent', color=NA))+
  theme_bw(base_size = 10) +
  theme(aspect.ratio = 1,panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"))
before_batch_fnf <- ggplot(fnf_pca_prep$pca_data, aes(x = PC1, y = PC2, color = Condition )) + #label = Donor, )) +
  geom_point() +
  #geom_text_repel(size = 2, box.padding = unit(0.1, "lines")) +
  scale_color_manual(values = c("CTL" = "#74CCEC", "FNF" = "#FAB394")) +
  labs(
    title = "After batch correction",
    x = paste0("PC1 Variance:", round(fnf_pca_prep$variance_explained[1] * 100, 2), "%"),
    y = paste0("PC2 Variance:", round(fnf_pca_prep$variance_explained[2] * 100, 2), "%")) +
  theme(plot.background = element_rect(fill='transparent', color=NA))+
  theme_bw(base_size = 10) +
  theme(aspect.ratio = 1,panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"))
supp_figure1_pca <- grid.arrange(before_batch_fnf, after_batch_fnf, ncol = 2, widths = c(2, 2))
ggsave(filename = "output/results_plots/Supplementary_figures/supp_fig1_pca.pdf",
       plot = supp_figure1_pca, width = 8, height = 4.5, units = "in")
save(supp_figure1_pca, file = "output/results_plots/Supplementary_figures/supp_figure1_pca.rda")


#-------------------------------------------------------------------------------

# Finding the background genes
split_names <- strsplit(rownames(counts_fnf), ":")
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
#write.table(subset_background_genes_ensg, file = "output/clu_fnf/background_gene_set.txt",sep='\t',quote=F,row.names=F,col.names=F)

#-------------------------------------------------------------------------------
#Using the new deltaPSI to combine into introns_fnf

limma_psi_batchremove.df <- calculate_delta_psi(limma_psi_batchremove.df, "CTL", "FNF")
introns_fnf_pval_include <- join_introns_deltapsi_fdr(limma_psi_batchremove.df,introns_fnf,"./output/clu_fnf/ctlvsfnf_ds_cluster_significance.txt")
#save(introns_fnf_pval_include, file="./output/clu_fnf/introns_fnf_joinAll")
#-------------------------------------------------------------------------------
#finding the significant
introns_fnf_sig <- introns_fnf_pval_include %>% dplyr::filter(p.adjust <= 0.05)
fnf_maxCluster <- introns_fnf_sig %>%
  mutate(abs_deltapsi = abs(deltapsi_batch)) %>%  # Add a new column for the absolute value of deltapsi
  group_by(clusterID) %>%
  # Use slice_max to select the row with the maximum absolute deltapsi value
  dplyr::slice(which.max(abs_deltapsi)) %>%
  # Ensure that ensemblID is not "."
  dplyr::filter(ensemblID != ".") %>%
  # Optionally, remove the abs_deltapsi column if it's no longer needed
  dplyr::select(-abs_deltapsi)

sig_psi_maxCluster <- fnf_maxCluster[abs(fnf_maxCluster$deltapsi_batch) >= 0.15,]
#Gene for the cluster Max
sig_psi_maxCluster_nodot <- sig_psi_maxCluster %>%
  mutate(ensemblID = sub("\\..*$", "", ensemblID))
sig_gene_20_percent_diff <- sig_psi_maxCluster_nodot$ensemblID |> unique()  #Save it to ENSG

sig_psi_maxCluste_psi15 <- fnf_maxCluster[abs(fnf_maxCluster$deltapsi_batch) >= 0.15,]
#Gene for the cluster Max
sig_psi_maxCluster_nodot_psi15 <- sig_psi_maxCluste_psi15 %>%
  mutate(ensemblID = sub("\\..*$", "", ensemblID))
sig_gene_15_percent_diff <- sig_psi_maxCluster_nodot_psi15$ensemblID |> unique()  #Save it to ENSG
write.table(sig_gene_20_percent_diff, file = "output/clu_fnf/sig_gene_20_percent_diff.txt",sep='\t',quote=F,row.names=F,col.names=F)
write.table(sig_gene_15_percent_diff, file = "output/clu_fnf/sig_gene_15_percent_diff.txt",sep='\t',quote=F,row.names=F,col.names=F)


#-------------------------------------------------------------------------------
#This is finding all the intron junctions that are FDR < 0.05 and deltaPSI > 20%

introns_fnf_all_sig <- introns_fnf_sig %>% dplyr::filter(abs(deltapsi_batch) >= 0.15)

#read the meta data for the Age, sex, Ancestry
load("output/combined_meta_data.RData") # It is loading name is combined_data

selected_columns <- c("Donor","ID","Condition","Sex", "Age") #select column needed it based
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
filtered_meta_samples <- meta_cqtl %>%
  dplyr::filter(!grepl("OA", ID))
meta_samples <-filtered_meta_samples[,2:6] |> as.data.frame()

#-------------------------------------------------------------------------------
#KEGG GO
#-------------------------------------------------------------------------------

## all significant differential splicing genes.
#system("scripts/Differential_splicing/run_homer.sh /work/users/s/e/seyoun/CQTL_sQTL/output/clu_fnf/sig_gene_20_percent_diff.txt /work/users/s/e/seyoun/CQTL_sQTL/output/clu_fnf/background_gene_set.txt /work/users/s/e/seyoun/CQTL_sQTL/output/clu_fnf/homer/homer_sig_diffsplicing_all_fdr05_psi2")
#system("scripts/Differential_splicing/run_homer.sh /work/users/s/e/seyoun/CQTL_sQTL/output/clu_fnf/sig_gene_15_percent_diff.txt /work/users/s/e/seyoun/CQTL_sQTL/output/clu_fnf/background_gene_set.txt /work/users/s/e/seyoun/CQTL_sQTL/output/clu_fnf/homer/homer_sig_diffsplicing_all_fdr05_psi15")

#GO

go_data <- read_delim("output/clu_fnf/homer/homer_sig_diffsplicing_all_fdr05_psi15/biological_process.txt") |>
  mutate(pval = exp(1)^logP) |>
  dplyr::filter(pval < 0.01)
sig_go <- reduceGO(go_data,
                   category = "GO")

go_table <- sig_go |>
  dplyr::select(-`Entrez Gene IDs`, -pval, -logP, -size, -termUniqueness, -termUniquenessWithinCluster,
                -termDispensability, -category) |>
  relocate(`-log10pval`, .after = Enrichment) |>
  arrange(desc(`-log10pval`))

write_csv(go_table, file = "output/clu_fnf/table/GO_sig.csv")


Siggo_plotting <- sig_go |>
  dplyr::filter(Term == parentTerm) |>
  dplyr::filter(parentTerm %in%  unique(go_table$parentTerm)[1:10]) |>
  arrange(`-log10pval`)

Siggo_plotting$parentTerm <- factor(Siggo_plotting$parentTerm, levels = Siggo_plotting$parentTerm)



GO_barplots <- ggplot(Siggo_plotting, aes(x = `-log10pval`, y = parentTerm, fill = category)) +
  geom_vline(xintercept = 2, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 4, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 5, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 8, color = "grey75", alpha = 0.4) +
  #geom_vline(xintercept =15 , color = "grey75", alpha = 0.4) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 8), expand = c(0, 0), name = "-log~10~pval",
                     breaks = seq(0, 8, 2)) +
  scale_fill_manual(values = "#c1daf3") +
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

ggsave(filename = "output/results_plots/Figure1_differntial_splicing/GO_barplots.pdf",
       plot = GO_barplots, width = 5, height = 5, units = "in")
save(GO_barplots, file = "output/results_plots/Figure1_differntial_splicing/GO_barplots.rda")

##pathway Reacome + kegg
# Read in from Homer
reactome_data <- read_delim("output/clu_fnf/homer/homer_sig_diffsplicing_all_fdr05_psi15/reactome.txt") |>
  mutate(pval = exp(1)^logP) |>
  dplyr::filter(pval < 0.01) |>
  distinct(Term, .keep_all = TRUE) |>
  mutate(`-log10pval` = -log10(pval)) |>
  mutate(category = "Reactome Pathway")

kegg_data <- read_delim("output/clu_fnf/homer/homer_sig_diffsplicing_all_fdr05_psi15/kegg.txt") |>
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
write_csv(pathway_table, file = "output/clu_fnf/table/pathway_table_sig.csv")

# Plot top 10 significant for each category
pathway_plotting <- pathway_table |>
  dplyr::filter(Term %in% c("Signal Transduction", 
                            "Signaling by FGFR1 in disease",
                            "DCC mediated attractive signaling",
                            "Fluid shear stress and atherosclerosis",
                            "PTK6 Regulates RHO GTPases, RAS GTPase and MAP kinases",
                            "Interaction With The Zona Pellucida",
                            "Insulin signaling pathway",
                            "IkBA variant leads to EDA-ID",
                            "Aldosterone-regulated sodium reabsorption",
                            "Osteoclast differentiation"
                            )) |>
  arrange(`-log10pval`)

pathway_plotting$Term <- factor(pathway_plotting$Term, levels = pathway_plotting$Term)

pathway_barplots <- ggplot(pathway_plotting, aes(x = `-log10pval`, y = Term)) +
  geom_vline(xintercept = 1, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 2, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 3, color = "grey75", alpha = 0.4) +
  geom_vline(xintercept = 4, color = "grey75", alpha = 0.4) +
 # geom_vline(xintercept = 5, color = "grey75", alpha = 0.4) +
  geom_bar(stat = "identity", fill="#c1daf3") +
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
ggsave(filename = "output/results_plots/Figure1_differntial_splicing/pathway_barplot.pdf",
       plot = pathway_barplots, width = 5, height = 5, units = "in")
save(pathway_barplots, file = "output/results_plots/Figure1_differntial_splicing/pathway_barplots.rda")
#-------------------------------------------------------------------------------

#Heatmap
rownames(meta_ctl_fnf) <- meta_ctl_fnf$ID
psi_fnf_limma_corr_num <- limma_psi_batchremove.df[,-203]

psi_sig_deltapsi15 <- psi_fnf_limma_corr_num[rownames(psi_fnf_limma_corr_num) %in% sig_psi_maxCluste_psi15$phe_id,]


cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_subset_norm <- t(apply(psi_sig_deltapsi15, 1, cal_z_score))

#z-score
brks <- seq(min(data_subset_norm, na.rm = TRUE), max(data_subset_norm, na.rm = TRUE), length = 51)
brks <- seq(-2,2,length.out=50) 

# Define the age group intervals
age_breaks <- c(30, 40, 50, 60, 70, 80, Inf)

# Labels for the age groups
age_labels <- c("30-39", "40-49", "50-59", "60-69", "70-79", "80+")

# Create age groups using the cut function
age_groups <- cut(meta_ctl_fnf$Age, breaks = age_breaks, labels = age_labels, right = FALSE)
meta_ctl_fnf$Age_range <- age_groups

age_palette <- c("#ffcba4","#cca283","#997a62","#806652","#665142","#332921")

# Create a named vector of colors for the age groups
age_colors <- setNames(age_palette, age_labels)


my_colour = list(
  Condition  = c("CTL" = "#9FCCE4", "FNF" = "#FAB394"),
  Sex  = c( "M" =  "#0075B0", "F" = "#F8B7CD"),
  Predicted_Ancestry= c("AMR"="#F5BC9F", "EUR"="#FAF1D2" ,"AFR"="#86CBB7","SAS"="#6B7EA4"),
  Age_range = age_colors)

pheatmap::pheatmap(data_subset_norm,
                   annotation_col = meta_ctl_fnf,
                   annotation_colors = my_colour,
                   show_rownames = FALSE,
                   clustering_method = "complete",
                   show_colnames = FALSE,
                   treeheight_row = 0,
                   treeheight_col = 0,
                   fontsize = 5,
                   width = 1,
                   height = 1,
                   #cutree_rows = 2,
                   #cutree_cols = 2,
                   #labels_row = labels,
                   color = colorRampPalette(c("#097EA4", "black", "#BFA527"))(50),
                   breaks = brks
)


heatmapColors <- colorRampPalette(c("#097EA4", "black", "#BFA527"))(7)

ID <- colnames(data_subset_norm)
meta_ctl_fnf <- meta_ctl_fnf %>% arrange(factor(ID, levels=ID)) %>%
  arrange(Condition)

sample_info <- meta_ctl_fnf[, c("Age_range","Predicted_Ancestry","Sex", "Condition")]
colnames(sample_info) <- c("Age","Ancestry","Sex","Condition")

my_colour = list(
  Condition  = c("CTL" = "#B8B8B8", "FNF" = "#4A4A4A"),
  Sex  = c("F" = "#DD8492", "M" = "#4788BA"),
  Ancestry= c("AFR"="#C74A53", "AMR"="#EFB06E","EUR"= "#5BAD58", "SAS"="#177F97"),
  Age = age_colors)



colAnn <- HeatmapAnnotation(df = sample_info,
                            col = my_colour,
                            gap = unit(0.5, 'mm'),
                            annotation_name_gp= gpar(fontsize = 7),
                            simple_anno_size = unit(3, "mm"),
                            annotation_legend_param  = list(
                              Age = list(
                                title = "Age",
                                title_position = "leftcenter",
                                title_gp = gpar(fontsize = 7),
                                labels_gp = gpar(fontsize = 6),
                                grid_width = unit(2, "mm"),
                                grid_height = unit(1, "mm"),
                                at = c("30-39","40-49","50-59","60-69","70-79","80+"),
                                labels = c("30","40","50","60","70","80+"),
                                ncol = 1
                              ),
                              Condition  = list(
                                title = "Condition",
                                title_gp = gpar(fontsize = 7),
                                labels_gp = gpar(fontsize = 6),
                                title_position = "leftcenter",
                                grid_width = unit(2, "mm"),
                                grid_height = unit(1, "mm"),
                                at = c("CTL", "FNF"),
                                labels = c("PBS","FNF"),
                                ncol = 1
                              ),
                              Sex = list(
                                title = "Sex",
                                title_position = "leftcenter",
                                title_gp = gpar(fontsize = 7),
                                labels_gp = gpar(fontsize = 6),
                                grid_width = unit(2, "mm"),
                                grid_height = unit(1, "mm"),
                                at = c("M","F"),
                                labels = c("Male","Female"),
                                ncol = 1
                              ),
                              Ancestry= list(
                                title = "Ancestry",
                                title_position = "leftcenter",
                                title_gp = gpar(fontsize = 7),
                                labels_gp = gpar(fontsize = 6),
                                grid_width = unit(2, "mm"),
                                grid_height = unit(1, "mm"),
                                at = c("AFR","AMR","EUR","SAS"),
                                labels = c("AFR","AMR","EUR","SAS"),
                                ncol = 1
                                
                              )))


column_order <- meta_ctl_fnf$ID
data_subset_norm_ordered <- data_subset_norm[, column_order,drop = FALSE]

hmap <- Heatmap(
  data_subset_norm_ordered,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,  # Disable column clustering
  column_order = column_order,  # Specify column order
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  row_dend_reorder = FALSE,
  column_dend_reorder = FALSE,
  top_annotation = colAnn,
  col = colorRamp2(seq(-3, 3), heatmapColors)
)



heatmapLegend <- Legend(at = c(-3, 3),
                        col_fun = colorRamp2(breaks = seq(-3, 3),
                                             colors = heatmapColors),
                        border = NA,
                        title_gp = gpar(fontsize = 0),
                        labels_gp = gpar(fontfamily = "Helvetica",
                                         fontsize = 8),
                        legend_width = unit(4.325, "in"),
                        grid_height = unit(0.11, "in"),
                        direction = "horizontal")

heatmapGrob <- grid.grabExpr(draw(hmap,
                                  show_annotation_legend = FALSE,
                                  show_heatmap_legend = FALSE,
                                  background = "transparent"))
heatmapLegendGrob <- grid.grabExpr(draw(heatmapLegend))

save(heatmapGrob, file = "output/results_plots/Figure1_differntial_splicing/heatmapGrob.rda")
save(heatmapLegendGrob, file = "output/results_plots/Figure1_differntial_splicing/heatmapLegendGrob.rda")
#-------------------------------------------------------------------------------
# plotting figure1
#-------------------------------------------------------------------------------

pdf(file = "output/results_plots/Figure1_differntial_splicing/figure1.pdf",   # The directory you want to save the file in
    width = 11.5, # The width of the plot in inches
    height = 9.7)
pageCreate(width = 11.5, height = 9.7, showGuides = FALSE)
plotText("A", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
load(file="output/results_plots/Figure1_differntial_splicing/heatmapGrob.rda")
load(file="output/results_plots/Figure1_differntial_splicing/heatmapLegendGrob.rda")

plotGG(plot = heatmapGrob, x = 0.25, y = 0.2, height = 4.5, width = 5)

# Colorbar
plotGG(plot = heatmapLegendGrob, x = 0.35, y = 4.7,
       width = 4.325, height = 0.11)


# Colorbar title
plotText(label = "Relative Percent Spliced In (PSI)", fontfamily = "Helvetica",
         fontsize = 8, x = 2.45, y = 4.8, just = "top")

y=0.3
x=5.1
# Age Legend
for (i in 1:length(age_colors)) { 
  plotRect(x = unit(5, "in") + unit(3*i, "mm"),
           y = 0.3, width = unit(3, "mm"), 
           height = unit(1, "mm"), linecolor = NA, fill = age_colors[i],
           just = c("top" ,"left"))
}

ageText <- c("30","40","50","60","70","80+")

for(i in 1:length(ageText)){
  plotText(label = ageText[i], x = unit(5.05, "in") + unit(3.1*(i), "mm"),
           y =y+0.07 ,
           fontsize = 5, fontfamily = "Helvetica")
}


# Race Legend
Ancestry= c("AFR"="#C74A53", "AMR"="#EFB06E","EUR"= "#5BAD58", "SAS"="#177F97")

for (i in 1:4) {
  plotRect(x = unit(5.05, "in") + unit(4 * i, "mm"),
           y = y+0.07*2, width = unit(4, "mm"),
           height = unit(1, "mm"), linecolor = NA,
           fill = Ancestry[i], just = c("top", "left"))
}


raceText <- c("AFR","AMR","EUR","SAS")

for(i in 1:4){
  plotText(label = raceText[i], x = unit(5.13, "in") + unit(4*(i), "mm"),
           y = y+0.07*2+0.08,
           fontsize = 4, fontfamily = "Helvetica")
}

# Sex legend
plotRect(x = unit(5.1, "in") + unit(3*6, "mm"), 
         y = y+0.07*4, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = my_colour[['Sex']]['M'],
         just = "right")
plotRect(x = unit(5.1, "in") + unit(3*5, "mm"), 
         y = y+0.07*4, width = unit(3, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = my_colour[['Sex']]['F'],
         just = "right")
plotText(label = "M", x = unit(x, "in") + unit(3*5.5, "mm"),
         y = y+0.07*5,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "F", x = unit(x, "in") + unit(3*4.5, "mm"),
         y =y+0.07*5,fontsize = 5, fontfamily = "Helvetica")

# Condition legend
plotRect(x = unit(x, "in") + unit(3*6, "mm"), 
         y =  y+0.07*6, width = unit(4, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = "#4A4A4A",
         just = "right")
plotRect(x = unit(x, "in") + unit(3*4.7, "mm"), 
         y =y+0.07*6, width = unit(4, "mm"), 
         height = unit(1, "mm"), linecolor = NA, fill = "#B8B8B8",
         just = "right")
plotText(label = "FNF", x = unit(x, "in") + unit(3*5.4, "mm"),
         y =  y+0.07*7,
         fontsize = 5, fontfamily = "Helvetica")
plotText(label = "PBS", x = unit(x, "in") + unit(3*4, "mm"),
         y =y+0.07*7,
         fontsize = 5, fontfamily = "Helvetica")





#------------------------------------------------------------------------------
#Plot B
plotText("B", x = 6, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

for (i in 1:5) {
  if(i == 4){
    plotSegments(
      x0 = unit(6.5, "in"),
      y0 = 0.2+(0.15 * 4),
      x1 = unit(7.5, "in")+unit(7, "mm"),
      y1 = 0.2+(0.15 * 4),
      default.units = "inches",
      lwd = 1.5, lty = 1,"darkgrey",
      just ="left"
    )
    plotRect(x = unit(6.5, "in") +unit(7, "mm") ,
             y = 0.2+(0.15 * i), width = unit(3, "mm"),
             height = unit(2, "mm"), linecolor = NA, fill = "#2473b5",
             just = "left")
    plotRect(x = unit(7.5, "in") ,
             y = 0.2+(0.15 * i), width = unit(3, "mm"),
             height = unit(2, "mm"), linecolor = NA, fill = "#2473b5",
             just = "right")
  }
  
  else if(i == 5){
    plotRect(x = unit(6.5, "in") ,
             y =0.2+(0.15 * i), width = unit(7, "mm"),
             height = unit(2, "mm"), linecolor = NA, fill = "#434e57",
             just = "left")
    plotRect(x = unit(7.5, "in") ,
             y =0.2+(0.15 * i), width = unit(7, "mm"),
             height = unit(2, "mm"), linecolor = NA, fill = "#96bce3",
             just = "left")
    plotSegments(
      x0 = unit(6.5, "in") +unit(7, "mm"),
      y0 = 0.2+(0.15 * i),
      x1 = 7.5,
      y1 = 0.2+(0.15 * i),
      default.units = "inches",
      lwd = 1.5, lty = 1,"darkgrey",
      just ="left"
    )
    
  }
  
  else{
    plotRect(x = unit(6.5, "in") ,
             y = 0.2+(0.15 * i), width = unit(7, "mm"),
             height = unit(2, "mm"), linecolor = NA, fill = "#434e57",
             just = "left")
    plotRect(x = unit(7.5, "in") ,
             y = 0.2+(0.15 * i), width = unit(7, "mm"),
             height = unit(2, "mm"), linecolor = NA, fill = "#434e57",
             just = "left")
    plotSegments(
      x0 = unit(6.5, "in") +unit(7, "mm"),
      y0 =0.2+(0.15 * i),
      x1 = 7.5,
      y1 = 0.2+ (0.15 * i),
      default.units = "inches",
      lwd = 1.5, lty = 1,"darkgrey",
      just ="left"
    )
  }
}


sp_text <- c('annotated',"cryptic_5'","cryptic_3'","cryptic_unanchored","novel annotated pair")

for (i in 1:5) {
  plotText(label = sp_text[i], x = unit(7.85, "in"),
           y = 0.2+(0.15 * i),
           fontsize = 7, fontfamily = "Helvetica",
           just="left"
           #fontface = "bold"
  )
}


anno_text <- c(': annotated SS',": unannotated SS",": unpaired SS")
annocolor <- c("#434e57" ,"#2473b5","#96bce3")
for (i in 1:3) {
  plotRect(x = unit(5.75+(0.75*i),"in") ,
           y = 0.2+0.95, width = unit(1.5, "mm"),
           height = unit(1.5, "mm"), linecolor = NA, fill = annocolor[i],
           just = "left")
  
}

plotText(label = 'annotated SS', x = 6.58,
         y = 0.2+0.95,
         fontsize = 5, fontfamily = "Helvetica",
         just ="left", fontcolor = "#434e57")
plotText(label = "unannotated SS", x = 6.58+0.77,
         y =0.2+ 0.95,
         fontsize = 5, fontfamily = "Helvetica",
         just ="left", fontcolor = "#2473b5")
plotText(label = "unpaired SS", x = 6.58+0.77*2,
         y =0.2+ 0.95,
         fontsize = 5, fontfamily = "Helvetica",
         just ="left",fontcolor ="#96bce3")


dchart <- table(introns_fnf_all_sig$verdict) |> as.data.frame()
dchart$Percentage <- round((dchart$Freq / sum(dchart$Freq)),4)*100
dchart$ymax <- cumsum(dchart$Freq / sum(dchart$Freq))
dchart$ymin <- c(0, head(dchart$ymax, n=-1))

# Compute label position
dchart$labelPosition <- (dchart$ymax + dchart$ymin) / 2

# Compute a good label
lab.pos = cumsum(dchart$Percentage)-.5*(dchart$Percentage)
colnames(dchart)[1] <- "Verdict"
pieColors <- c("annotated"="lightgrey",
               "cryptic_fiveprime"="#96bce3" ,
               "cryptic_threeprime"="#3f85cc",
               "cryptic_unanchored"="#085099",
               "novel annotated pair"="#043363")
pieVerdict <-ggplot(data = dchart, aes(x = "", y = Percentage, fill = Verdict)) +
  geom_bar(stat = "identity") +
  coord_polar("y", direction = -1,start=55) +
  theme_void() +
  scale_fill_manual(values = pieColors) +
  theme(legend.position = "none", plot.background = element_rect(fill = 'transparent', color = NA))

plotGG(plot = pieVerdict, x = 8.75, y = 0.1, height = 1.3, width = 1.3)

plotText(label = sp_text[1], x = unit(10.2, "in"),
         y = 0.45,
         fontsize = 7, fontfamily = "Helvetica",
         just="left",
         fontcolor = "lightgrey",
         fontface = "bold")
plotText(label = "90.86%", x = unit(10.15, "in"),
         y = 0.45,
         fontsize = 7, fontfamily = "Helvetica",
         just="right",
         fontcolor = "lightgrey",
         fontface = "bold")


plotText(label = sp_text[2], x = unit(10.2, "in"),
         y = 0.45+0.15,
         fontsize = 7, fontfamily = "Helvetica",
         just="left",
         fontcolor = "#96bce3",
         fontface = "bold")
plotText(label = "6.88%", x = unit(10.15, "in"),
         y = 0.45+0.15,
         fontsize = 7, fontfamily = "Helvetica",
         just="right",
         fontcolor = "#96bce3",
         fontface = "bold")


plotText(label = sp_text[3], x = unit(10.2, "in"),
         y = 0.45+0.15*2,
         fontsize = 7, fontfamily = "Helvetica",
         just="left",
         fontcolor = "#3f85cc",
         fontface = "bold")
plotText(label = "1.03%", x = unit(10.15, "in"),
         y = 0.45+0.15*2,
         fontsize = 7, fontfamily = "Helvetica",
         just="right",
         fontcolor = "#3f85cc",
         fontface = "bold")


plotText(label =  sp_text[4], x = unit(10.2, "in"),
         y =  0.45+0.15*3,
         fontsize = 7, fontfamily = "Helvetica",
         just="left",
         fontcolor = "#085099",
         fontface = "bold")
plotText(label = "0.62%", x = unit(10.15, "in"),
         y =  0.45+0.15*3,
         fontsize = 7, fontfamily = "Helvetica",
         just="right",
         fontcolor = "#085099",
         fontface = "bold")

plotText(label =  sp_text[5], x = unit(10.2, "in"),
         y =  0.45+0.15*4,
         fontsize = 7, fontfamily = "Helvetica",
         just="left",
         fontcolor = "#043363",
         fontface = "bold")
plotText(label = "0.62%", x = unit(10.15, "in"),
         y =  0.45+0.15*4,
         fontsize = 7, fontfamily = "Helvetica",
         just="right",
         fontcolor = "#043363",
         fontface = "bold")

#------------------------------------------------------------------------------
#Plot C
#------------------------------------------------------------------------------
plotText("C", x = 6, y = 1.5, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

load(file="output/results_plots/Figure1_differntial_splicing/pathway_barplots.rda")

plotGG(pathway_barplots, x = 6.5, y = 1.65, width = 4, height = 2.5)

load(file="output/results_plots/Figure1_differntial_splicing/GO_barplots.rda")
plotGG(GO_barplots, x = 6.5, y = 4.1, width = 4, height = 2.5)


#-------------------------------------------------------------------------------
#Plot D
#PPI network
library(png)
plotText("D", x = 0.1, y = 5.25, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

PPI_test <- readPNG("output/results_plots/Figure1_differntial_splicing/CTL_FNF_deltapsi15_added_name.png")
plotRaster(image= PPI_test,
           x=0.1, y=5.4,width=6,height=4,
           just = c("top","left"))

#-------------------------------------------------------------------------------
plotText("E", x = 6, y = 6.5, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

load(file="output/results_plots/Supplementary_figures/Fig1E_OA_boxplot.rda")
plotGG(oa_boxplots, x = 7.0, y = 6.7, width = 3.7, height =2.8 )

dev.off()




