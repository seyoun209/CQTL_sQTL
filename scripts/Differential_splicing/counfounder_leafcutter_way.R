# This is the leafcutter finding thecounfouders
## Author: Seyoun Byun
## Date: 06.11.2024
## Edited:
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

# This is the leafcutter finding the counfouder 

# I am running the psi branch of leafcutter
# for this I use r/3.6
#module purge
#module add r/3.6.0

#help="One+K column file: 1. sample names (must match column names in counts_file), 
#2. Additional columns are used to specify confounders, e.g. batch/sex/age. Numeric columns will be treated as continuous, so use e.g. batch1, batch2, batch3 rather than 1, 2, 3 if you want a categorical variable.")
#awk '{for (i=2; i<=NF; i++) {$i = "batch" $i} print}' CTLvsFNF_group.txt | awk '{$2=""; print $0}' > CTLvsFNF_confounder.txt
#awk '{for (i=2; i<=NF; i++) {$i = "batch" $i} print}' CTLvsOA_group.txt | awk '{$2=""; print $0}'  > CTLvsOA_confounder.txt

#.libPaths(c("/nas/longleaf/home/seyoun/R/x86_64-pc-linux-gnu-library/3.6", "/nas/longleaf/rhel8/apps/r/3.6.0/lib64/R/library"))
#Rscript tools/psi_leafcutter/leafcutter/scripts/leafcutter_quantify_psi.R output/clu_fnf/ctlvsfnf_perind_numers.counts.gz -c output/CTLvsFNF_confounder.txt -o output/clu_fnf/ratio_fnf_batch_corrected.txt.gz > output/logs/qunatify_psiFNF.err 2>&1 &
#Rscript tools/psi_leafcutter/leafcutter/scripts/leafcutter_quantify_psi.R output/clu_oa/ctlvsoa_perind_numers.counts.gz -c output/CTLvsOA_confounder.txt -o output/clu_oa/ratio_oa_batch_corrected.txt.gz > output/logs/qunatify_psiOA.err 2>&1 

#This is before correcting the confounders

for (name in pbs_fnf) {
  # Construct the new name by prefixing with "fnf_"
  new_name <- paste0(name,"_fnf")
  
  # Assign the object to the new name in the global environment
  assign(new_name, get(name))
  remove(name)
}
normalize_column <- function(x) {
  x / sum(x, na.rm = TRUE)
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
#corrected confounder
#-------------------------------------------------------------------------------
ctl_fnf_ratio <- fread("output/clu_fnf/ratio_fnf_counfounder.txt.gz") |> as.data.frame()
rownames(ctl_fnf_ratio) <- ctl_fnf_ratio$V1
colnames(ctl_fnf_ratio)[1] <- "Junction"

#-------------------------------------------------------------------------------
#pca_before batch remove
config <- yaml::read_yaml("config/rna_prcoess.yaml")
conditions_to_include <- unlist(strsplit("CTL FNF", " "))

meta_cqtl <- fread("output/clu_fnf/meta_cqtl")
# Omit samples and filter rows based on configuration and conditions
meta_ctl_fnf <- meta_cqtl[!meta_cqtl$Donor %in% config$samples_to_omit, ] %>%
  dplyr::filter(Condition %in% conditions_to_include)

#In this part, the confounder already made it to  calculated as it is. Therefore, remove from the exist for the 40% below. 
filtered_ratios_fnf_adj <- ctl_fnf_ratio %>%
  filter(Junction %in% ratios_fnf_qc$Junction)
ratios_fnf_adj <- filtered_ratios_fnf_adj[,-1]
psi_fnf_transposed <- t(ratios_fnf_adj)

pca_pre_fnf <- merge(psi_fnf_transposed ,meta_ctl_fnf, by.x = 'row.names',by.y="ID", all = TRUE)
pca_pre_fnf_subset <- pca_pre_fnf[, !(colnames(pca_pre_fnf) %in% c("Row.names","Donor","ID", "Condition","Sex","Age",
                                                                   "Predicted_Ancestry","FragmentBatch",
                                                                   "RNAextractionKitBatch","RIN","RNAshippedDate"))]

zero_variance_columns <- sapply(pca_pre_fnf_nobatch, function(col) var(col) == 0)

# Print zero-variance columns if you want to see them
print(names(pca_pre_fnf_nobatch)[zero_variance_columns])
pca_pre_fnf_prcomp <- prcomp(pca_pre_fnf_subset, scale. = TRUE)

z.df <- as.data.frame(pca_pre_fnf_prcomp$x)
z.df$Condition <- pca_pre_fnf$Condition
z.df$Donor <- pca_pre_fnf$Donor

variance_explained <- summary(pca_pre_fnf_prcomp)$importance[2, ]


#leafcutter fixed batch
# pdf(file = "output/results_plots/PCA_leafcutter_correction.pdf",   # The directory you want to save the file in
#     width = 9, # The width of the plot in inches
#     height = 6)
after_batch_fnf<-  ggplot(z.df, aes(x = PC1, y = PC2, color = Condition )) + #label = Donor, )) +
  geom_point() +
  #geom_text_repel(size = 2, box.padding = unit(0.1, "lines")) +
  scale_color_manual(values = c("CTL" = "#74CCEC", "FNF" = "#FAB394")) +
  labs(
    title = "After batch correction",
    x = paste0("PC1 Variance:", round(variance_explained[1] * 100, 2), "%"),
    y = paste0("PC2 Variance:", round(variance_explained[2] * 100, 2), "%")) +
  theme(plot.background = element_rect(fill='transparent', color=NA))+
  theme_bw(base_size = 10) +
  theme(panel.background = element_rect(fill = "white")) +
  theme(aspect.ratio = 1,panel.grid = element_blank()) 



#no batch

ratios_fnf_qc_num <- ratios_fnf_qc[,-1]
psi_fnf_t_nobatch <- t(ratios_fnf_qc_num)

pca_pre_fnf_nobatch <- merge(psi_fnf_t_nobatch ,meta_ctl_fnf, by.x = 'row.names',by.y="ID", all = TRUE)
pca_pre_fnf_nobatch_subset <- pca_pre_fnf_nobatch[, !(colnames(pca_pre_fnf_nobatch) %in% c("Row.names","Donor","ID", "Condition","Sex","Age",
                                                                                           "Predicted_Ancestry","FragmentBatch",
                                                                                           "RNAextractionKitBatch","RIN","RNAshippedDate"))]
zero_variance_columns <- which(apply(pca_pre_fnf_nobatch_subset, 2, var) == 0)
pca_pre_fnf_nobatch_subset_no0variance <- pca_pre_fnf_nobatch_subset[, -zero_variance_columns]
pca_pre_fnf_nobatch_prcomp <- prcomp(pca_pre_fnf_nobatch_subset_no0variance, scale. = TRUE)

z.df <- as.data.frame(pca_pre_fnf_nobatch_prcomp$x)
z.df$Condition <- pca_pre_fnf_nobatch$Condition
z.df$Donor <- pca_pre_fnf_nobatch$Donor

variance_explained <- summary(pca_pre_fnf_nobatch_prcomp)$importance[2, ]


before_batch_fnf <- ggplot(z.df, aes(x = PC1, y = PC2, color = Condition )) + #label = Donor, )) +
  geom_point() +
  #geom_text_repel(size = 2, box.padding = unit(0.1, "lines")) +
  scale_color_manual(values = c("CTL" = "#74CCEC", "FNF" = "#FAB394")) +
  labs(
    title = "Before batch correction",
    x = paste0("PC1 Variance:", round(variance_explained[1] * 100, 2), "%"),
    y = paste0("PC2 Variance:", round(variance_explained[2] * 100, 2), "%")) +
  theme(plot.background = element_rect(fill='transparent', color=NA))+
  theme_bw(base_size = 10) +
  theme(panel.background = element_rect(fill = "white")) +
  theme(aspect.ratio = 1,panel.grid = element_blank())

pdf(file = "output/results_plots/PCA_Plot_leafcutter_Confounder.pdf",   # The directory you want to save the file in
    width = 9, # The width of the plot in inches
    height = 6)
grid.arrange(before_batch_fnf, after_batch_fnf, ncol = 2, widths = c(2, 2))
dev.off()