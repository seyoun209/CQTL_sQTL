# Making boxplots for the QTL for shared, PBS or FNF treated
## Author: Seyoun Byun
## Date: 03.08.2024
## Edited:03.28.2024
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
library(softmaxreg)

leafviz("output/clu_fnf/PBSvsFNF.Rdata")
pbs_fnf <- load("output/clu_fnf/PBSvsFNF.Rdata")
introns_gencode <- fread("tools/leafcutter/gencode_hg38_all_introns.bed.gz")
introns_gencode_renamed <- introns_gencode %>% dplyr::rename(chr = "V1", start = "V2", end = "V3")
#make normalized counts for the PSI

# I am running the psi branch of leafcutter
# Rscript tools/psi_leafcutter/leafcutter/scripts/leafcutter_quantify_psi.R output/clu_fnf/ctlvsfnf_perind_numers.counts.gz -c output/CTLvsFNF_group.txt -o output/clu_fnf/ratio_fnf_counfounder.txt.gz > output/logs/qunatify_psiFNF.err 2>&1

#This is before correcting the confounders

ratios_fnf <- counts_fnf %>%
  mutate(clu = str_split_fixed(rownames(counts_fnf), ":", 4)[,4]) %>%
  group_by(clu) %>%
  mutate(across(everything(), ~./sum(.))) %>%
  ungroup() %>%
  as.data.frame() %>%
  set_rownames(rownames(counts_fnf)) %>%
  dplyr::select(-clu)

ratios_fnf = ratios_fnf[rowMeans(is.na(ratios_fnf)) <= 0.4,,drop=F ] #Try to remove 40% or less NA values in each rows
#147127 to 144656 --> 2471 rows dropped
row_means_fnf = rowMeans(ratios_fnf, na.rm = T) #calculate the mean without NAs in the rows. 
row_means_outer_fnf = outer(row_means_fnf, rep(1,ncol(ratios_fnf)))# making outlier to the NAs
ratios_fnf[is.na(ratios_fnf)] = row_means_outer_fnf[is.na(ratios_fnf)] # instead of Na, add the rowmean
colnames(ratios_fnf) <- colnames(counts_fnf)
ratios_fnf <- cbind(rownames(ratios_fnf), ratios_fnf)
colnames(ratios_fnf)[1] <- c('Junction')
#write.table(ratios_fnf, file = "output/clu_fnf/ratio_fnf.txt",sep='\t',quote=F,row.names=F,col.names=T)

meta=read.table("output/CTLvsFNF_group.txt", header=F, stringsAsFactors = F)

if (!is.null("output/CTLvsFNF_group.txt")) {
  cat("Loading counfounders from","output/CTLvsFNF_group.txt","\n")
  if (!file.exists("output/CTLvsFNF_group.txt")) stop("File ","output/CTLvsFNF_group.txt"," does not exist")
  meta=read.table("output/CTLvsFNF_group.txt", header=F, stringsAsFactors = F)
  colnames(meta)[1]="sample"
  counts=counts[,meta$sample]
  
  confounders=meta[,2:ncol(meta),drop=F]
  # scale continuous confounders
  for (i in seq_len(ncol(confounders)))
    if (is.numeric(confounders[,i]))
      confounders[,i]=scale(confounders[,i])
  # convert factors to one-of-K encoding
  confounders=model.matrix( ~., data=confounders )
  confounders=confounders[,2:ncol(confounders),drop=F] # remove intercept
  
} else {
  confounders = matrix(0,nrow=ncol(counts),ncol=0)
}



#-------------------------------------------------------------------------------
#pca_before batch remove
config <- yaml::read_yaml("config/rna_prcoess.yaml")
conditions_to_include <- unlist(strsplit("CTL FNF", " "))

meta_cqtl <- fread("output/clu_fnf/meta_cqtl")
# Omit samples and filter rows based on configuration and conditions
meta_ctl_fnf <- meta_cqtl[!meta_cqtl$Donor %in% config$samples_to_omit, ] %>%
  filter(Condition %in% conditions_to_include)
psi_fnf_transposed <- t(ratios_fnf)
psi_fnf_transposed <-psi_fnf_transposed[-1,]
temp_df  <- as.data.frame(psi_fnf_transposed)
# Convert each column to numeric
temp_df[] <- lapply(temp_df, as.numeric)

# Optionally, convert back to a matrix if needed
psi_fnf_transposed <- as.matrix(temp_df)

pca_pre_fnf <- merge(psi_fnf_transposed ,meta_ctl_fnf, by.x = 'row.names',by.y="ID", all = TRUE)
pca_pre_fnf_prcomp <- prcomp(pca_pre_fnf[, !(colnames(pca_pre_fnf) %in% c("Row.names","Donor","ID", "Condition","Sex","Age","Predicted_Ancestry","FragmentBatch","RNAextractionKitBatch","RIN","RNAshippedDate"))], scale. = TRUE)

z.df <- as.data.frame(pca_pre_fnf_prcomp$x)
z.df$Condition <- pca_pre_fnf$Condition
z.df$Donor <- pca_pre_fnf$Donor

variance_explained <- summary(pca_pre_fnf_prcomp)$importance[2, ]

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
  


#------------------------------------------------------------------------------
#corrected PCA

psi_fnf_batch_corrected_transposed <- t(limma_psi_batchremove.df[,-ncol(limma_psi_batchremove.df)])
temp_fnf_batch_df  <- as.data.frame(psi_fnf_batch_corrected_transposed)

# Convert each column to numeric
temp_fnf_batch_df[] <- lapply(temp_fnf_batch_df, as.numeric)

# Optionally, convert back to a matrix if needed
psi_fnf_batch_transposed <- as.matrix(temp_fnf_batch_df)

pca_pre_fnf_batch <- merge(psi_fnf_batch_transposed ,meta_ctl_fnf, by.x = 'row.names',by.y="ID", all = TRUE)
pca_pre_fnf_batch_prcomp <- prcomp(pca_pre_fnf_batch[, !(colnames(pca_pre_fnf_batch) %in% c("Row.names","Donor","ID", "Condition","Sex","Age","Predicted_Ancestry","FragmentBatch","RNAextractionKitBatch","RIN","RNAshippedDate"))], scale. = TRUE)

z_batch.df <- as.data.frame(pca_pre_fnf_batch_prcomp$x)
z_batch.df$Condition <- pca_pre_fnf_batch$Condition
z_batch.df$Donor <- pca_pre_fnf_batch$Donor

variance_explained <- summary(pca_pre_fnf_batch_prcomp)$importance[2, ]

after_batch_fnf <- ggplot(z_batch.df, aes(x = PC1, y = PC2, color = Condition )) + #label = Donor, )) +
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
pdf(file = "output/results_plots/PCA_Plot_batchcorrection.pdf",   # The directory you want to save the file in
    width = 9, # The width of the plot in inches
    height = 6)
grid.arrange(before_batch_fnf, after_batch_fnf, ncol = 2, widths = c(2, 2))
dev.off()




#corrected PCA

ratios_fnf_batch_corrected_transposed <- t(ratios_fnf_batch_corrected[,-1])
ratios_fnf_batch_corrected_df  <- as.data.frame(ratios_fnf_batch_corrected_transposed)

# Convert each column to numeric
ratios_fnf_batch_corrected_df[] <- lapply(ratios_fnf_batch_corrected_df, as.numeric)

# Optionally, convert back to a matrix if needed
ratios_fnf_batch_corrected_transposed <- as.matrix(ratios_fnf_batch_corrected_df)
constant_or_zero_cols <- apply(pca_pre_fnf_batch[, !(colnames(pca_pre_fnf_batch) %in% c("Row.names","Donor","ID", "Condition","Sex","Age","Predicted_Ancestry","FragmentBatch","RNAextractionKitBatch","RIN","RNAshippedDate"))], 2, function(x) var(x) == 0 || all(x == 0))

# Names of columns to remove
cols_to_remove <- names(constant_or_zero_cols[constant_or_zero_cols])
pca_pre_fnf_batch <- merge(ratios_fnf_batch_corrected_transposed ,meta_ctl_fnf, by.x = 'row.names',by.y="ID", all = TRUE)
pca_pre_fnf_batch_prcomp <- prcomp(pca_pre_fnf_batch[, !(colnames(pca_pre_fnf_batch) %in% c("Row.names","Donor","ID", "Condition","Sex","Age","Predicted_Ancestry","FragmentBatch","RNAextractionKitBatch","RIN","RNAshippedDate", cols_to_remove))], scale. = TRUE)

z_batch.df <- as.data.frame(pca_pre_fnf_batch_prcomp$x)
z_batch.df$Condition <- pca_pre_fnf_batch$Condition
z_batch.df$Donor <- pca_pre_fnf_batch$Donor

variance_explained <- summary(pca_pre_fnf_batch_prcomp)$importance[2, ]

after_batch_fnf <- ggplot(z_batch.df, aes(x = PC1, y = PC2, color = Condition )) + #label = Donor, )) +
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
pdf(file = "output/results_plots/PCA_Plot_batchcorrection.pdf",   # The directory you want to save the file in
    width = 9, # The width of the plot in inches
    height = 6)
grid.arrange(before_batch_fnf, after_batch_fnf, ncol = 2, widths = c(2, 2))
dev.off()


#-------------------------------------------------------------------------------

# Finding the background genes
split_names <- strsplit(rownames(counts), ":")
split_df <- matrix(unlist(split_names), ncol = 4, byrow = TRUE) |> as.data.frame()
split_df <- split_df %>%
  mutate(start = as.integer(start),
         end = as.integer(end))
colnames(split_df) <- c("chr", "start", "end", "clusterID")
combined_df <- inner_join(split_df, introns_gencode_renamed, by = c("chr", "start", "end"))
combined_df_nodot <- combined_df %>%
  mutate(V5 = sub("\\..*$", "", V5))

subset_background_genes_ensg <- combined_df_nodot $V5 |> unique() 
#write.table(subset_background_genes_ensg, file = "output/clu_fnf/background_gene_set.txt",sep='\t',quote=F,row.names=F,col.names=F)


#Finding the Differential gene expression and plot (Cluster based)
# Finding the one introns per cluster that is the most different.


for (name in pbs_fnf) {
  # Construct the new name by prefixing with "fnf_"
  new_name <- paste0(name,"_fnf")
  
  # Assign the object to the new name in the global environment
  assign(new_name, get(name))
  remove(name)
}


fnf_maxCluster <- introns_fnf %>%
  mutate(abs_deltapsi = abs(deltapsi)) %>%  # Add a new column for the absolute value of deltapsi
  group_by(clusterID) %>%
  # Use slice_max to select the row with the maximum absolute deltapsi value
  slice(which.max(abs_deltapsi)) %>%
  # Ensure that ensemblID is not "."
  filter(ensemblID != ".") %>%
  # Optionally, remove the abs_deltapsi column if it's no longer needed
  dplyr::select(-abs_deltapsi)

sig_psi_maxCluster <- fnf_maxCluster[abs(fnf_maxCluster$deltapsi ) >= 0.2,]
sig_psi_maxCluster$loc <- paste0(sig_psi_maxCluster$chr,":",sig_psi_maxCluster$start,":",sig_psi_maxCluster$end,":",sig_psi_maxCluster$clusterID)
#Gene for the cluster Max
sig_psi_maxCluster_nodot <- sig_psi_maxCluster %>%
  mutate(ensemblID = sub("\\..*$", "", ensemblID))
sig_gene_20_percent_diff <- sig_psi_maxCluster_nodot$ensemblID |> unique()  #Save it to ENSG
#write.table(sig_gene_20_percent_diff, file = "output/clu_fnf/sig_gene_20_percent_diff.txt",sep='\t',quote=F,row.names=F,col.names=F)

#read the meta data for the Age, sex, Ancestry
load("output/combined_meta_data.RData") # It is loading name is combined_data

selected_columns <- c("Donor","ID","Condition","Sex", "Age") #select column needed it based
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
filtered_meta_samples <- meta_cqtl %>%
  filter(!grepl("OA", ID))
meta_samples <-filtered_meta_samples[,2:6] |> as.data.frame()
#-------------------------------------------------------------------------------
#Heatmap
rownames(meta_samples) <- meta_samples$ID
psi_subset <- ratios[rownames(ratios) %in% sig_psi_maxCluster$loc,]
psi_df <- cbind(rownames(psi_subset),psi_subset ) |> as.data.frame()
colnames(psi_df)[1] <- c("Junction")

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_subset_norm <- t(apply(psi_subset, 1, cal_z_score))

#z-score
brks <- seq(min(data_subset_norm, na.rm = TRUE), max(data_subset_norm, na.rm = TRUE), length = 51)
brks <- seq(-2,2,length.out=50) 

#data_subset_norm <- t(apply(psi_df[,-1], 1, scale))
#colnames(data_subset_norm) <- colnames(psi_df[,-1])
# Define the age group intervals
age_breaks <- c(30, 40, 50, 60, 70, 80, Inf)

# Labels for the age groups
age_labels <- c("30-39", "40-49", "50-59", "60-69", "70-79", "80+")

# Create age groups using the cut function
age_groups <- cut(meta_samples$Age, breaks = age_breaks, labels = age_labels, right = FALSE)
meta_samples$Age_range <- age_groups

age_palette <- c("#ffcba4","#cca283","#997a62","#806652","#665142","#332921")

# Create a named vector of colors for the age groups
age_colors <- setNames(age_palette, age_labels)


my_colour = list(
  Condition  = c("CTL" = "#9FCCE4", "FNF" = "#FAB394"),
  Sex  = c( "M" =  "#0075B0", "F" = "#F8B7CD"),
  Predicted_Ancestry= c("AMR"="#F5BC9F", "EUR"="#FAF1D2" ,"AFR"="#86CBB7","SAS"="#6B7EA4"),
  Age_range = age_colors)

pheatmap::pheatmap(data_subset_norm,
                   annotation_col = meta_samples,
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


sample_info <- meta_samples[, c(6,5, 3, 2)]
colnames(sample_info) <- c("Age","Ancestry","Sex","Condition")
my_colour = list(
  Condition  = c("CTL" = "#B8B8B8", "FNF" = "#4A4A4A"),
  Sex  = c("F" = "#DD8492", "M" = "#4788BA"),
  Ancestry= c("AFR"="#C74A53", "AMR"="#EFB06E","EUR"= "#5BAD58", "SAS"="#177F97"),
  Age = age_colors)


sample_info_ordered <- sample_info[colnames(data_subset_norm), ]

colAnn <- HeatmapAnnotation(df = sample_info_ordered,
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


column_order <- colnames(data_subset_norm)[order(sample_info_ordered$Condition)]

hmap <- Heatmap(
  data_subset_norm,
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

pdf(file = "output/results_plots/figure1.pdf",   # The directory you want to save the file in
    width = 11, # The width of the plot in inches
    height = 6.5)
pageCreate(width = 11.5, height = 9.5, showGuides = TRUE)
plotText("A", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

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





#---------------------------

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


dchart <- table(sig_psi_maxCluster$verdict) |> as.data.frame()
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
plotText(label = "92.2%", x = unit(10.15, "in"),
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
plotText(label = "5.1%", x = unit(10.15, "in"),
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
plotText(label = "1.5%", x = unit(10.15, "in"),
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
plotText(label = "0.6%", x = unit(10.15, "in"),
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
plotText(label = "0.6%", x = unit(10.15, "in"),
         y =  0.45+0.15*4,
         fontsize = 7, fontfamily = "Helvetica",
         just="right",
         fontcolor = "#043363",
         fontface = "bold")


#------------------------------------------------------------------------------
#Pathway



#-------------------------------------------------------------------------------
#GO



dev.off()




