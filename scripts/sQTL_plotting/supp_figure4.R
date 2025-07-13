# Supp figure 4 (Correlation plot) 
setwd("/work/users/s/e/seyoun/CQTL_sQTL/output")
library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table)
library(yaml)
library(corrplot)
library(grid)
library(ggpubr)
library(ggforce)

# Function for the correlation 
create_correlation_plot <- function(pca_df, title) {
  # Separate the non-PC columns (predictors) and the PC columns (responses)
  predic <- pca_df[, grep("PC", colnames(pca_df), invert = TRUE)]
  PCs <- pca_df[, grep("PC", colnames(pca_df))]
  
  # Calculate p-values
  conflevel = cor.mtest(pca_df, conf.level = 0.95)
  
  # Calculate correlations
  correlations <- sapply(1:ncol(predic), function(i) {
    cor(PCs, pca_df[, i], use = "complete.obs", method = "pearson")
  })
  colnames(correlations) <- colnames(predic)
  
  # Reshape the correlation matrix
  melted_corr <- melt(correlations)
  names(melted_corr) <- c("x", "y", "cor")
  
  # Convert x to PC format to match melted_p
  melted_corr$x <- paste0("PC", melted_corr$x)
  
  # Reshape the p-value matrix
  melted_p <- melt(conflevel$p[grep("PC", colnames(pca_df)),
                               grep("PC", colnames(pca_df), invert = TRUE)])
  names(melted_p) <- c("x", "y", "pval")
  
  # Combine correlation and p-value data
  all_data <- merge(melted_corr, melted_p, by = c("x", "y"))
  
  # Add significance labels
  all_data$pval_sig <- ifelse(all_data$pval < 0.001, "***",
                              ifelse(all_data$pval < 0.01, "**",
                                     ifelse(all_data$pval < 0.05, "*", "")))
  
  # Format correlation values
  all_data$cor_text <- sprintf("%.2f", all_data$cor)
  
  # Ensure correct ordering of PCs for y-axis only
  pc_order <- paste0("PC", 10:1)
  all_data$x <- factor(all_data$x, levels = pc_order)
  
  # Create the plot
  plot <- ggplot(all_data, aes(x = y, y = x, fill = cor)) +
    geom_tile() +
    scale_fill_gradient2(limits = c(-1, 1), low = "#418C82", mid = "white", high = "#C55E2D") +
    geom_text(aes(label = cor_text), family = "Helvetica", size = 2) +
    geom_text(aes(label = pval_sig), family = "Helvetica", size = 2, fontface = "bold",
              vjust = -0.75) +
    guides(fill = guide_colorbar(title = "",
                                 title.position = "top",
                                 title.hjust = 1, 
                                 direction = "vertical")) +
    theme(axis.text.x = element_text(angle = 45, size = 7, vjust = 1, hjust = 0.8,
                                     color = "black"),
          axis.text.y = element_text(size = 7, color = "black"),
          panel.background = element_rect(fill = 'transparent', color = "transparent"),
          plot.background = element_rect(fill = 'transparent', color = "transparent"),
          strip.background = element_blank(),
          text = element_text(family = "Helvetica"),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.title = element_text(size = 4, color = "black"),
          legend.text = element_text(size = 4, color = "black"),
          legend.position = "right",
          legend.key.width = unit(0.6, 'mm'),
          legend.margin = margin(t = -35),
          strip.text = element_text(size = 10),
          title = element_blank())
    #ggtitle(paste("Pearson Correlation -", title))
  
  return(plot)
}



# samplesheet
aligned_samplesheet_path <- "../aligned_samplesheet.txt"
donor_samples_path <- "../donor_samples.txt"
rna_extraction_path <- "../rna_extraction.txt"
conditions_to_include <- unlist(strsplit("CTL FNF OA", " ")) # Split the first argument into condition values
#use_wasp_id <- tail(args, 1) == "wasp" # Check if the last argument is "wasp"


# Load configuration from YAML file
config <- yaml::read_yaml("../config/rna_prcoess.yaml")
aligned_samples <- fread(aligned_samplesheet_path)
donor_samples <- fread(donor_samples_path)
rna_extraction <- fread(rna_extraction_path)

# Identify common columns between aligned_samples and donor_samples, excluding 'Donor'
common_cols_donor <- setdiff(intersect(colnames(aligned_samples), colnames(donor_samples)), "Donor")

# Remove these common columns from donor_samples before the join
donor_samples_clean <- dplyr::select(donor_samples, -all_of(common_cols_donor))

# Perform the first left join
combined_data <- aligned_samples %>%
  left_join(donor_samples_clean, by = "Donor")

common_cols_rna <- setdiff(intersect(colnames(combined_data), colnames(rna_extraction)), "Read2")
rna_extraction_clean <- dplyr::select(rna_extraction, -all_of(common_cols_rna))
combined_data <- combined_data %>%
  left_join(rna_extraction_clean, by = "Read2")

# Process ID column with and without "_wasp"
combined_data <- combined_data %>%
  mutate(ID = paste(Donor, Condition, Tech_Rep, Sex, sep = "_"),
         ID_wasp = paste(Donor, Condition, Tech_Rep, Sex, "wasp", sep = "_"))

# Conditionally adjust FragmentBatch for 'OA' condition
combined_data <- combined_data %>%
  mutate(FragmentBatch = ifelse(Condition == "OA", 0, FragmentBatch))

# Omit samples and filter rows based on configuration and conditions
combined_data <- combined_data[!combined_data$Donor %in% config$samples_to_omit, ] %>%
  dplyr::filter(Condition %in% conditions_to_include)
save(combined_data, file = "combined_meta_data.RData")
# Select and rename columns dynamically, based on the use of wasp ID
selected_columns <- c("ID","Condition","Sex", "Age","Race","OAGradeAvg","CauseOfDeath","FragmentBatch","RIN","RNAextractionKitBatch","RNAshippedDate")
selected_columns_wasp <- c("ID_wasp","Condition","Sex", "Age","Race","OAGradeAvg","CauseOfDeath","FragmentBatch","RIN","RNAextractionKitBatch","RNAshippedDate")
final_meta_data <- combined_data %>% dplyr::select(all_of(selected_columns))
meta_cqtl <- fread("clu_fnf/meta_cqtl")
meta_cqtl_ancetry_only <- meta_cqtl %>% dplyr::select(ID,Predicted_Ancestry)
final_meta_all <- left_join(final_meta_data,meta_cqtl_ancetry_only,by="ID")
final_meta_data_wasp <- combined_data %>% dplyr::select(all_of(selected_columns_wasp))

#pca for ctl vs fnf
pca_raw <- fread("./clu_fnf/ctlvsfnf_perind.counts.gz.PCs") |> 
  t() |> 
  as.data.frame()
pca_raw <- pca_raw[-1,]
colnames(pca_raw) <- c(sprintf("PC%s",seq(1:20)))
pca_raw_nm <- cbind(rownames(pca_raw),pca_raw)
colnames(pca_raw_nm)[1] <-c("ID")
pc10_ctl_fnf <- pca_raw_nm[,1:11]

final_meta_all_fixed <- final_meta_all %>% dplyr::select(-c("Race"))
#Merge with the 
splicingPCA_df <- merge(final_meta_all_fixed, pc10_ctl_fnf, by="ID" , all=FALSE)

#Change to factor

# Assuming splicingPCA_df is your dataframe
for (col in colnames(splicingPCA_df)) {
  # Convert only columns starting with "PC" to numeric
  if (!grepl("^PC", col)) {
    splicingPCA_df[[col]] <-  as.factor(splicingPCA_df[[col]])
  }
}
colnames(splicingPCA_df)[which(colnames(splicingPCA_df) == 'RNAshippedDate')] <- c("SeqeuncingBatch")
splicingPCA_numeric <- sapply(splicingPCA_df,as.numeric) |> as.data.frame() # Convert columns to numeric, assuming first column is ID

# This step is done when there is pbs_corrPlot and then save the legend and save the plot without the 
#legend <- cowplot::get_legend(pbs_corrPlot)
#legend_corr <- cowplot::ggdraw() + cowplot::draw_grob(legend)
#save(legend_corr, file="output/results_plots/sqtl_plots/legend_corr.rda")

# PBS

# For PBS condition
ctl_pca_df <- splicingPCA_numeric %>% 
  dplyr::filter(Condition == "1") %>% 
  dplyr::select(-Condition, -ID)

pbs_corrPlot <- create_correlation_plot(ctl_pca_df, "PBS")


# For FNF condition
fnf_pca_df <- splicingPCA_numeric %>% 
  dplyr::filter(Condition == "2") %>% 
  dplyr::select(-Condition, -ID)

# Save the full plot with legend
#legend <- cowplot::get_legend(pbs_corrPlot)
#legend_corr <- cowplot::ggdraw() + cowplot::draw_grob(legend)
#save(legend_corr, file="results_plots/sqtl_plots/legend_corr.rda")

pbs_corrPlot <- create_correlation_plot(ctl_pca_df, "PBS") 
#+ theme(legend.position = "none")
save(pbs_corrPlot, file = "results_plots/sqtl_plots/pbs_corrPlot.rda")

fnf_corrPlot <- create_correlation_plot(fnf_pca_df, "FNF") 
save(fnf_corrPlot, file = "results_plots/sqtl_plots/fnf_corrPlot.rda")


# This is for the Supp figure 4d (rank plot)

response_pbs_results <- readRDS("01.qtltools_re/conditional_pbs/response_pbs_resInv.rds")
response_fnf_results <- readRDS("01.qtltools_re/conditional_pbs/response_fnf_resInv.rds")

rank_counts <- function(counts, dataset_name) {
  df <- counts |> table() |> as_tibble()  |> dplyr::rename(rank="counts") |> dplyr::rename( counts ="n") %>%
    dplyr::mutate(dataset_name=dataset_name) %>% dplyr::rename(group="dataset_name")
  return(df)
}

# Prepare data for both PBS and FNF
all_cond_count_rank.df <- rbind(rank_counts(response_pbs_results$rank, "PBS"), rank_counts(response_fnf_results$rank, "FN-f"))
all_cond_count_rank <- all_cond_count_rank.df %>%
  mutate(rank = paste0(rank, "°"))
all_cond_count_rank <- all_cond_count_rank %>%
  mutate(rank = factor(rank, levels = c("3°", "2°", "1°", "0°"), ordered = TRUE)) %>%
  mutate(group = factor(group, levels = c("PBS","FN-f"), ordered = TRUE))
#dfm_only_cond <- dfm %>% dplyr::filter(rank != "0°")


replace_zero <- function(x) {
  ifelse(is.na(x) | x == 0, 0.01, x)
}

 all_cond_count_rank_prepared <- all_cond_count_rank %>%
  complete(rank, group, fill = list(counts = NA)) %>%
  mutate(counts = replace_zero(counts),
         group = factor(group, levels = c("PBS", "FN-f")))  # Ensure PBS comes first

rankslog2_barplot <-ggplot(all_cond_count_rank_prepared, aes(x = rank, y = counts, fill = group)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  scale_fill_manual(values = c("PBS" = "#BFDDFF", "FN-f" = "#FFDDA2")) +
  scale_y_continuous(
    name = "Counts of sQTL-splice intron junction",
    trans = scales::pseudo_log_trans(base = 2),
    breaks = c(0, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096),
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_x_discrete(
    name = "Number of Independent Signals",
    limits = c("0°", "1°", "2°", "3°")
  ) +
  coord_cartesian(clip = "off") +
  theme(
    strip.placement = "outside",
    axis.line = element_line(linewidth = 0.25),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black", linewidth = 0.25),
    axis.ticks.length.y = unit(-0.1, "cm"),
    axis.title.x = element_markdown(size = 8, family = "Helvetica", margin = margin(t = 5)),
    axis.title.y = element_markdown(size = 8, family = "Helvetica", margin = margin(r = 5)),
    text = element_text(family = "Helvetica"),
    axis.text.y = element_text(color = "black", size = 6),
    axis.text.x = element_text(color = "black", size = 6, margin = margin(t = 5)),
    strip.background = element_blank(),
    strip.text.x.top = element_text(size = 8, margin = margin(b = 5)),
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.9, 0.9),
    panel.spacing.y = unit(0.5, "cm")
  )

save(rankslog2_barplot, file = "results_plots/sqtl_plots/rankslog2_barplot.rda")
#------------------------------------------------------------------------------
# Venndiagram for the response QTL 

# Subset significant reponse sQTL only
pbs_resQTL <- response_pbs_results %>%
  dplyr::filter(interaction_pval  < 0.05) |>
  arrange(interaction_pval)

fnf_resQTL <- response_fnf_results %>%
  dplyr::filter(interaction_pval  < 0.05) |>
  arrange(interaction_pval)

# Subset high confidence re-sQTls ####

pbs_highConf_resQtL <- response_pbs_results %>% 
  dplyr::filter(interaction_pval  < 0.05) %>%
  dplyr::filter(abs(delta_beta ) > 0.2) |>
  dplyr::filter(minor_alle_count >= 5) |>
  arrange(interaction_pval)

fnf_highConf_resQtL <- response_fnf_results %>% 
  dplyr::filter(interaction_pval  < 0.05) %>%
  dplyr::filter(abs(delta_beta ) > 0.2) |>
  dplyr::filter(minor_alle_count >= 5) |>
  arrange(interaction_pval)

pbs_specific_resQTL <- pbs_highConf_resQtL |> 
  dplyr::filter(!clusterID %in% fnf_highConf_resQtL$clusterID ) 

fnf_specific_resQTL <- fnf_highConf_resQtL |> 
  dplyr::filter(!clusterID %in% pbs_highConf_resQtL$clusterID )  

# pieChart  fro the significant response sQTL   

resqtl_tibble <- tibble(values = unique(c(pbs_resQTL$ensg,fnf_resQTL$ensg))) %>%
  mutate(PBS = values %in% pbs_resQTL$ensg,
         FNF = values %in% fnf_resQTL$ensg)

resqtl_intron_tibble <- tibble(values = unique(c(pbs_resQTL$phe_id,fnf_resQTL$phe_id))) %>%
  mutate(PBS = values %in% pbs_resQTL$phe_id,
         FNF = values %in% fnf_resQTL$phe_id)

# sqtl cluster counts
resqtl_cluster_pbs <-  pbs_resQTL  %>% dplyr::filter(!is.na(clusterID)) %>% distinct(clusterID)  %>% unlist() 
resqtl_cluster_fnf <-  fnf_resQTL  %>% dplyr::filter(!is.na(clusterID)) %>% distinct(clusterID)  %>% unlist() 

#sQTL cluste tibble

resqtl_clusters <- tibble (values = unique(c(resqtl_cluster_pbs,resqtl_cluster_fnf ))) %>%
  mutate(PBS = values %in% resqtl_cluster_pbs,
         FNF = values %in% resqtl_cluster_fnf)

#------------------------------------------------------------------------------
#Re-make the sQtl venn diagram for the sIntrons 
#Intron junctions -2

pbs_sGene_introns <- response_pbs_results %>%
  dplyr::filter(!is.na(genomicLoc)) %>%
  dplyr::select(genomicLoc) %>%
  unique() %>%
  unlist()

fnf_sGene_introns <- response_fnf_results %>%
  dplyr::filter(!is.na(genomicLoc)) %>%
  dplyr::select(genomicLoc) %>%
  unique() %>%
  unlist()



sig_introns_all_tibble <- tibble(values = unique(c(pbs_sGene_introns, fnf_sGene_introns))) %>%
  mutate(PBS = values %in% pbs_sGene_introns,
         FNF = values %in% fnf_sGene_introns)


pbs_sGene_clusters <- response_pbs_results %>%
  dplyr::filter(!is.na(clusterID)) %>%
  dplyr::select(clusterID) %>%
  unique() %>%
  unlist()

fnf_sGene_clusters <- response_fnf_results %>%
  dplyr::filter(!is.na(clusterID)) %>%
  dplyr::select(clusterID) %>%
  unique() %>%
  unlist()

create_venn_diagram <- function(total_PBS, total_FNF, overlap, only_PBS, only_FNF, max_count) {
  ggplot() +
    geom_circle(aes(x0 = -0.6, y0 = 0, r = sqrt(total_PBS/max_count)), 
                fill = "#BFDDFF", color = NA, alpha = 0.5) +
    geom_circle(aes(x0 = 0.6, y0 = 0, r = sqrt(total_FNF/max_count)), 
                fill = "#FFDDA2", color = NA, alpha = 0.5) +
    geom_text(aes(x = -0.8, y = 0, label = only_PBS), size = 3) +
    geom_text(aes(x = 0.8, y = 0, label = only_FNF), size = 3) +
    geom_text(aes(x = 0.1, y = 0, label = overlap), size = 3) +
    coord_fixed() +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    )
}


#-------------------------------------------------------------------------------
pdf(file = "./results_plots/sqtl_plots/fig4_supp_edited.pdf", width = 10.5, height = 9)
pageCreate(width = 10.5, height= 9 , default.units = "inches", showGuides = TRUE)

plotText("a", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
load("results_plots/sqtl_plots/pbs_corrPlot.rda")
load("results_plots/sqtl_plots/fnf_corrPlot.rda")

plotGG(pbs_corrPlot, x = 0.4, y = 0.5, width = 3.35, height = 3.5)
plotGG(fnf_corrPlot, x = 0.4, y = 4, width = 3.35, height = 3.5)

plotText("PBS", x = 2, y = 0.5, just = "center", fontfamily = "Helvetica",
         fontsize = 8, fontcolor ="#005587",fontface="bold" )

plotText("FN-f", x = 2, y = 4, just = "center", fontfamily = "Helvetica",
         fontsize = 8, fontcolor ='darkorange',fontface="bold" )

#load("results_plots/sqtl_plots/legend_corr.rda")

#plotGG(legend_corr, x = 0.5, y = 7.75, width = 2.5, height = 0.3)

# Figure 2b 

plotText("b", x = 4.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

load("results_plots/sqtl_plots/counts_dotplot.rda")

plotGG(counts_dotPlot,x = 4.45, y = 0.5, width = 5, height = 2.75)


# Figure 3c  venndiagrams
plotText("c", x = 4.1, y =3.65 , just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

plotText(label =  "sIntron-junctions", 
         x = 5.3,
         y= 3.9,
         fontsize = 8, fontfamily = "Helvetica",
         just="center", fontface = "bold",
         fontcolor = "black")

total_PBS <- length(pbs_sGene_introns)
total_FNF <- length(fnf_sGene_introns)
overlap <- length(intersect(pbs_sGene_introns, fnf_sGene_introns))
only_PBS <- total_PBS - overlap
only_FNF <- total_FNF - overlap
max_count <- max(total_PBS, total_FNF)

sIntrons_venn <- create_venn_diagram(total_PBS, total_FNF, overlap, only_PBS, only_FNF, max_count)

plotGG(plot =  sIntrons_venn, x =4, y = 3.75, height = 2, width = 2.25)

plotText(label =  "PBS", 
         x = 5.1,
         y= 4.1,
         fontsize = 8, fontfamily = "Helvetica",
         just="right",
         fontcolor = "#2057A7")

plotText(label =  "Fn-f", 
         x =5.45,
         y= 4.1,
         fontsize = 8, fontfamily = "Helvetica",
         just="left",
         fontcolor = "#F2BC40")

#------------------------------------------------------------------------------
total_PBS <- length(pbs_sGene_clusters)
total_FNF <- length(fnf_sGene_clusters)
overlap <- length(intersect(pbs_sGene_clusters, fnf_sGene_clusters))
only_PBS <- total_PBS - overlap
only_FNF <- total_FNF - overlap
max_count <- max(total_PBS, total_FNF)

sClusters_venn <- create_venn_diagram(total_PBS, total_FNF, overlap, only_PBS, only_FNF, max_count)
plotGG(plot =  sClusters_venn, x = 4, y =5.55 , height = 2, width = 2.25)
plotText(label = "sClusters", 
         x = 5.3,
         y=5.7,
         fontsize = 8, fontfamily = "Helvetica",
         just="center",fontface = "bold",
         fontcolor = "black")

plotText(label =  "PBS", 
         x = 5.05,
         y= 5.85,
         fontsize = 8, fontfamily = "Helvetica",
         just="right",
         fontcolor = "#2057A7")

plotText(label =  "Fn-f", 
         x =5.45,
         y= 5.85,
         fontsize = 8, fontfamily = "Helvetica",
         just="left",
         fontcolor = "#F2BC40")

# Figure 3d  venndiagrams
#Response Venndiagram ---------------------------------------------------------
plotText("d", x = 6.25, y =3.65 , just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

total_PBS <- sum(resqtl_tibble$PBS)
total_FNF <- sum(resqtl_tibble$FNF)
overlap <- sum(resqtl_tibble$PBS & resqtl_tibble$FNF)
only_PBS <- total_PBS - overlap
only_FNF <- total_FNF - overlap
max_count <- max(total_PBS, total_FNF)

resQTL_venn <- create_venn_diagram(total_PBS, total_FNF, overlap, only_PBS, only_FNF, max_count)
plotGG(plot =  resQTL_venn, x = 6.2, y =3.85 , height = 2, width = 2.25)
plotText(label = "re-sGene",
         x = 7.35,
         y=3.9,
         fontsize = 8, fontfamily = "Helvetica",
         just="center",fontface = "bold",
         fontcolor = "black")

plotText(label =  "PBS",
         x = 7.1,
         y= 4.2,
         fontsize = 8, fontfamily = "Helvetica",
         just="right",
         fontcolor = "#2057A7")

plotText(label =  "Fn-f",
         x =7.5,
         y= 4.2,
         fontsize = 8, fontfamily = "Helvetica",
         just="left",
         fontcolor = "#F2BC40")

#-------------------------------------------------------------------------------

total_PBS <- sum(resqtl_intron_tibble$PBS)
total_FNF <- sum(resqtl_intron_tibble$FNF)
overlap <- sum(resqtl_intron_tibble$PBS & resqtl_intron_tibble$FNF)
only_PBS <- total_PBS - overlap
only_FNF <- total_FNF - overlap
max_count <- max(total_PBS, total_FNF)

resQTL_intron_venn <- create_venn_diagram(total_PBS, total_FNF, overlap, only_PBS, only_FNF, max_count)

plotGG(plot =  resQTL_intron_venn, x = 6.2, y =5.55 , height = 2, width = 2.25)
plotText(label = "re-sIntrons", 
         x = 7.35,
         y=5.7,
         fontsize = 8, fontfamily = "Helvetica",
         just="center",fontface = "bold",
         fontcolor = "black")

plotText(label =  "PBS", 
         x = 7.1,
         y= 5.85,
         fontsize = 8, fontfamily = "Helvetica",
         just="right",
         fontcolor = "#2057A7")

plotText(label =  "Fn-f", 
         x =7.5,
         y= 5.85,
         fontsize = 8, fontfamily = "Helvetica",
         just="left",
         fontcolor = "#F2BC40")
#-------------------------------------------------------------------------------
# re-sQTL cluster counts

total_PBS <- sum(resqtl_clusters$PBS)
total_FNF <- sum(resqtl_clusters$FNF)
overlap <- sum(resqtl_clusters$PBS & resqtl_clusters$FNF)
only_PBS <- total_PBS - overlap
only_FNF <- total_FNF - overlap
max_count <- max(total_PBS, total_FNF)

resQTL_cluster_venn <- create_venn_diagram(total_PBS, total_FNF, overlap, only_PBS, only_FNF, max_count)

plotGG(plot =  resQTL_cluster_venn, x=8.2, y =3.85 , height = 2, width = 2.25)
plotText(label = "re-sClusters", 
         x = 9.5,
         y=3.9,
         fontsize = 8, fontfamily = "Helvetica",
         just="center",fontface = "bold",
         fontcolor = "black")

plotText(label =  "PBS", 
         x = 9.0,
         y= 4.2,
         fontsize = 8, fontfamily = "Helvetica",
         just="right",
         fontcolor = "#2057A7")

plotText(label =  "Fn-f", 
         x =9.75,
         y= 4.2,
         fontsize = 8, fontfamily = "Helvetica",
         just="left",
         fontcolor = "#F2BC40")

#-------------------------------------------------------------------------------
plotText("e", x = 8.1, y =3.65 , just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
load("results_plots/sqtl_plots/rankslog2_barplot.rda")
plotGG(plot =  rankslog2_barplot, x = 8.2, y =4.5, height = 2.5, width = 2.25)
dev.off()
