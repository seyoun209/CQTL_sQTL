setwd("/work/users/s/e/seyoun/CQTL_sQTL/output")
library(GGally)
library(ggplot2)
library(corrplot)
library(dplyr)
library(yaml)
library(readr)
library(data.table)
library(grid)
library(gridExtra)
library(cowplot)
library(gridGraphics)
source("../scripts/utils/utils.R")

#For the corrplot star as significant do this(https://stackoverflow.com/questions/63227830/r-corrplot-plot-correlation-coefficients-along-with-significance-stars_:
#trace(corrplot, edit=TRUE)
# and the find the line 443
# adjust text(X,Y ...) according to your needs, here +0.25 is added to the Y-position    
#place_points = function(sig.locs, point) {
#  text(pos.pNew[, 1][sig.locs], (pos.pNew[, 2][sig.locs])+0.25, 
#       labels = point, col = pch.col, cex = pch.cex, 
#       lwd = 2)
# and then Save


#args <- commandArgs(trailingOnly = TRUE)
#aligned_samplesheet_path <- args[1]
#donor_samples_path <- args[2]
#rna_extraction_path <- args[3]
#conditions_to_include <- unlist(strsplit(args[4], " ")) # Split the first argument into condition values
#use_wasp_id <- tail(args, 1) == "wasp" # Check if the last argument is "wasp"


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


#Seperate out the condition
# Separate the dataframe into two based on the Condition column

colnames(splicingPCA_df)[which(colnames(splicingPCA_df) == 'RNAshippedDate')] <- c("SeqeuncingBatch")
splicingPCA_numeric <- sapply(splicingPCA_df,as.numeric) |> as.data.frame() # Convert columns to numeric, assuming first column is ID

# Subset for CTL condition and drop Condition column
ctl_pca_df <- splicingPCA_numeric %>% subset(Condition == "1") %>% dplyr::select(-Condition)

# Subset for FNF condition and drop Condition column
fnf_pca_df <- splicingPCA_numeric %>% subset(Condition == "2") %>% dplyr::select(-Condition)

# Separate the non-PC columns (predictors) and the PC columns (responses)
predic_ctl <- ctl_pca_df[, grep("PC", colnames(ctl_pca_df), invert = TRUE)]
PCs_ctl <- ctl_pca_df[, grep("PC", colnames(ctl_pca_df))]

predict_fnf <- fnf_pca_df[,  grep("PC", colnames(fnf_pca_df), invert = TRUE)]
Pcs_fnf <- fnf_pca_df[, grep("PC", colnames(fnf_pca_df))]

conflevel.ctl = cor.mtest(ctl_pca_df, conf.level = 0.95) #taking out condition and calculate pvalue
correlations_ctl <- NULL
for (i in 1:10) {
  print(i)
  cor_calc <- cor(PCs_ctl, ctl_pca_df[, i], use = "complete.obs",method = "pearson")
  correlations_ctl <- cbind(correlations_ctl,cor_calc)
}

colnames(correlations_ctl) <- c(colnames(predic_ctl))





conflevel.fnf = cor.mtest(fnf_pca_df, conf.level = 0.95) #taking out condition and calculate pvalue
correlations_fnf <- NULL
for (i in 1:10) {
  print(i)
  cor_calc <- cor(Pcs_fnf, fnf_pca_df[, i], use = "complete.obs",method = "pearson")
  correlations_fnf <- cbind(correlations_fnf,cor_calc)
}

colnames(correlations_fnf) <- c(colnames(predict_fnf))

pdf(file = "results_plots/Figure1_differntial_splicing/pearsoncorrelation_fnf.pdf",   # The directory you want to save the file in
    width = 11, # The width of the plot in inches
    height = 6.5) # The height of the plot in inches
# Set up plotting area to handle two plots
par(mfrow = c(1, 2))

# First plot for PBS
corrplot(
  correlations_ctl,
  p.mat = conflevel.ctl$p[grep("PC", colnames(ctl_pca_df)),grep("PC", colnames(ctl_pca_df), invert = TRUE)],
  method = "color",
  tl.col = 'black',
  cl.ratio = 0.2,
  tl.srt = 45,
  insig = 'label_sig',
  sig.level = c(0.001),
  #sig.level = c(0.001, 0.01, 0.05),
  pch.cex = 0.9,
  pch.col = 'grey20',
  addCoef.col = "black",
  number.cex = 0.5,
  col = COL2('BrBG', 10)
)
title("Pearson Correlation - PBS")
#mtext("PBS", side = 1, line = 4, cex = 1.2)

# Second plot for FNF
corrplot(
  correlations_fnf,
  p.mat = conflevel.fnf$p[grep("PC", colnames(fnf_pca_df)),grep("PC", colnames(fnf_pca_df), invert = TRUE)],
  method = "color",
  tl.col = 'black',
  cl.ratio = 0.2,
  tl.srt = 45,
  insig = 'label_sig',
  sig.level = c(0.001),
  #sig.level = c(0.001, 0.01, 0.05),
  pch.cex = 0.9,
  pch.col = 'grey20',
  addCoef.col = "black",
  number.cex = 0.5,
  col = COL2('BrBG', 10)
)
title("Pearson Correlation - FN-f")
#mtext("FNF", side = 1, line = 4, cex = 1.2)
dev.off()
#------------------------------------------------------------------------------
#PBS VS OA
#This is calculating splicing PC for the PBS vs OA
#system("../scripts/Differential_splicing/run_pca.sh /work/users/s/e/seyoun/CQTL_sQTL/output//clu_oa/ctlvsoa_perind.counts.gz 20")

pca_raw_oa <- fread("clu_oa/ctlvsoa_perind.counts.gz.PCs" ) |>
  t() |> 
  as.data.frame()
pca_raw_oa <- pca_raw_oa[-1,]

colnames(pca_raw_oa) <- c(sprintf("PC%s",seq(1:20)))
pca_raw_oa_nm <- cbind(rownames(pca_raw_oa),pca_raw_oa)
colnames(pca_raw_oa_nm)[1] <-c("ID")
pc10_ctl_oa <- pca_raw_oa_nm[,1:11]

#oa_meta_data <- final_meta_data[ ,-c("Race","OAGradeAvg" , "CauseOfDeath","RNAshippedDate","RNAextractionKitBatch","FragmentBatch")]

splicingPCA_oa_df <- merge(final_meta_all, pc10_ctl_oa, by="ID" , all=FALSE)

#Change to factor

# Assuming splicingPCA_df is your dataframe
for (col in colnames(splicingPCA_oa_df)) {
  # Convert only columns starting with "PC" to numeric
  if (!grepl("^PC", col)) {
    splicingPCA_oa_df[[col]] <-  as.factor(splicingPCA_oa_df[[col]])
  }
}

#Seperate out the condition
# Separate the dataframe into two based on the Condition column

colnames(splicingPCA_oa_df)[which(colnames(splicingPCA_oa_df) == 'RNAshippedDate')] <- c("SeqeuncingBatch")
splicingPCA_numeric <- sapply(splicingPCA_oa_df,as.numeric) |> as.data.frame() # Convert columns to numeric, assuming first column is ID

# Subset for CTL condition and drop Condition column
ctl_pca_df <- splicingPCA_numeric %>% subset(Condition == "1") %>% dplyr::select(-Condition)
ctl_pca_df <- ctl_pca_df %>% dplyr::select(-c("Race"))

# Subset for FNF condition and drop Condition column
oa_pca_df <- splicingPCA_numeric %>% subset(Condition == "2") %>% dplyr::select(-Condition)
oc_pca_df_noNA <- oa_pca_df %>% dplyr::select(-c("Race", "OAGradeAvg","CauseOfDeath","FragmentBatch","RNAextractionKitBatch","SeqeuncingBatch"))

# Separate the non-PC columns (predictors) and the PC columns (responses)
predic_ctl <- ctl_pca_df[, grep("PC", colnames(ctl_pca_df), invert = TRUE)]
PCs_ctl <- ctl_pca_df[, grep("PC", colnames(ctl_pca_df))]

predict_oa <- oc_pca_df_noNA[,  grep("PC", colnames(oc_pca_df_noNA), invert = TRUE)]
Pcs_oa <- oc_pca_df_noNA[, grep("PC", colnames(oc_pca_df_noNA))]

conflevel.ctl = cor.mtest(ctl_pca_df, conf.level = 0.95) #taking out condition and calculate pvalue
correlations_ctl <- NULL
for (i in 1:10) {
  print(i)
  cor_calc <- cor(PCs_ctl, ctl_pca_df[, i], use = "complete.obs",method = "pearson")
  correlations_ctl <- cbind(correlations_ctl,cor_calc)
}

colnames(correlations_ctl) <- c(colnames(predic_ctl))


conflevel.oa = cor.mtest(oc_pca_df_noNA, conf.level = 0.95) #taking out condition and calculate pvalue
correlations_oa <- NULL
for (i in 1:5) {
  print(i)
  cor_calc <- cor(Pcs_oa, oc_pca_df_noNA[, i], use = "complete.obs",method = "pearson")
  correlations_oa <- cbind(correlations_oa,cor_calc)
}

colnames(correlations_oa) <- c(colnames(predict_oa))

pdf(file = "results_plots/Figure1_differntial_splicing/pearsoncorrelation_OA.pdf",   # The directory you want to save the file in
    width = 11, # The width of the plot in inches
    height = 7) # The height of the plot in inches
# Set up plotting area to handle two plots
# First plot for PBS
# Set graphics parameters for the desired aspect ratio and margins
layout(matrix(c(1, 2), nrow = 1), widths = c(2.3, 1))
corrplot(
  correlations_ctl,
  p.mat = conflevel.ctl$p[grep("PC", colnames(ctl_pca_df)), grep("PC", colnames(ctl_pca_df), invert = TRUE)],
  method = "color",
  tl.col = 'black',
  cl.ratio = 0.2,
  tl.srt = 45,
  insig = 'label_sig',
  sig.level = c(0.001),
  pch.cex = 0.9,
  pch.col = 'grey20',
  addCoef.col = "black",
  number.cex = 0.5,
  tl.cex = 0.6,
  col = COL2('BrBG', 10),
  bg="transparent",
  mar=c(0,0,2,0)
)
title(main = "Pearson Correlation - CTL",outer = TRUE, line = -1.5, adj = 0.3)

corrplot(
  correlations_oa,
  #p.mat = conflevel.oa$p[grep("PC", colnames(oc_pca_df_noNA)), grep("PC", colnames(oc_pca_df_noNA), invert = TRUE)],
  method = "color",
  tl.col = 'black',
  cl.ratio = 0.2,
  tl.srt = 45,
  insig = 'label_sig',
  sig.level = c(0.001),
  pch.cex = 0.9,
  pch.col = 'grey20',
  addCoef.col = "black",
  number.cex = 0.5,
  tl.cex = 0.6,
  col = COL2('RdBu', 10),
  bg = "transparent",
  mar=c(0,0,5,0)
)
title(main = "Pearson Correlation - OA", outer = TRUE, line = -1.5, adj = 0.95)
dev.off()

# pbs vs OA

splicingPCA_numeric_noNA <- splicingPCA_numeric %>% dplyr::select(-c("Race", "OAGradeAvg","CauseOfDeath","Condition","FragmentBatch"))
predic_all <- splicingPCA_numeric_noNA[, grep("PC", colnames(splicingPCA_numeric_noNA), invert = TRUE)]
PCs_all <- splicingPCA_numeric_noNA[, grep("PC", colnames(splicingPCA_numeric_noNA))]


conflevel.all = cor.mtest(splicingPCA_numeric_noNA, conf.level = 0.99) #taking out condition and calculate pvalue
correlations_all <- NULL
for (i in 1:7) {
  print(i)
  cor_calc <- cor(PCs_all, predic_all[, i], use = "complete.obs",method = "pearson")
  correlations_all <- cbind(correlations_all,cor_calc)
}

colnames(correlations_all) <- c(colnames(predic_all))


corrplot(
  correlations_all,
  p.mat = conflevel.all$p[grep("PC", colnames(splicingPCA_numeric_noNA)),grep("PC", colnames(splicingPCA_numeric_noNA), invert = TRUE)],
  method = "color",
  tl.col = 'black',
  cl.ratio = 0.2,
  tl.srt = 45,
  insig = 'label_sig',
  sig.level = c(0.001),
  #sig.level = c(0.001, 0.01, 0.05),
  pch.cex = 0.9,
  pch.col = 'grey20',
  addCoef.col = "black",
  number.cex = 0.5,
  col = COL2('RdBu', 10)
)



#Wasp---------------------------------------------------------------------------
#pca for ctl vs fnf (WASP)
pca_raw_wasp <- fread("./clu_fnf_wasp/ctlvsfnf_perind.counts.gz.PCs") |> 
  t() |> 
  as.data.frame()
pca_raw_wasp <- pca_raw_wasp[-1,]
colnames(pca_raw_wasp) <- c(sprintf("PC%s",seq(1:20)))
pca_raw_nm <- cbind(rownames(pca_raw_wasp),pca_raw_wasp)
colnames(pca_raw_nm)[1] <-c("ID_wasp")
pc10_ctl_fnf_wasp <- pca_raw_nm[,1:11]

#Merge with the 
splicingPCA_wasp_df <- merge(final_meta_data_wasp, pc10_ctl_fnf_wasp, by="ID_wasp" , all=FALSE)

#Change to factor

# Assuming splicingPCA_df is your dataframe
for (col in colnames(splicingPCA_wasp_df)) {
  # Convert only columns starting with "PC" to numeric
  if (!grepl("^PC", col)) {
    splicingPCA_wasp_df[[col]] <-  as.factor(splicingPCA_wasp_df[[col]])
  }
}


#Seperate out the condition
# Separate the dataframe into two based on the Condition column

colnames(splicingPCA_wasp_df)[which(colnames(splicingPCA_wasp_df) == 'RNAshippedDate')] <- c("SeqeuncingBatch")
splicingPCA_numeric <- sapply(splicingPCA_wasp_df,as.numeric) |> as.data.frame() # Convert columns to numeric, assuming first column is ID

# Subset for CTL condition and drop Condition column
ctl_pca_df <- splicingPCA_numeric %>% subset(Condition == "1") %>% select(-Condition)

# Subset for FNF condition and drop Condition column
fnf_pca_df <- splicingPCA_numeric %>% subset(Condition == "2") %>% select(-Condition)

# Separate the non-PC columns (predictors) and the PC columns (responses)
predic_ctl <- ctl_pca_df[, grep("PC", colnames(ctl_pca_df), invert = TRUE)]
PCs_ctl <- ctl_pca_df[, grep("PC", colnames(ctl_pca_df))]

predict_fnf <- fnf_pca_df[,  grep("PC", colnames(fnf_pca_df), invert = TRUE)]
Pcs_fnf <- fnf_pca_df[, grep("PC", colnames(fnf_pca_df))]

conflevel.ctl = cor.mtest(ctl_pca_df, conf.level = 0.95) #taking out condition and calculate pvalue
correlations_ctl <- NULL
for (i in 1:10) {
  print(i)
  cor_calc <- cor(PCs_ctl, ctl_pca_df[, i], use = "complete.obs",method = "pearson")
  correlations_ctl <- cbind(correlations_ctl,cor_calc)
}

colnames(correlations_ctl) <- c(colnames(predic_ctl))





conflevel.fnf = cor.mtest(fnf_pca_df, conf.level = 0.95) #taking out condition and calculate pvalue
correlations_fnf <- NULL
for (i in 1:10) {
  print(i)
  cor_calc <- cor(Pcs_fnf, fnf_pca_df[, i], use = "complete.obs",method = "pearson")
  correlations_fnf <- cbind(correlations_fnf,cor_calc)
}

colnames(correlations_fnf) <- c(colnames(predict_fnf))

pdf(file = "results_plots/pearsoncorrelation_WASP.pdf",   # The directory you want to save the file in
    width = 11, # The width of the plot in inches
    height = 6.5) # The height of the plot in inches
# Set up plotting area to handle two plots
par(mfrow = c(1, 2))

# First plot for PBS
corrplot(
  correlations_ctl,
  p.mat = conflevel.ctl$p[grep("PC", colnames(ctl_pca_df)),grep("PC", colnames(ctl_pca_df), invert = TRUE)],
  method = "color",
  tl.col = 'black',
  cl.ratio = 0.2,
  tl.srt = 45,
  insig = 'label_sig',
  sig.level = c(0.001),
  #sig.level = c(0.001, 0.01, 0.05),
  pch.cex = 0.9,
  pch.col = 'grey20',
  addCoef.col = "black",
  number.cex = 0.5,
  col = COL2('BrBG', 10)
)
title("Pearson Correlation - PBS")
#mtext("PBS", side = 1, line = 4, cex = 1.2)

# Second plot for FNF
corrplot(
  correlations_fnf,
  p.mat = conflevel.fnf$p[grep("PC", colnames(fnf_pca_df)),grep("PC", colnames(fnf_pca_df), invert = TRUE)],
  method = "color",
  tl.col = 'black',
  cl.ratio = 0.2,
  tl.srt = 45,
  insig = 'label_sig',
  sig.level = c(0.001),
  #sig.level = c(0.001, 0.01, 0.05),
  pch.cex = 0.9,
  pch.col = 'grey20',
  addCoef.col = "black",
  number.cex = 0.5,
  col = COL2('RdBu', 10)
)
title("Pearson Correlation - FNF")
#mtext("FNF", side = 1, line = 4, cex = 1.2)
dev.off()


#-------------------------------------------------------------------------------

#To check Genotype PCA

#-------------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(gridExtra)

# read in data
pca_pbs <- fread("geno/pbs_geno/05.pca/CQTL.eigenvec")
pca_fnf <- fread("geno/fnf_geno/05.pca/CQTL.eigenvec")
eigenval_pbs <- scan("geno/pbs_geno/05.pca/CQTL.eigenval")
eigenval_fnf <- scan("geno/fnf_geno/05.pca/CQTL.eigenval")

pca_pbs <- pca_pbs[,-1]
# set names
names(pca_pbs)[1] <- "ID"
names(pca_pbs)[2:ncol(pca_pbs)] <- paste0("PC", 1:(ncol(pca_pbs)-1))

pca_fnf <- pca_fnf[,-1]
# set names
names(pca_fnf)[1] <- "ID"
names(pca_fnf)[2:ncol(pca_fnf)] <- paste0("PC", 1:(ncol(pca_fnf)-1))


# first convert to percentage variance explained
#https://speciationgenomics.github.io/pca/
pve_pbs <- data.frame(PC = 1:20, pve = eigenval_pbs/sum(eigenval_pbs)*100)
pve_fnf <- data.frame(PC = 1:20, pve = eigenval_fnf/sum(eigenval_fnf)*100)


geno_pca_pbs <- merge(final_meta_data,pca_pbs, by="ID",all.x=FALSE)
geno_pca_pbs$Race <- as.factor(geno_pca_pbs$Race)

pdf(file = "results_plots/Genotype_pc_eigenvalue.pdf",   # The directory you want to save the file in
    width = 11, # The width of the plot in inches
    height = 6) # The height of the plot in inches

# Assuming pve_pbs and geno_pca_pbs data frames are already prepared

# Plot 1: Variance Explained by PCs (PBS)
p1 <- ggplot(pve_pbs[1:10,], aes(x = PC, y = pve)) + 
  geom_line(color = "lightblue") + 
  geom_point(shape = 2) + 
  ylab("Percentage variance explained") + 
  scale_x_continuous(breaks = 1:10, limits = c(1, 10)) + 
  theme_light() +
  ggtitle("A")

# Plot 2: PCA Plot (PBS)
p2 <- ggplot(geno_pca_pbs, aes(x = PC1, y = PC2, color = Race)) +
  geom_point() +
  coord_equal() +
  theme_minimal() +
  labs(title = "B", x = "PC1", y = "PC2", color = "Race") +
  xlab(paste0("PC1 (", signif(pve_pbs$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve_pbs$pve[2], 3), "%)"))

# Arrange plots side by side
grid.arrange(p1, p2, ncol = 2)

dev.off()



#Prepare the covariates input text

#geno merge
geno_pcs <-  rbind(pca_pbs,pca_fnf)
#splice_pcs <- fread("./clu_fnf/ctlvsfnf_perind.counts.gz.PCs") |> t() |> as.data.frame()

ctl_meta_subset <- splicingPCA_df %>% subset(Condition == "CTL") %>% select(ID ,RNAextractionKitBatch, FragmentBatch ,SeqeuncingBatch)
fnf_mete_subset <- splicingPCA_df %>% subset(Condition == "FNF") %>% select(ID ,RNAextractionKitBatch, FragmentBatch ,SeqeuncingBatch)


for (i in 1:length(splice_pcs)) {
  print(i)
  splice_pc_subset <- splice_pcs[,1:i] |> as.data.frame()
  col_n <-ncol(splice_pc_subset)
  colnames(splice_pc_subset) <-  c(sprintf("SplicePC%s",seq(1:col_n)))
  splice_id_df <- cbind(rownames(splice_pcs),splice_pc_subset)
  colnames(splice_id_df)[1] <- c("ID")
  splice_id_df <- splice_id_df[-1,]
  spliceMeta_df <- merge(splice_id_df ,ctl_meta_subset, by="ID" ,all.x=F)
  cov_df <- merge(geno_pcs[,1:11],spliceMeta_df, by="ID" ,all.x=F)
  cov_df_num <- sapply(cov_df,as.numeric) |> as.data.frame()
  cov_df_num$ID <- c(cov_df$ID)
  cov_df_t <- cov_df_num |> t()
  colnames(cov_df_t) <- cov_df_t[1, ]
  cov_df_t <- cov_df_t[-1, ]
  cov_fi <- cbind(rownames(cov_df_t),cov_df_t)
  colnames(cov_fi)[1] <- c("id")
  write.table(cov_fi,file=paste0("/work/users/s/e/seyoun/CQTL_sQTL/output/clu_fnf/Covariates/","PBS_covariates_PC",col_n),sep='\t',quote=F,row.names=F,col.names=T)
}


for (i in 1:length(splice_pcs)) {
  print(i)
  splice_pc_subset <- splice_pcs[,1:i]|> as.data.frame()
  col_n <-ncol(splice_pc_subset)
  colnames(splice_pc_subset) <-  c(sprintf("SplicePC%s",seq(1:col_n)))
  splice_id_df <- cbind(rownames(splice_pc_subset),splice_pc_subset)
  colnames(splice_id_df)[1] <- c("ID")
  splice_id_df <- splice_id_df[-1,]
  spliceMeta_df <- merge(splice_id_df ,fnf_mete_subset, by="ID" ,all.x=F)
  cov_df <- merge(geno_pcs[,1:11],spliceMeta_df, by="ID" ,all.x=F)
  cov_df_num <- sapply(cov_df,as.numeric) |> as.data.frame()
  cov_df_num$ID <- c(cov_df$ID)
  cov_df_t <- cov_df_num |> t()
  colnames(cov_df_t) <- cov_df_t[1, ]
  cov_df_t <- cov_df_t[-1, ]
  cov_fi <- cbind(rownames(cov_df_t),cov_df_t)
  colnames(cov_fi)[1] <- c("id")
  write.table(cov_fi,file=paste0("/work/users/s/e/seyoun/CQTL_sQTL/output/clu_fnf/Covariates/","FNF_covariates_PC",col_n),sep='\t',quote=F,row.names=F,col.names=T)
}


#------------------------------------------------------------------
#qqnorm file changes for qtltools
#------------------------------------------------------------------

li_qqnorm <- list.files(path ="/work/users/s/e/seyoun/CQTL_sQTL/output/clu_fnf/",
                 pattern = "*qqnorm_chr([1-9]|1[0-9]|2[0-2]|X|Y)$",
                 full.names = TRUE)


for (i in 1:length(li_qqnorm)){
  print(i)
  df <- fread(li_qqnorm[i])
  sp_number <- gsub(".gz","",strsplit(strsplit(li_qqnorm[i],"/")[[1]][13],"_")[[1]][3])
  df$strand <- do.call(rbind,strsplit(df$ID, "_"))[,3]
  df$pid <-do.call(rbind,strsplit(df$ID, ":"))[,4]
  new_dt <- df[,c(1:4,212,211,5:210)]
  write.table(new_dt,file=paste0("ctrvsfnf_qqnorm_",sp_number,".bed"),sep='\t',quote=F,row.names=F,col.names=T)
}
  



            
