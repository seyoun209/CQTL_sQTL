# Make a covariate plot to makre sure I am choosing the right covariates
setwd("/work/users/s/e/seyoun/CQTL_sQTL")
library(ggplot2)
library(corrplot)

#For the corrplot star as significant do this(https://stackoverflow.com/questions/63227830/r-corrplot-plot-correlation-coefficients-along-with-significance-stars_:
#trace(corrplot, edit=TRUE)
# and the find the line 443
# adjust text(X,Y ...) according to your needs, here +0.25 is added to the Y-position    
#place_points = function(sig.locs, point) {
#  text(pos.pNew[, 1][sig.locs], (pos.pNew[, 2][sig.locs])+0.25, 
#       labels = point, col = pch.col, cex = pch.cex, 
#       lwd = 2)
# and then Save
#-------------------------------------------------------------------------------

load("output/combined_meta_data.RData")
# Select and rename columns dynamically, based on the use of wasp ID
selected_columns <- c("ID","Condition","Sex", "Age","Race","OAGradeAvg","CauseOfDeath","FragmentBatch",
                      "RIN","RNAextractionKitBatch","RNAshippedDate")

final_meta_all <- combined_data %>% dplyr::select(all_of(selected_columns))


#------------------------------------------------------------------------------
#PBS vs OA
pc_dir  <- "/work/users/s/e/seyoun/CQTL_sQTL/output/gtex_cluster_oa/pbs_oa.leafcutter.PCs.txt"
# PCA for PBS vs OA
pca_raw <- fread(pc_dir) |>
  t() |>
  as.data.frame()
pca_raw <- pca_raw[-1,]
colnames(pca_raw) <- c(sprintf("PC%s",seq(1:20)))
pca_raw_nm <- cbind(rownames(pca_raw),pca_raw)
colnames(pca_raw_nm)[1] <-c("ID")
pc10_pbs_oa <- pca_raw_nm[,1:11]

splicingPCA_oa_df <- merge(final_meta_all, pc10_pbs_oa, by="ID" , all=FALSE)

#Chage to a factor

# Assuming splicingPCA_df is your dataframe
for (col in colnames(splicingPCA_oa_df)) {
  # Convert only columns starting with "PC" to numeric
  if (!grepl("^PC", col)) {
    splicingPCA_oa_df[[col]] <-  as.factor(splicingPCA_oa_df[[col]])
  }
}

colnames(splicingPCA_oa_df)[which(colnames(splicingPCA_oa_df) == 'RNAshippedDate')] <- c("SeqeuncingBatch")
splicingPCA_numeric <- sapply(splicingPCA_oa_df,as.numeric) |> as.data.frame() # Convert columns to numeric, assuming first column is ID

# Subset for CTL condition and drop Condition column
ctl_pca_df <- splicingPCA_numeric %>% subset(Condition == "1") %>% dplyr::select(-Condition)
ctl_pca_df <- ctl_pca_df %>% dplyr::select(-c("Race"))

# Subset for FNF condition and drop Condition column
oa_pca_df <- splicingPCA_numeric %>% subset(Condition == "2") %>% dplyr::select(-Condition)
oa_pca_df_noNA <- oa_pca_df %>% dplyr::select(-c("Race", "OAGradeAvg","CauseOfDeath","FragmentBatch","RNAextractionKitBatch","SeqeuncingBatch"))


# Separate the non-PC columns (predictors) and the PC columns (responses)
predic_ctl <- ctl_pca_df[, grep("PC", colnames(ctl_pca_df), invert = TRUE)]
PCs_ctl <- ctl_pca_df[, grep("PC", colnames(ctl_pca_df))]

predict_oa <- oa_pca_df_noNA[,  grep("PC", colnames(oa_pca_df_noNA), invert = TRUE)]
Pcs_oa <- oa_pca_df_noNA[, grep("PC", colnames(oa_pca_df_noNA))]

conflevel.ctl = cor.mtest(ctl_pca_df, conf.level = 0.95) #taking out condition and calculate pvalue
correlations_ctl <- NULL
for (i in 1:9) {
  print(i)
  cor_calc <- cor(PCs_ctl, ctl_pca_df[, i], use = "complete.obs",method = "pearson")
  correlations_ctl <- cbind(correlations_ctl,cor_calc)
}

colnames(correlations_ctl) <- c(colnames(predic_ctl))


conflevel.oa = cor.mtest(oa_pca_df_noNA, conf.level = 0.95) #taking out condition and calculate pvalue
correlations_oa <- NULL
for (i in 1:4) {
  print(i)
  cor_calc <- cor(Pcs_oa, oa_pca_df_noNA[, i], use = "complete.obs",method = "pearson")
  correlations_oa <- cbind(correlations_oa,cor_calc)
}

colnames(correlations_oa) <- c(colnames(predict_oa))


#-------------------------------------------------------------------------------
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
  p.mat = conflevel.oa$p[grep("PC", colnames(oa_pca_df_noNA)), grep("PC", colnames(oa_pca_df_noNA), invert = TRUE)],
  method = "color",
  tl.col = 'black',
  cl.ratio = 0.2,
  tl.srt = 45,
  insig = 'label_sig',
  sig.level = c(0.005),
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


