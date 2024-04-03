# Making OA donors compared plot with the current CTL vs FNF 
## Author: Seyoun Byun
## Date: 03.28.2024
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

pbs_oa <- load("output/clu_oa/PBSvsOA.Rdata")

ratios_oa <- counts %>%
  mutate(clu = str_split_fixed(rownames(counts), ":", 4)[,4]) %>%
  group_by(clu) %>%
  mutate(across(everything(), ~./sum(.))) %>%
  ungroup() %>%
  as.data.frame() %>%
  set_rownames(rownames(counts)) %>%
  dplyr::select(-clu)

ratios_oa = ratios_oa[rowMeans(is.na(ratios_oa)) <= 0.4,,drop=F ] #Try to remove 40% or less NA values in each rows
#From 134376 to 133944 --> it is dropped 432
row_means = rowMeans(ratios_oa, na.rm = T) #calculate the mean without NAs in the rows. 
row_means_outer = outer(row_means, rep(1,ncol(ratios_oa)))# making outlier to the NAs
ratios_oa[is.na(ratios_oa)] = row_means_outer[is.na(ratios_oa)] # instead of Na, add the rowmean
colnames(ratios_oa) <- colnames(counts_oa)
ratios_oa <- cbind(rownames(ratios_oa), ratios_oa)
colnames(ratios_oa)[1] <- c('Junction')
write.table(ratios_oa, file = "output/clu_oa/ratio_oa.txt",sep='\t',quote=F,row.names=F,col.names=T)



#Finding the Differential gene expression and plot (Cluster based)
# Finding the one introns per cluster that is the most different.


for (name in pbs_oa) {
  # Construct the new name by prefixing with "fnf_"
  new_name <- paste0(name,"_oa")
  
  # Assign the object to the new name in the global environment
  assign(new_name, get(name))
  remove(name)
}

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





#-------------------------------------------------------------------------------
#reading the meta and find the significant introns

meta_cqtl <- fread("output/clu_fnf/meta_cqtl")
# Omit samples and filter rows based on configuration and conditions
meta_ctl_fnf <- meta_cqtl[!meta_cqtl$Donor %in% config$samples_to_omit, ] %>%
  filter(Condition %in% conditions_to_include)
psi_fnf <- read.table("output/clu_fnf/ratio_fnf.txt",sep="\t",header=T)
rownames(psi_fnf) <- psi_fnf[,1]
psi_fnf <- psi_fnf[, -1]
model_fnf <-model.matrix(~Condition,data=meta_ctl_fnf)
model_cov_fnf <-model.matrix(~Donor,data=meta_ctl_fnf)
limma_psi_batchremove <- limma::removeBatchEffect(psi_fnf[, -ncol(psi_fnf)], batch = factor(meta_ctl_fnf$RNAextractionKitBatch),
                                                  batch2 = factor(meta_ctl_fnf$FragmentBatch),
                                                  covariates = as.numeric(factor(meta_ctl_fnf$Donor)))
#limma_psi_batchremove.df <- limma_psi_batchremove |> as.data.frame()
fit_fnf <- lmFit(limma_psi_batchremove,design = model_fnf)
fiteBays_fnf<- eBayes(fit_fnf)
table_fnf <- topTable(fiteBays_fnf, coef = "ConditionFNF",n=Inf,adjust.method='fdr')
sig_fnf_lima <- table_fnf[table_fnf$adj.P.Val <=0.05 ,]
sig_fnf_lima <- sig_fnf_lima %>%
  mutate(chr = str_split(rownames(sig_fnf_lima), ":", simplify = TRUE)[, 1],
         start = as.double(str_split(rownames(sig_fnf_lima), ":", simplify = TRUE)[, 2]),
         end = as.double(str_split(rownames(sig_fnf_lima), ":", simplify = TRUE)[, 3]),
         clusterID = str_split(rownames(sig_fnf_lima), ":", simplify = TRUE)[, 4])

merged_sig <- merge(sig_fnf_lima, introns, by ="row.names",all.x=TRUE)
sig_fnf_lima_fc0.2 <- sig_fnf_lima[(abs(sig_fnf_lima$logFC) >= 0.2),]


limma_psi_batchremove <- limma::removeBatchEffect(counts_fnf, batch = factor(meta_ctl_fnf$RNAextractionKitBatch),
                                                  batch2 = factor(meta_ctl_fnf$FragmentBatch),
                                                  design =model_fnf,
                                                  covariates = as.numeric(factor(meta_ctl_fnf$Donor)))
limma_psi_batchremove.df <- limma_psi_batchremove |> as.data.frame()
limma_psi_batchremove.df$deltaPSI <- apply(limma_psi_batchremove.df, 1, function(row) {
  ctl_cols <- grep("CTL", colnames(limma_psi_batchremove.df))
  fnf_cols <- grep("FNF", colnames(limma_psi_batchremove.df))
  
  ctl_mean <- mean(as.numeric(row[ctl_cols]), na.rm = TRUE)
  fnf_mean <- mean(as.numeric(row[fnf_cols]), na.rm = TRUE)
  
  if (is.na(ctl_mean) || is.na(fnf_mean)) {
    NA
  } else {
    fnf_mean - ctl_mean
  }
})



psi_fnf$location <- sapply(strsplit(psi_fnf$Junction, ":"), function(x) paste0(x[1:3], collapse = ":"))
#Find the deltapsi 20% differenfce in FNF
introns_fnf_deltapsi2 <-introns_fnf[abs(introns_fnf$deltapsi) >=  0.2,] 
introns_fnf_deltapsi2$location <- paste0(introns_fnf_deltapsi2$chr, ":", introns_fnf_deltapsi2$start, ":", introns_fnf_deltapsi2$end)
introns_fnf_deltapsi2$Junction <- paste0(introns_fnf_deltapsi2$chr, ":", introns_fnf_deltapsi2$start, ":", introns_fnf_deltapsi2$end, ":", introns_fnf_deltapsi2$clusterID)


#Find the deltapsi 20% differenfce in OA
psi_oa <- fread("output/clu_oa/ratio_oa.txt",sep="\t",header=T)
psi_oa$location <- sapply(strsplit(psi_oa$Junction, ":"), function(x) paste0(x[1:3], collapse = ":"))
introns_oa_deltapsi2 <-introns_oa[abs(introns_oa$deltapsi) >=  0.2,] 

# Create the location column
introns_oa_deltapsi2$location <- paste0(introns_oa_deltapsi2$chr, ":", introns_oa_deltapsi2$start, ":", introns_oa_deltapsi2$end)

# Create the Junction column
introns_oa_deltapsi2$Junction <- paste0(introns_oa_deltapsi2$chr, ":", introns_oa_deltapsi2$start, ":", introns_oa_deltapsi2$end, ":", introns_oa_deltapsi2$clusterID)

up_introns_oa_deltapsi2 <- introns_oa_deltapsi2[introns_oa_deltapsi2$deltapsi > 0,]
down_introns_oa_deltapsi2  <- introns_oa_deltapsi2[introns_oa_deltapsi2$deltapsi < 0,]


psi_fnf$deltaPSI <- apply(psi_fnf, 1, function(row) {
  ctl_cols <- grep("CTL", colnames(psi_fnf))
  fnf_cols <- grep("FNF", colnames(psi_fnf))
  
  ctl_mean <- mean(as.numeric(row[ctl_cols]), na.rm = TRUE)
  fnf_mean <- mean(as.numeric(row[fnf_cols]), na.rm = TRUE)
  
  if (is.na(ctl_mean) || is.na(fnf_mean)) {
    NA
  } else {
    fnf_mean - ctl_mean
  }
})

fnf_up_oa_subset <- psi_fnf |>
  filter(location %in% up_introns_oa_deltapsi2$location )

fnf_down_oa_subset <- psi_fnf |>
  filter(location %in% down_introns_oa_deltapsi2$location )


fnf_all_oa_subset <- bind_rows(fnf_up_oa_subset |> mutate(group = "Up in OA"),
                               fnf_down_oa_subset |> mutate(group = "Down in OA")) |>
  mutate(group = factor(group, levels = c("Up in OA", "Down in OA")))







up_test_oa <- wilcox.test(x = fnf_up_oa_subset$deltaPSI,
                          y = psi_fnf |>
                            filter(! location %in% fnf_up_oa_subset$location) |>
                            pull(deltaPSI))


down_test_oa <- wilcox.test(x = fnf_down_oa_subset$deltaPSI,
                          y = psi_fnf |>
                            filter(! location %in% fnf_down_oa_subset$location) |>
                            pull(deltaPSI), alternative = "less")

