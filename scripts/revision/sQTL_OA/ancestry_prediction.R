# Ancestry for OA sample
library(data.table)
library(dplyr)
library(readr)
library(ggplot2)
library(e1071)
setwd("/proj/phanstiel_lab/Data/processed/CQTL/geno/COA8_OA")

## Read in pca data
pcaData <- fread("ancestry/CQTL_COA8_ref_merge.pca.evec", data.table = FALSE)

## Grab sample names and first 2 PC's
pcaData <- pcaData[,c(1:3)]
colnames(pcaData) <- c("sample", "PC1","PC2")

## Read in panel data
panel <- fread("/proj/phanstiel_lab/References/genomes/1000G/GRCh38/igsr-1000 genomes on grch38.tsv", data.table = FALSE) |>
  dplyr::select(sample ="Sample name",
         pop ="Population code",
         super_pop ="Superpopulation code",
         gender ="Sex")

#write.table(panel,file="/proj/phanstiel_lab/References/genomes/1000G/GRCh38/1000G_phase3.Updatedpanel", row.names = F, col.names = T,
#            sep = ",", quote=F)
## Match ref panel to pca data based on "sample" column

pca_panel <- left_join(pcaData, panel, by = c("sample"))


## SVN to infer super_pop
pca_panel_train <- pca_panel %>% filter(!is.na(super_pop))
pca_panel_train$super_pop <- as.factor(pca_panel_train$super_pop)
pca_panel_test <- pca_panel %>% filter(is.na(super_pop))
svm_ancestry <- svm(super_pop~PC1+PC2, data = pca_panel_train,
                    type = "C-classification", kernel = "radial")
prediction <- predict(svm_ancestry, pca_panel_test[,c("PC1", "PC2")])
pca_panel_test$super_pop <- prediction

pca_panel_test %>% dplyr::select(sample, super_pop) %>%
  dplyr::rename("Donor" = sample) %>%
  dplyr::rename("Predicted_Ancestry" = super_pop) %>%
  write_csv(file = "Current_study_predictedAncestry.csv")
