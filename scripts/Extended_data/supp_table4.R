# Supp table4 -  after edits of SNRNP70 alt.HET-KD ex8, differentially expressed genes and exon level counts (?)
setwd("/work/users/s/e/seyoun/CQTL_sQTL/")
library(data.table)
library(openxlsx)
library(dplyr)
library(tidyr)

#Function
write_formatted_sheet <- function(wb, sheet_name, data) {
  # Add worksheet
  addWorksheet(wb, sheet_name)
  
  # Write data
  writeData(wb, sheet_name, data, startRow = 1, startCol = 1)
  
  # Create styles
  headerStyle <- createStyle(textDecoration = "bold", halign = "center", valign = "center", fgFill = "grey90")
  dataStyle <- createStyle(wrapText = TRUE, halign = "left", valign = "center")
  
  
  # Apply header style
  addStyle(wb, sheet_name, headerStyle, rows = 1, cols = 1:ncol(data), gridExpand = TRUE)
  
  # Apply data style to all cells except header
  addStyle(wb, sheet_name, dataStyle, rows = 2:(nrow(data)+1), cols = 1:ncol(data), gridExpand = TRUE)
  
  
  # Set column widths
  setColWidths(wb, sheet_name, cols = 1:ncol(data), widths = "auto")
  
  # Freeze the top row
  freezePane(wb, sheet_name, firstRow = TRUE)
}
# First, getting the exon counts.

# Function to process exon tables
process_exon_table <- function(exon_table, chr_pos_data) {
  exon_table %>%
    mutate(ExonID = gsub("SNRNP70:", "", GeneExon)) %>%
    select(ExonID, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
    left_join(chr_pos_data, by = c("ExonID" = "exonid")) %>%
    select(ExonID, Chr, Start, End, baseMean, log2FoldChange, padj) %>%
    arrange(Start)
}

# Chromosome and position data
chr_pos_data <- data.frame(
  exonid = c("ENSE00001555728.1", "ENSE00003535913.1", "ENSE00003567241.1", 
             "ENSE00003499363.1", "ENSE00003552141.1", "ENSE00003683920.1", 
             "ENSE00002503116.1", "ENSE00001549580.1", "ENSE00003593260.1", 
             "ENSE00003614741.1", "ENSE00003179801.2"),
  Chr = rep("chr19", 11),
  Start = c(49085451, 49086405, 49090291, 49090466, 49098427, 49098642, 
            49101390, 49102114, 49104634, 49107625, 49107795),
  End = c(49085451, 49086405, 49090291, 49090466, 49098427, 49098642, 
          49101390, 49102114, 49104634, 49107625, 49107795)
)

load("output/featurecounts_snrnp70/exon.table_oa")
load("output/featurecounts_snrnp70/exon.table.fnf")
load("output/featurecounts_snrnp70/exon.table_edited")

# Process each table
fnf_processed <- process_exon_table(exon.table_fnf, chr_pos_data)
oa_processed <- process_exon_table(exon.table_oa, chr_pos_data)
wd_hetkd_processed <- process_exon_table(exon.table, chr_pos_data)



#-------------------------------------------------------------------------------
# Differential gene expression edited WD vs. HET-KD

# Read the down and up regulated genes (padj 0.01 and log2fodchange 2)
downsig_deGenes_pval01_l2fc2 <- fread("/work/users/s/e/seyoun/crispr/02.test_seq/condition_de/downsig_deGenes_pval01_l2fc2.csv")
upsig_deGenes_pval01_l2fc2 <- fread("/work/users/s/e/seyoun/crispr/02.test_seq/condition_de/upsig_deGenes_pval01_l2fc2.csv")

# Function to process gene tables
process_gene_table <- function(gene_table, regulation) {
  gene_table %>%
    dplyr::rename(
      Gene = SYMBOL,
      `Ensembl ID` = gene_id,
      Chr = seqnames
    ) %>%
    dplyr::select(Gene, `Ensembl ID`, Chr, start, end, log2FoldChange, padj) %>%
    mutate(Direction = regulation) %>%
    arrange(Chr, start)
}

# Process upregulated and downregulated tables
up_processed <- process_gene_table(upsig_deGenes_pval01_l2fc2, "Upregulated")
down_processed <- process_gene_table(downsig_deGenes_pval01_l2fc2, "Downregulated")

# Combine the tables
combined_table <- bind_rows(up_processed, down_processed) %>%
  arrange(Chr, start)


#Pathway and GO term


GO_up_edited <- fread("/work/users/s/e/seyoun/crispr/02.test_seq/condition_de/table/GO_Upsig.csv")
Go_down_edited <- fread("/work/users/s/e/seyoun/crispr/02.test_seq/condition_de/table/GO_Downsig.csv")
Path_up_edited <-  fread("/work/users/s/e/seyoun/crispr/02.test_seq/condition_de/table/pathway_up.csv",drop="category")
Path_down_edited <-  fread("/work/users/s/e/seyoun/crispr/02.test_seq/condition_de/table/pathway_down.csv",drop="category")




# Create a new workbook
wb <- createWorkbook()
write_formatted_sheet(wb, "Diff-Exon PBS vs FN-f", fnf_processed)
write_formatted_sheet(wb, "Diff-Exon PBS vs OA", oa_processed)
write_formatted_sheet(wb, "Diff-Exon WD vs HET-KD", wd_hetkd_processed)
write_formatted_sheet(wb, "Diff-GeneExp", combined_table)
write_formatted_sheet(wb, "Pathway-up diffgenes", Path_up_edited)
write_formatted_sheet(wb, "Pathway-down diffgenes", Path_down_edited)
write_formatted_sheet(wb, "GO-up diffgenes", GO_up_edited)
write_formatted_sheet(wb, "GO-down diffgene", Go_down_edited)
saveWorkbook(wb, "output/results_plots/Supp_Data/Supplementary_Data4.xlsx", overwrite = TRUE)
