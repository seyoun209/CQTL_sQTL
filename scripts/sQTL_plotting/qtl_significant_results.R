# Finding the total counts and find the significant QTL
## Author: Seyoun Byun
## Date: 03.08.2024
## Edited: 06.11.2024
#-------------------------------------------------------------------------------
setwd("/work/users/s/e/seyoun/CQTL_sQTL")
library(tidyr)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(ensembldb)
library(AnnotationHub)
library(AnnotationFilter)
library(org.Hs.eg.db)
library(ggtext)
source("scripts/sQTL_rscripts/utils.R")

# need to run the sh qtltools_merge.sh 


# Define the possible arguments files, significant p-values etc
adjusted_beta_p <- 0.05
pbs_perm_dir <- "output/01.qtltools_re/perm_pbs/"
fnf_perm_dir <- "output/01.qtltools_re/perm_fnf/"
pbs_wasp_perm_dir <- "output/01.qtltools_re/perm_pbs_wasp"
fnf_wasp_perm_dir <- "output/01.qtltools_re/perm_fnf_wasp"


# Function for reading and processing PBS files --------------------------------
process_pbs_files <- function(filepath_pattern, pattern_suffix) {
filepaths <- list.files(filepath_pattern, pattern = pattern_suffix, full.names = TRUE)   # List files based on the pattern

# Sort filepaths numerically based on the PC number
sorted_filepaths <- filepaths[order(as.numeric(gsub(".*pc(\\d+)_.*", "\\1", filepaths)))]

# Read files and process them
list_processed <- lapply(sorted_filepaths, function(filepath) {
  df <- read.table(filepath, header = FALSE, stringsAsFactors = FALSE)
  
# Filter out rows with "chrX" or "chrY" from the specified column
df_filtered <- df[!(df[, 2] %in% c("chrX", "chrY")), ]
    
# Add adjusted p-values using FDR
df_filtered$p_adjusted <- p.adjust(df_filtered[, 20], method = "fdr")
    
# Filter based on significance threshold
df_sig <- df_filtered[df_filtered$p_adjusted <= adjusted_beta_p, ]

return(df_sig)
})
  
  return(list_processed)
}


pbs_list_processed <- process_pbs_files(pbs_perm_dir, "_allchr.pbs.perm")
fnf_list_processed <- process_pbs_files(fnf_perm_dir, "_allchr.fnf.perm")

pbs_wasp_list_processed <- process_pbs_files(pbs_wasp_perm_dir, "_allchr.pbs.wasp.perm")
fnf_wasp_list_processed <- process_pbs_files(fnf_wasp_perm_dir, "_allchr.fnf.wasp.perm")



#dotplot to see the counts of snps only for the ---------------------------------

# Count unique values in column 1 of each list element
pbs_unique_counts <- lapply(pbs_list_processed, function(x) {
  length(unique(x$V1))-1
})

fnf_unique_counts <- lapply(fnf_list_processed, function(x) {
  length(unique(x$V1))-1
})


# Create a data frame with the list elements and names
df <- data.frame(pc = paste0("pc", 1:20), pbs_counts = unlist(pbs_unique_counts), fnf_counts=unlist(fnf_unique_counts))
df$pc <- factor(df$pc, levels = paste0("pc", 1:20))
df_long <- df %>%
  pivot_longer(cols = c(pbs_counts, fnf_counts),
               names_to = "type",
               values_to = "counts") %>%
  mutate(type = dplyr::recode(type, pbs_counts = "PBS", fnf_counts = "FNF"))

df_long$type <- factor(df_long$type, levels = c("PBS","FNF"))

# Calculate the max for each group
max_counts <- df_long %>%
  dplyr::group_by(type) %>%
  dplyr::summarize(max_count = max(counts), pc_max = pc[which.max(counts)]) %>%
  ungroup()

max_counts$color <- ifelse(max_counts$type == "PBS", "#18417c", "#F2AA40")


#pdf(file = "output/results_plots/Significant_count_introns.pdf",   # The directory you want to save the file in
#    width = 11, # The width of the plot in inches
#    height = 6.5) # The height of the plot in inches

# Create the dot plot
counts_dotPlot <- ggplot(df_long, aes(x = pc, y = counts, color = type)) +
  geom_point(size = 2) +
  labs(x = "PC", y = "Number of significant introns", color = "Type") +
  scale_color_manual(values = c("PBS" = "#2057A7", "FNF" = "#F2BC40")) +
  scale_y_continuous(limits = c(2500, 6000)) +
  theme_classic() +
  # geom_segment(data = max_counts, aes(x = "pc6", xend = pc_max, y = max_count +150, yend = max_count), 
  #              arrow = arrow(type = "closed", length = unit(0.2, "cm")), color = "black") +
  geom_text(data = max_counts, aes(x = pc_max, y = max_count + 100, 
                                   label = paste0(pc_max, ": ", max_count)),
            color = max_counts$color, vjust = 0, size = 3.5)+
  theme(text = element_text(family = "Helvetica"),
        panel.background = element_rect(fill = 'transparent', color = "transparent"),
        plot.background = element_rect(fill = 'transparent', color = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.background = element_blank(),
        legend.position=c(0.9,0.9),
        axis.title = element_text(size = 9),
        axis.line = element_line(color = "black", 
                                 linewidth = 0.25),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.ticks.y = element_line(color = "black",
                                    linewidth = 0.25),
        axis.text = element_text(color = "black", size = 8),
        axis.ticks.x = element_line(color = "black",
                                    linewidth = 0.25),
        axis.text.x = element_text(color = "black", size = 8,angle=45),
        plot.margin = margin(t = 5, r = 10, b = 5, l = 5))
#dev.off()
save(counts_dotPlot, file = "output/results_plots/sqtl_plots/counts_dotplot.rda")
# #dotplot to see WASP counts------------------- ---------------------------------
# 
# # Count unique values in column 1 of each list element
# pbs_wasp_counts <- lapply(pbs_wasp_list_processed, function(x) {
#   length(unique(x$V1))
# })
# 
# fnf_wasp_counts <- lapply(fnf_wasp_list_processed, function(x) {
#   length(unique(x$V1))
# })
# 
# 
# # Create a data frame with the list elements and names
# df_wasp <- data.frame(pc = paste0("pc", 1:20), pbs_counts = unlist(pbs_wasp_counts), fnf_counts=unlist(fnf_wasp_counts))
# df_wasp$pc <- factor(df_wasp$pc, levels = paste0("pc", 1:20))
# df_wasp_long <- df_wasp %>%
#   pivot_longer(cols = c(pbs_counts, fnf_counts),
#                names_to = "type",
#                values_to = "counts") %>%
#   mutate(type = recode(type, pbs_counts = "PBS", fnf_counts = "FNF"))
# 
# df_wasp_long$type <- factor(df_wasp_long$type, levels = c("PBS","FNF"))
# 
# # Calculate the max for each group
# max_counts_wasp <-df_wasp_long %>%
#   group_by(type) %>%
#   summarize(max_count = max(counts), pc_max = pc[which.max(counts)]) %>%
#   ungroup()
# 
# max_counts_wasp$color <- ifelse(max_counts_wasp$type == "PBS", "#2ea7e8", "#f77943")
# 
# pdf(file = "output/results_plots/wasp_Significant_count_introns.pdf",   # The directory you want to save the file in
#     width = 11, # The width of the plot in inches
#     height = 6.5)
# # Create the dot plot
# ggplot(df_wasp_long, aes(x = pc, y = counts, color = type)) +
#   geom_point(size = 3) +
#   labs(x = "PC", y = "Number of significant introns", color = "Type") +
#   scale_color_manual(values = c("PBS" = "#9FCCE4", "FNF" = "#FAB394")) 
#   scale_y_continuous(limits = c(300, 700)) +
#   theme_classic() +
#   # geom_segment(data = max_counts, aes(x = "pc6", xend = pc_max, y = max_count +150, yend = max_count), 
#   #              arrow = arrow(type = "closed", length = unit(0.2, "cm")), color = "black") +
#   geom_text(data = max_counts_wasp, aes(x = pc_max, y = max_count + 10, label = max_count), 
#             color = max_counts_wasp$color, vjust = 0, size = 3.5)
# 
# 
# dev.off()


#check pvalue are okay ---------------------------------------------------------
perm_pbs_pc5 <- read.table(paste0(pbs_perm_dir,"pc5_allchr.pbs.perm"), header = FALSE, stringsAsFactors = FALSE)
perm_fnf_pc4 <- read.table(paste0(fnf_perm_dir,"pc4_allchr.fnf.perm"), header = FALSE, stringsAsFactors = FALSE)
header <- read.table("scripts/sQTL_rscripts/header.txt")
colnames(perm_pbs_pc5) <- header
colnames(perm_fnf_pc4) <- header

## check first beta approximated P-values are okay

pdf(file = "output/results_plots/beta_approximation_p_value_check_Both_cond.pdf",   # The directory you want to save the file in
    width = 11, # The width of the plot in inches
    height = 6.5) # The height of the plot in inches

par(mfrow=c(1,2))
plot(perm_pbs_pc5$adj_emp_pval, perm_pbs_pc5$adj_beta_pval, xlab="Direct method", ylab="Beta approximation", main="PBS-Beta approximation check (PC5)")
abline(0, 1, col="red")

plot(perm_fnf_pc4$adj_emp_pval, perm_fnf_pc4$adj_beta_pval, 
     xlab="Direct method", 
     ylab="Beta approximation", 
     main="FNF-Beta approximation check (PC4)")
abline(0, 1, col="red") # Add a red line to this plot as well
par(mfrow=c(1,1))

dev.off()


#-------------------------------------------------------------------------------
#Finding the significant 
header <- read.table("scripts/sQTL_rscripts/header.txt")
new_header <- cbind(header, "qval","therehod")
sigQTL_pbs <- read.table("output/01.qtltools_re/01.significant/pbs_0.05_pc5.significant.txt") |> as.data.frame()
colnames(sigQTL_pbs) <- new_header


hg38_intron <- fread("/work/users/s/e/seyoun/CQTL_sQTL/tools/leafcutter/gencode_hg38_all_introns.bed.gz")
colnames(hg38_intron) <- c("chr","start","end","symbol","ensg","strand","enst","id","region","info")
#hg38_intron_filtered <- hg38_intron %>%
#  filter(!startsWith(symbol, "ENSG"))

hg38_intron_sub <- hg38_intron %>%
  mutate(genomicLoc = paste(chr, start, end, sep = ":")) %>%
  dplyr::select(genomicLoc, ensg)

hg38_intron_sub_select_first <- hg38_intron_sub %>%
  group_by(genomicLoc) %>%
  arrange(ensg) %>%  # Ensure that the data is sorted by ensg within each genomicLoc
  dplyr::slice(1) %>%       # Select the first ensg for each genomicLoc
  ungroup()


sig_sGene_table_pbs <- sigQTL_pbs %>%
  mutate(
    clusterID=gsub("\\..*$", "", sapply(strsplit(as.character(phe_id), ":"), `[`, 4)),
    genomicLoc=sapply(strsplit(as.character(phe_id), ":"), function(x) paste(x[1:3], collapse = ":"))
  )


sig_sGene_table_pbs_gene_annot <- sig_sGene_table_pbs %>%
  left_join(hg38_intron_sub_select_first, by = "genomicLoc") 

sig_sGene_table_pbs_noensg <- sig_sGene_table_pbs_gene_annot[is.na(sig_sGene_table_pbs_gene_annot$ensg),]

# second trial to find the NA genes

leafcutter_pheno <- fread("output/gtex_cluster/ctl_fnf.leafcutter.phenotype_groups.txt",header=F)
colnames(leafcutter_pheno) <- c("id_ensg","ensg")
leafcutter_pheno_mutated <- leafcutter_pheno %>%
  mutate(
    chr=gsub("\\..*$", "", sapply(strsplit(as.character(id_ensg), ":"), `[`, 1)),
    start=gsub("\\..*$", "", sapply(strsplit(as.character(id_ensg), ":"), `[`, 2)),
    end=gsub("\\..*$", "", sapply(strsplit(as.character(id_ensg), ":"), `[`, 3)),
    clusterID=gsub("\\..*$", "", sapply(strsplit(as.character(id_ensg), ":"), `[`, 4)),
    phe_id=sapply(strsplit(as.character(id_ensg), ":"), function(x) paste(x[1:4], collapse = ":")),
    genomicLoc=sapply(strsplit(as.character(id_ensg), ":"), function(x) paste(x[1:3], collapse = ":"))
  )
leafcutter_pheno_subset <- leafcutter_pheno_mutated %>%
  dplyr::select(ensg, clusterID) %>%
  group_by(clusterID) %>%
  arrange(ensg) %>%  # Ensure that the data is sorted by ensg within each genomicLoc
  dplyr::slice(1) %>%       # Select the first ensg for each genomicLoc
  ungroup()

sig_sGene_table_pbs_add_gene_anno_v2 <- sig_sGene_table_pbs_noensg %>%
  left_join(leafcutter_pheno_subset,by = "clusterID") %>%
  mutate(ensg = coalesce(ensg.x, ensg.y)) %>%
  dplyr::select(-ensg.x, -ensg.y)

sig_sGene_table_pbs_add_nogene_v3 <- sig_sGene_table_pbs_add_gene_anno_v2[is.na(sig_sGene_table_pbs_add_gene_anno_v2$ensg),]

#Third trial to find the genes
sig_noname_gr <- GRanges(
  seqnames = sig_sGene_table_pbs_add_nogene_v3$phe_chr,
  ranges = IRanges(start = sig_sGene_table_pbs_add_nogene_v3$phe_from-1, end = sig_sGene_table_pbs_add_nogene_v3$phe_to)
)
#txdb <- makeTxDbFromGFF(file="/work/users/s/e/seyoun/Ref/genome/gencode.v45.annotation.gtf")
#saveDb(x=txdb, file = "gencode.v45.annotation.TxDb")
txdb <- loadDb("../crispr/02.test_seq/gencode.v45.annotation.TxDb")
txdb_genes <- genes(txdb)

overlaps <- mergeByOverlaps(sig_noname_gr, txdb_genes,type="within") 

overlaps_df <- data.frame(
  #phe_chr = seqnames(overlaps$sig_noname_gr),
  #phe_from = start(overlaps$sig_noname_gr),
  #phe_to = end(overlaps$sig_noname_gr),
  gene_id = overlaps$gene_id,
  genomicLoc = paste0(seqnames(overlaps$sig_noname_gr), ":", start(overlaps$sig_noname_gr), ":", end(overlaps$sig_noname_gr))
)

overlaps_df_selected_first <- overlaps_df %>% 
  group_by(genomicLoc) %>%
  arrange(gene_id) %>%  # Ensure that the data is sorted by ensg within each genomicLoc
  dplyr::slice(1) %>%       # Select the first ensg for each genomicLoc
  ungroup()

sig_sGene_table_pbs_add_gene_v4 <- sig_sGene_table_pbs_add_nogene_v3 %>%
  left_join(overlaps_df_selected_first, by = "genomicLoc") %>%
  mutate(ensg = coalesce(gene_id, ensg)) %>%
  dplyr::select(-gene_id) 


# merge all 
sig_sGene_annot_final <- rbind(sig_sGene_table_pbs_gene_annot[!is.na(sig_sGene_table_pbs_gene_annot$ensg),],
sig_sGene_table_pbs_add_gene_anno_v2[!is.na(sig_sGene_table_pbs_add_gene_anno_v2$ensg),], 
sig_sGene_table_pbs_add_gene_v4)

##fnf
sigQTL_fnf <- read.table("output/01.qtltools_re/01.significant/fnf_0.05_pc4.significant.txt") |> as.data.frame()
colnames(sigQTL_fnf) <- new_header

sig_sGene_table_fnf <- sigQTL_fnf %>%
  mutate(
    clusterID = gsub("\\..*$", "", sapply(strsplit(as.character(phe_id), ":"), `[`, 4)),
    genomicLoc = sapply(strsplit(as.character(phe_id), ":"), function(x) paste(x[1:3], collapse = ":"))
  )


sig_sGene_fnf_annot_v1 <- sig_sGene_table_fnf %>%
  left_join(hg38_intron_sub_select_first, by = "genomicLoc") 

sig_sGene_fnf_NOensg_v1 <- sig_sGene_fnf_annot_v1[is.na(sig_sGene_fnf_annot_v1$ensg),]

#second trial-fnf
sig_sGene_fnf_anno_v2 <- sig_sGene_fnf_NOensg_v1 %>%
  left_join(leafcutter_pheno_subset,by = "clusterID") %>%
  mutate(ensg = coalesce(ensg.x, ensg.y)) %>%
  dplyr::select(-ensg.x, -ensg.y)

sig_sGene_fnf_Noensg_v2 <- sig_sGene_fnf_anno_v2[is.na(sig_sGene_fnf_anno_v2$ensg),]

#third trial-fnf
sig_noname_fnf_gr <- GRanges(
  seqnames = sig_sGene_fnf_Noensg_v2$phe_chr,
  ranges = IRanges(start = sig_sGene_fnf_Noensg_v2$phe_from-1, end = sig_sGene_fnf_Noensg_v2$phe_to)
)
overlaps_fnf <- mergeByOverlaps(sig_noname_fnf_gr, txdb_genes,type="within") 

overlaps_fnf_df <- data.frame(
  #phe_chr = seqnames(overlaps$sig_noname_gr),
  #phe_from = start(overlaps$sig_noname_gr),
  #phe_to = end(overlaps$sig_noname_gr),
  gene_id = overlaps_fnf$gene_id,
  genomicLoc = paste0(seqnames(overlaps_fnf$sig_noname_fnf_gr), ":", start(overlaps_fnf$sig_noname_fnf_gr), ":", end(overlaps_fnf$sig_noname_fnf_gr))
)

overlaps_fnf_df_select_first <- overlaps_fnf_df %>% 
  group_by(genomicLoc) %>%
  arrange(gene_id) %>%  # Ensure that the data is sorted by ensg within each genomicLoc
  dplyr::slice(1) %>%       # Select the first ensg for each genomicLoc
  ungroup()

sig_sGene_fnf_annot_v3 <- sig_sGene_fnf_Noensg_v2 %>%
  left_join(overlaps_fnf_df_select_first, by = "genomicLoc") %>%
  mutate(ensg = coalesce(gene_id, ensg)) %>%
  dplyr::select(-gene_id) 


# merge all 
sig_sGene_fnf_anno_final <- rbind(sig_sGene_fnf_annot_v1[!is.na(sig_sGene_fnf_annot_v1$ensg),],
                               sig_sGene_fnf_anno_v2[!is.na(sig_sGene_fnf_anno_v2$ensg),], 
                               sig_sGene_fnf_annot_v3)


#now annotate the HGNC for the ENSG

library(rtracklayer)
gtf_path <- "/work/users/s/e/seyoun/Ref/genome/gencode.v45.annotation.gtf"
gtf_data <- import(gtf_path)
gene_info <- gtf_data[gtf_data$type == "gene"]

genes_df <- data.frame(
  gene_id = mcols(gtf_data)$gene_id,
  gene_name = mcols(gtf_data)$gene_name,
  stringsAsFactors = FALSE
)

genes_df_unique <- genes_df %>%
  dplyr::rename(ensg = gene_id) %>%
  dplyr::distinct(ensg, gene_name, .keep_all = TRUE)  


sig_sGene_pbs_anno_HGNC_final <- sig_sGene_annot_final %>%
  left_join(genes_df_unique, by="ensg")

sig_sGene_fnf_anno_HGNC_final <- sig_sGene_fnf_anno_final %>%
  left_join(genes_df_unique, by="ensg")



# sig_sGene_table_pbs_multiple_ENSG <- sig_sGene_table_pbs_add_nogene_v3 %>%
#   group_by(phe_id) %>%
#   summarize(distinct_ENSG_count = n_distinct(ensg), .groups = 'drop') %>%
#   dplyr::filter(distinct_ENSG_count > 1)
# 
# sig_sGene_table_pbs_multiple_ENSG_rows <- sig_sGene_table_pbs_add_nogene_v3 %>%
#   dplyr::filter(phe_id %in% sig_sGene_table_pbs_multiple_ENSG$phe_id)

#-----

pbs_sGene_sig <- sig_sGene_pbs_anno_HGNC_final %>%
  dplyr::filter(!is.na(ensg)) %>%
  dplyr::select(ensg) %>%
  unique() %>%
  unlist()

fnf_sGene_sig <- sig_sGene_fnf_anno_HGNC_final %>%
  dplyr::filter(!is.na(ensg)) %>%
  dplyr::select(ensg) %>%
  unique() %>%
  unlist()

library(VennDiagram) 
library(ggvenn)
sgene_list <- list("sGene_PBS"=pbs_sGene_sig, "sGene_FN-f"=fnf_sGene_sig)


  
pbs_fnf_sig_sGenes <- tibble(values = unique(c(pbs_sGene_sig, fnf_sGene_sig))) %>%
  mutate(PBS = values %in% pbs_sGene_sig,
         FNF = values %in% fnf_sGene_sig)

venn_diagram <- ggplot(pbs_fnf_sig_sGenes, aes(A = PBS, B = FNF)) +
  geom_venn(set_names = c("PBS", "FN-f"), 
            fill_color = c("#CBD5E8","#FDCDAC"), 
            stroke_color = NA, auto_scale = TRUE, show_percentage = TRUE,
            text_size = 3, set_name_size = 0) +
  coord_fixed()  +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"))
ggsave(filename = "output/results_plots/sqtl_plots/venn_sGenes.pdf", 
       width = 4, height = 4, units = "in")
save(venn_diagram, file = "output/results_plots/sqtl_plots/venn_sGenes.rda")
venn_font(venn_diagram, font = "Helvetica")


#intron junction 

intron_junction_pbs_unique <-sig_sGene_pbs_anno_HGNC_final$phe_id  %>%
  unique()
intron_junction_fnf_unique <-sig_sGene_fnf_anno_HGNC_final$phe_id  %>%
  unique() 

pbs_fnf_sig_introns <- tibble(values = unique(c(intron_junction_pbs_unique, intron_junction_fnf_unique))) %>%
  mutate(PBS = values %in% intron_junction_pbs_unique,
         FNF = values %in% intron_junction_fnf_unique)

venn_introns_diagram <- ggplot(pbs_fnf_sig_introns, aes(A = PBS, B = FNF)) +
  geom_venn(set_names = c("PBS", "FN-f"), 
            fill_color = c("#CBD5E8","#FDCDAC"), 
            stroke_color = NA, auto_scale = TRUE, show_percentage = FALSE,
            text_size = 3, set_name_size = 0) +
  coord_fixed()  +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"))

ggsave(filename = "output/results_plots/sqtl_plots/venn_sintrons.pdf", 
       width = 4, height = 4, units = "in")
save(venn_introns_diagram, file = "output/results_plots/sqtl_plots/venn_sintrons.rda")
venn_font(venn_introns_diagram, font = "Helvetica")

#-------------------------------------------------------------------------------
#Conditional plots
library(ggpubr)
library(gridExtra)

header_cond <- c("phe_id","phe_chr","phe_from",
                 "phe_to","phe_strd","n_var_in_cis","dist_phe_var","var_id",
                 "var_chr","var_from","var_to","rank",
                 "fwd_pval","fwd_r_squared","fwd_slope","fwd_best_hit","fwd_sig",
                 "bwd_pval","bwd_r_squared","bwd_slope","bwd_best_hit","bwd_sig")
cond_pbs <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/conditional_pbs/conditional_pbs_top_variants.txt")
colnames(cond_pbs) <- header_cond
cond_fnf <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/01.qtltools_re/conditional_fnf/conditional_fnf_top_variants.txt")
colnames(cond_fnf) <- header_cond
# Create a function to prepare the data
prepare_data <- function(counts, dataset_name) {
  ranks <- as.character(0:3)
  df <- data.frame(rank = ranks, count = ifelse(ranks %in% names(counts), as.integer(counts[match(ranks, names(counts))]), 0))
  df$dataset <- dataset_name
  return(df)
}

# Prepare data for both PBS and FNF
dfm <- rbind(prepare_data(rank_pbs_counts, "PBS"), prepare_data(rank_fnf_counts, "FNF"))
dfm <- dfm %>%
  mutate(rank = paste0(rank, "°"))
dfm$rank <- factor(dfm$rank, levels = c("0°", "1°", "2°", "3°"))
dfm_only_cond <- dfm %>% dplyr::filter(rank != "0°")



# Correctly assigning colors to the levels
palette_colors <- c("PBS" = "#74CCEC", "FNF" = "#FAB394")

# Now call ggdotchart
ggdot_conditional <- ggdotchart(dfm_only_cond, x = "rank", y = "count",
           color = "dataset",                          # Column to define color groups
           palette = c("PBS"="#74CCEC", "FNF"="#FAB394"),  # Define colors for each group
           sorting = "descending",                      # Sort dots in descending order
           add = "segments",                            # Add segments from zero to the dot's value
           position = position_dodge(0.3),              # Dodge positions for clarity
           rotate = TRUE,                               # Horizontal plot
           group = "dataset",                           # Group data by 'dataset' for plotting
           dot.size = 3,                                # Size of dots
           ggtheme = theme_pubr()) +  
  labs(title = "Number of intron junctions") + # Use ggpubr theme for aesthetics
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_markdown(size = 6, family = "Helvetica", margin = margin(r = -15)),
        text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.text.x =  element_text(color = "black", size = 6),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        #panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, family = "Helvetica", size = 10, margin = margin(b = -3)))

save(ggdot_conditional, file = "output/results_plots/sqtl_plots/ggdot_conditional.rda")
ggsave(filename =  "output/results_plots/sqtl_plots/ggdot_conditional.pdf",
       plot = ggdot_conditional, width = 5, height = 5, units = "in")


#intron junction counts



#-------------------------------------------------------------------------------
#distance density plot 
library(plyr)

sigQTL_pbs <- sigQTL_pbs %>%
  mutate(distance_from_var = dist_phe_var / 1000)


sigQTL_pbs$group <- ifelse(sigQTL_pbs$distance_from_var < 0, "Below 0", "Above 0")
zero_distances <- sigQTL_pbs %>% 
  dplyr::filter(distance_from_var == 0) %>%
  mutate(group = "Below 0") # Assign "Below 0" for one copy
zero_distances_above <- mutate(zero_distances, group = "Above 0")

sigQTL_pbs_adjusted <- sigQTL_pbs %>%
  dplyr::filter(distance_from_var != 0) %>%
  mutate(group = ifelse(distance_from_var < 0, "Below 0", "Above 0"))

sigQTL_pbs_adjusted <- bind_rows(sigQTL_pbs_adjusted, zero_distances, zero_distances_above)

mu <- ddply(sigQTL_pbs_adjusted, "group", summarise, grp.mean=mean(distance_from_var))

mu$group <- factor(mu$group, levels = c("Below 0","Above 0"))

sigQTL_pbs_adjusted$group <- factor(sigQTL_pbs_adjusted$group, levels = c("Below 0","Above 0"))

# breaks_below <- seq(-100, 20, by = 20)
# breaks_above <- seq(-20, 100, by = 20)
# labels_below <- as.character(breaks_below)
# labels_above <- as.character(breaks_above)
# labels_below[breaks_below > 0] <- NA  # Hide labels from 0 to 20 for 'Below 0'
# labels_above[breaks_above < 0] <- NA  # Hide labels from -20 to 0 for 'Above 0'
# 
# scale_x_below <- scale_x_continuous(limits = c(-100, 20), breaks = breaks_below, labels = labels_below)
# scale_x_above <- scale_x_continuous(limits = c(-20, 100), breaks = breaks_above, labels = labels_above)

library(colorspace)
density_PBS_sqtl <- ggplot(sigQTL_pbs_adjusted, aes(x = distance_from_var, fill = group)) +
  geom_density(alpha = 0.5, na.rm = TRUE,color = NA) +
  scale_fill_manual(values = c("Below 0" = "#0067B9", "Above 0" = "#FCCE52")) +
  labs(title = "PBS-distance between SD/SA  and a lead sQTL", y = "Density", x = "Distance (in kb) to splice site") +
  facet_wrap(~ group, scales = "free_x", ncol = 2) +
  geom_vline(data = mu, aes(xintercept = grp.mean), color = c("#FCCE52", "#0067B9"), linetype = "dashed", size = 0.5) +
  geom_text(data = mu, aes(x = grp.mean, y = 0.075, label = sprintf("Median = %.2f kb", grp.mean)), vjust = -1, color = darken(c("#FCCE52", "#0067B9")),size = 2) +
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
        axis.title.x = element_markdown(size = 6, family = "Helvetica", margin = margin(r = -15)),
        axis.ticks = element_blank(),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_markdown(size = 6, family = "Helvetica", margin = margin(r = -15)),
        text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.text.x =  element_text(color = "black", size = 6),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        #panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, family = "Helvetica", size = 10, margin = margin(b = 3)))


save(density_PBS_sqtl, file = "output/results_plots/sqtl_plots/density_PBS_sqtl.rda")
ggsave(filename =  "output/results_plots/sqtl_plots/density_PBS_sqtl.pdf",
       plot = density_PBS_sqtl, width = 7, height = 4, units = "in")


# FNF

sigQTL_fnf <- sigQTL_fnf %>%
  mutate(distance_from_var = dist_phe_var / 1000)


sigQTL_fnf$group <- ifelse(sigQTL_fnf$distance_from_var < 0, "Below 0", "Above 0")
zero_distances <- sigQTL_fnf %>% 
  dplyr::filter(distance_from_var == 0) %>%
  mutate(group = "Below 0") # Assign "Below 0" for one copy
zero_distances_above <- mutate(zero_distances, group = "Above 0")

sigQTL_fnf_adjusted <- sigQTL_fnf %>%
  dplyr::filter(distance_from_var != 0) %>%
  mutate(group = ifelse(distance_from_var < 0, "Below 0", "Above 0"))

sigQTL_fnf_adjusted <- bind_rows(sigQTL_fnf_adjusted, zero_distances, zero_distances_above)

mu <- ddply(sigQTL_fnf_adjusted, "group", summarise, grp.mean=mean(distance_from_var))

mu$group <- factor(mu$group, levels = c("Below 0","Above 0"))

sigQTL_fnf_adjusted$group <- factor(sigQTL_fnf_adjusted$group, levels = c("Below 0","Above 0"))


density_FNF_sqtl <- ggplot(sigQTL_fnf_adjusted, aes(x = distance_from_var, fill = group)) +
  geom_density(alpha = 0.5, na.rm = TRUE,color = NA) +
  scale_fill_manual(values = c("Below 0" = "#287C6F", "Above 0" = "#E07653")) +
  labs(title = "FN-f-distance between SD/SA  and a lead sQTL", y = "Density", x = "Distance (in kb) to splice site") +
  facet_wrap(~ group, scales = "free_x", ncol = 2) +
  geom_vline(data = mu, aes(xintercept = grp.mean), color = c("#E07653", "#287C6F"), linetype = "dashed", size = 0.5) +
  geom_text(data = mu, aes(x = grp.mean, y = 0.075, label = sprintf("Median = %.2f kb", grp.mean)), vjust = -1, color = darken(c("#E07653", "#287C6F")),size = 2) +
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
        axis.title.x = element_markdown(size = 6, family = "Helvetica", margin = margin(r = -15)),
        axis.ticks = element_blank(),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_markdown(size = 6, family = "Helvetica", margin = margin(r = -15)),
        text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.text.x =  element_text(color = "black", size = 6),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        #panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, family = "Helvetica", size = 10, margin = margin(b = 3)))


save(density_FNF_sqtl, file = "output/results_plots/sqtl_plots/density_FNF_sqtl.rda")
ggsave(filename =  "output/results_plots/sqtl_plots/density_FNF_sqtl.pdf",
       plot = density_FNF_sqtl, width = 7, height = 4, units = "in")


#-------------------------------------------------------------------------------
#distance visualization_2 with boxplot


#------------------------------------------------------------------------------
#mkaing figures
library(plotgardener)
pdf(file = "results_plots/sqtl_plots/Fig3.pdf",   # The directory you want to save the file in
    width = 9, # The width of the plot in inches
    height = 5)
pageCreate(width = 9, height = 5, default.units = "inches",showGuides =FALSE)
plotText("A", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

plotSegments(
      x0 = unit(0.35, "in"),
      y0 = 0.75,
      x1 = unit(2.65, "in"),
      y1 = 0.75,
      default.units = "inches",
      lwd = 2.5, lty = 1,"darkgrey",
      just ="left"
    )
plotRect(x = unit(0.5, "in")  ,
             y = 0.75, width = unit(10, "mm"),
             height = unit(4, "mm"), linecolor = NA, fill = "#2473b5",
             just = "left")

plotRect(x = unit(1.5, "in")  ,
         y = 0.75, width = unit(10, "mm"),
         height = unit(4, "mm"), linecolor = NA, fill = "#2473b5",
         just = "center")

plotRect(x = unit(2.5, "in")  ,
         y = 0.75, width = unit(10, "mm"),
         height = unit(4, "mm"), linecolor = NA, fill = "#2473b5",
         just = "right")

plotSegments(
  x0 = unit(0.5, "in")+unit(10, "mm"),
  y0 = 0.75,
  x1 = unit(1.5, "in"),
  y1 = 0.45,
  default.units = "inches",
  lwd = 1, lty = 2,"#C43D4D",
  just ="left"
)

plotSegments(
  x0 = unit(2.5, "in")-unit(10, "mm"),
  y0 = 0.75,
  x1 = unit(1.5, "in"),
  y1 = 0.45,
  default.units = "inches",
  lwd = 1, lty = 2,"#C43D4D",
  just ="right"
)


plotSegments(
  x0 = unit(1.5, "in")-unit(5, "mm"),
  y0 = 0.75,
  x1 = unit(0.5, "in")+unit(15, "mm"),
  y1 = 0.65,
  default.units = "inches",
  lwd = 1, lty = 2,"#C43D4D",
  just ="right"
)

plotSegments(
  x0 = unit(1.5, "in")+unit(5, "mm"),
  y0 = 0.75,
  x1 = unit(2.5, "in")-unit(15, "mm"),
  y1 = 0.65,
  default.units = "inches",
  lwd = 1, lty = 2,"#C43D4D",
  just ="left"
)

plotSegments(
  x0 = 0.2, y0 = 1.1, x1 = 2.8, y1 = 1.1,
  default.units = "inches",
  arrow = arrow(ends="both",type = "closed",angle = 20,length = unit(0.075, "inches")),
  lwd = 0.5, "#C43D4D",fill="#C43D4D", lty=2
)

# plotSegments(
#   x0 = unit(0.5, "in")+unit(10, "mm"), y0 = 0.35, 
#   x1 =unit(0.5, "in")+unit(10, "mm"), y1 = 1.9,
#   default.units = "inches",
#   lwd = 0.5, "#C43D4D",lty=2
# )
# plotSegments(
#   x0 = unit(2.5, "in")-unit(10, "mm"), y0 = 0.35, 
#   x1 =unit(2.5, "in")-unit(10, "mm"), y1 = 1.9,
#   default.units = "inches",
#   lwd = 0.5, "#C43D4D",lty=2
# )
plotText(label = 'SD', x = unit(0.5, "in")+unit(10, "mm"),
         y = 0.9,
         fontsize = 7, fontfamily = "Helvetica",
         fontface = "bold",
         just ="center", fontcolor = "#434e57")

plotText(label = 'SA', x = unit(1.5, "in")-unit(5, "mm"),
         y = 0.9,
         fontface = "bold",
         fontsize = 7, fontfamily = "Helvetica",
         just ="center", fontcolor = "#434e57")

plotText(label = 'SD', x = unit(1.5, "in")+unit(5, "mm"),
         y = 0.9,
         fontface = "bold",
         fontsize = 7, fontfamily = "Helvetica",
         just ="center", fontcolor = "#434e57")

plotText(label = 'SA', x = unit(2.5, "in")-unit(10, "mm"),
         y = 0.9,
         fontface = "bold",
         fontsize = 7, fontfamily = "Helvetica",
         just ="center", fontcolor = "#434e57")

plotText(label = 'Cluster', x = unit(1.5, "in"),
         y = 0.35,
         fontsize = 7, fontfamily = "Helvetica",
         just ="center", fontcolor = "#434e57")


plotText(label = '100 kb', x = unit(0.5, "in")+unit(10, "mm"),
         y = 1.1,
         fontsize = 7, fontfamily = "Helvetica",
         fontface = "bold",
         just ="right", fontcolor = "#C43D4D")

plotText(label = '100 kb', x = unit(2.5, "in")-unit(10, "mm"),
         y = 1.1,
         fontsize = 7, fontfamily = "Helvetica",
         fontface = "bold",
         just ="left", fontcolor = "#C43D4D")


plotSegments(
  x0 = 0.2, y0 = 1.1, 
  x1 =0.2, y1 = 1.25,
  default.units = "inches",
  lwd = 1, "grey70",lty=1
)

plotSegments(
  x0 = 2.8, y0 = 1.1, 
  x1 =2.8, y1 = 1.25,
  default.units = "inches",
  lwd = 1, "grey70",lty=1
)

plotSegments(
  x0 = 0.2, y0 = 1.25, 
  x1 =2.8, y1 = 1.25,
  default.units = "inches",
  lwd = 1, "grey70",lty=1
)

plotText(label = 'sQTL testing range 1', x = unit(1.5, "in"),
         y = 1.3,
         fontsize = 7, fontfamily = "Helvetica",
         just ="center", fontcolor = "#434e57")

#Test range2

plotSegments(
  x0 = 0.2, y0 = 1.5, x1 = 1.9, y1 = 1.5,
  default.units = "inches",
  arrow = arrow(ends="both",type = "closed",angle = 20,length = unit(0.075, "inches")),
  lwd = 0.5, "#C43D4D",fill="#C43D4D", lty=2
)


plotText(label = '100 kb', x = unit(0.5, "in")+unit(10, "mm"),
         y = 1.5,
         fontsize = 7, fontfamily = "Helvetica",
         fontface = "bold",
         just ="right", fontcolor = "#C43D4D")

plotText(label = '100 kb', x = unit(1.5, "in")-unit(5, "mm"),
         y = 1.5,
         fontsize = 7, fontfamily = "Helvetica",
         fontface = "bold",
         just ="left", fontcolor = "#C43D4D")


plotSegments(
  x0 = 0.2, y0 = 1.5, 
  x1 =0.2, y1 = 1.65,
  default.units = "inches",
  lwd = 1, "grey70",lty=1
)

plotSegments(
  x0 = 1.9, y0 = 1.5, 
  x1 =1.9, y1 = 1.65,
  default.units = "inches",
  lwd = 1, "grey70",lty=1
)

plotSegments(
  x0 = 0.2, y0 = 1.65, 
  x1 =1.9, y1 = 1.65,
  default.units = "inches",
  lwd = 1, "grey70",lty=1
)

plotText(label = 'sQTL testing range 2', x = unit(1.1, "in"),
         y = 1.7,
         fontsize = 7, fontfamily = "Helvetica",
         just ="center", fontcolor = "#434e57")

#Test range 3

plotSegments(
  x0 = 1, y0 = 1.9, x1 = 2.8, y1 = 1.9,
  default.units = "inches",
  arrow = arrow(ends="both",type = "closed",angle = 20,length = unit(0.075, "inches")),
  lwd = 0.5, "#C43D4D",fill="#C43D4D", lty=2
)


plotText(label = '100 kb', x = unit(1.5, "in")+unit(5, "mm"),
         y = 1.9,
         fontsize = 7, fontfamily = "Helvetica",
         fontface = "bold",
         just ="right", fontcolor = "#C43D4D")

plotText(label = '100 kb', x = unit(2.5, "in")-unit(10, "mm"),
         y = 1.9,
         fontsize = 7, fontfamily = "Helvetica",
         fontface = "bold",
         just ="left", fontcolor = "#C43D4D")


plotSegments(
  x0 = 1, y0 = 1.9, 
  x1 =1, y1 = 2.05,
  default.units = "inches",
  lwd = 1, "grey70",lty=1
)

plotSegments(
  x0 = 2.8, y0 = 1.9, 
  x1 =2.8, y1 = 2.05,
  default.units = "inches",
  lwd = 1, "grey70",lty=1
)

plotSegments(
  x0 = 1, y0 = 2.05, 
  x1 =2.8, y1 = 2.05,
  default.units = "inches",
  lwd = 1, "grey70",lty=1
)

plotText(label = 'sQTL testing range 3', x = unit(1.85, "in"),
         y = 2.1,
         fontsize = 7, fontfamily = "Helvetica",
         just ="center", fontcolor = "#434e57")

# plotSegments(
#   x0 = unit(0.8, "in"),
#   y0 = 1.1,
#   x1 = unit(1.3, "in"),
#   y1 = 1.1,
#   default.units = "inches",
#   lwd = 1, lty = 1,"darkgrey",
#   just ="left"
# )
# 
# plotRect(x = unit(0.9, "in")  ,
#          y = 1.1, width = unit(0.2, "in"),
#          height = unit(2, "mm"), linecolor = NA, fill = "#2473b5",
#          just = "right")
# 
# plotRect(x = unit(1.5, "in")-unit(5, "mm"),
#          y = 1.1, width = unit(0.2, "in"),
#          height = unit(2, "mm"), linecolor = NA, fill = "#2473b5",
#          just = "left")
# 
# plotSegments(
#   x0 = 0.7, y0 = 1.1, x1 = 0.4, y1 = 1.1,
#   default.units = "inches",
#   arrow = arrow(type = "closed",angle = 20,length = unit(0.075, "inches")),
#   lwd = 0.3, fill="#C43D4D",lty=2
# )
# 
# plotSegments(
#   x0 = 1.5, y0 = 1.1, x1 = 1.8, y1 = 1.1,
#   default.units = "inches",
#   arrow = arrow(type = "closed",angle = 20,length = unit(0.075, "inches")),
#   lwd = 0.3, fill="#C43D4D",lty=2
# )
# 
# 
# plotText(label = '100 kb', x = unit(0.4, "in"),
#          y = 1.2,
#          fontsize = 7, fontfamily = "Helvetica",
#          just ="left", fontcolor = "#C43D4D")
# 
# plotText(label = '100 kb', x = unit(1.5, "in"),
#          y = 1.2,
#          fontsize = 7, fontfamily = "Helvetica",
#          just ="left", fontcolor = "#C43D4D")
# 
# 
# plotSegments(
#   x0 = unit(0.8, "in"),
#   y0 = 1.3,
#   x1 = unit(2.2, "in"),
#   y1 = 1.3,
#   default.units = "inches",
#   lwd = 1, lty = 1,"darkgrey",
#   just ="left"
# )
# 
# plotRect(x = unit(0.9, "in")  ,
#          y = 1.3, width = unit(0.2, "in"),
#          height = unit(2, "mm"), linecolor = NA, fill = "#2473b5",
#          just = "right")
# 
# plotRect(x = unit(2.1, "in"),
#          y = 1.3, width = unit(0.2, "in"),
#          height = unit(2, "mm"), linecolor = NA, fill = "#2473b5",
#          just = "left")
# 
# plotSegments(
#   x0 = 0.7, y0 = 1.3, x1 = 0.4, y1 = 1.3,
#   default.units = "inches",
#   arrow = arrow(type = "closed",angle = 20,length = unit(0.075, "inches")),
#   lwd = 0.3, fill="#C43D4D",lty=2
# )
# 
# plotSegments(
#   x0 = 2.3, y0 = 1.3, x1 = 2.6, y1 = 1.3,
#   default.units = "inches",
#   arrow = arrow(type = "closed",angle = 20,length = unit(0.075, "inches")),
#   lwd = 0.3, fill="#C43D4D",lty=2
# )
# 
# 
# plotText(label = '100 kb', x = unit(0.4, "in"),
#          y = 1.4,
#          fontsize = 7, fontfamily = "Helvetica",
#          just ="left", fontcolor = "#C43D4D")
# 
# plotText(label = '100 kb', x = unit(2.3, "in"),
#          y = 1.4,
#          fontsize = 7, fontfamily = "Helvetica",
#          just ="left", fontcolor = "#C43D4D")

#-------------------------------------------------------------------------------

plotText("B", x = 3.45, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

dev.off()

load(file="output/results_plots/sqtl_plots/venn_sintrons.rda")
plotGG(venn_introns_diagram, x = 0, y = 2, width = 3, height = 2.2)
plotText("Intron junction", x = 1.5, y = 2.4, just = c("center", "top"), fontfamily = "Helvetica",
         fontsize = 8, fontface = "bold")

c("#CBD5E8","#FDCDAC")
plotText("PBS", x = 0.5, y = 3, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 8, fontface = "bold", fontcolor = darken("#CBD5E8"))
plotText("FN-f", x = 2.25, y = 3, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 8, fontface = "bold", fontcolor = darken("#FDCDAC"))
plotGG(venn_diagram, x = 0.5, y = 3.0, width = 2.5, height = 1.5)

#-------------------------------------------------------------------------------

plotText("C", x = 3, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")


load(file="output/results_plots/sqtl_plots/density_PBS_sqtl.rda")
#load(file="output/results_plots/sqtl_plots/density_FNF_sqtl.rda")
#plotGG(density_PBS_sqtl, x = 3.2, y = 0.25, width = 3.9, height = 2.3)
#plotGG(density_FNF_sqtl, x = 3.2, y = 2.6, width = 3.9, height = 2.3)

