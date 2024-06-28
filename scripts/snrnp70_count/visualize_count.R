#Count for SNRNP70 for 101 samples
## Author: Seyoun Byun
## Date: 06.11.2024
## Edited:
#-------------------------------------------------------------------------------
setwd("/work/users/s/e/seyoun/CQTL_sQTL/output")
library(data.table)
library(DESeq2)
library(EBSEA)
library(tidyr)
library(ggplot2)
library(ggtext)
library(colorspace)

#flatten the SNRNP70 GTF TO saf and then 
#awk 'NR==FNR{ids[$1]; next} FNR==1 || $1 in ids' geneid.txt snrnp70_flatten.gtf > snrnp70_flatten_subset.gtf

exon_count <- fread("featurecounts_snrnp70/snrnp70counts.txt")
colnames(exon_count) <- gsub("_sorted.bam_snrnp70.bam","",gsub("/work/users/s/e/seyoun/CQTL_sQTL/output/snrnp70_subset_bam/","",colnames(exon_count)))
colnames(exon_count)[1] <- c("exonid")


#exclude samples for  ['AM7352', 'AM7244']
exon_count_subset <- exon_count %>% dplyr::select(-matches("AM7352|AM7244")) 

exon_count_subset <- exon_count_subset %>%
  mutate(exonid = factor(exonid, levels  = c("ENSE00001555728.1", "ENSE00003535913.1", "ENSE00003567241.1",
                                             "ENSE00003499363.1", "ENSE00003552141.1", "ENSE00003683920.1",
                                             "ENSE00002503116.1", "ENSE00001549580.1", "ENSE00003593260.1",
                                             "ENSE00003614741.1", "ENSE00003179801.2")))

setorder(exon_count_subset, exonid)

# subset of CTL vs FNF and OA 
count_matrix <- exon_count_subset %>% dplyr::select(contains("CTL"), contains("FNF"),contains("OA")) |> 
  as.data.frame()
rownames(count_matrix) <- paste0("SNRNP70:",exon_count_subset$exonid)

col_data <- data.frame(cbind(colnames(count_matrix),do.call(rbind,strsplit(colnames(count_matrix),"_"))))
colnames(col_data) <- c("id","donor","condition","rep","sex")

col_data$condition <- factor(col_data$condition, levels = c("CTL","FNF","OA"))

original_order <- rownames(count_matrix)
design <- ~condition
ebsea.out <- EBSEA(count_matrix, col_data, design)
ebsea.out$Group <- factor(col_data$condition, levels = c("CTL","FNF","OA"))


# Dougle checking the order 
if (!identical(rownames(ebsea.out$NormCounts), original_order)) {
  ebsea.out$NormCounts <- ebsea.out$NormCounts[match(original_order, rownames(ebsea.out$NormCounts)), ]
}

if (!identical(rownames(ebsea.out$ExonTable), original_order)) {
  ebsea.out$ExonTable <- ebsea.out$ExonTable[match(original_order, rownames(ebsea.out$ExonTable)), ]
}


group <- as.factor(ebsea.out$Group)
exon.table <- ebsea.out$ExonTable
exon.table$GeneExon <- row.names(exon.table)
gene <- 'SNRNP70'


counts <- ebsea.out$NormCounts[grep(gene, row.names(ebsea.out$NormCounts),
                                    fixed = TRUE), , drop = FALSE]


# taking the mean and standard error of the rows of both samples
mean.count <- data.frame('Group1' = apply(counts[, which(group %in% levels(group)[1]), drop = FALSE], 1, mean),
                         'Group2' = apply(counts[, which(group %in% levels(group)[2]), drop = FALSE], 1, mean),
                         'Group3' = apply(counts[, which(group %in% levels(group)[3]), drop = FALSE], 1, mean))
se.count <- data.frame('Group1' = apply(counts[, group %in% levels(group)[1], drop = FALSE],
                                        1, function(x)sd(x)/sqrt(length(x))),
                       'Group2' = apply(counts[, group %in% levels(group)[2], drop = FALSE],
                                        1, function(x)sd(x)/sqrt(length(x))),
                       'Group3' = apply(counts[, group %in% levels(group)[3], drop = FALSE],
                                        1, function(x)sd(x)/sqrt(length(x))))
mean.count$Exon <- c("ex1", "ex2", "ex3", "ex4", "ex5", "ex6", "ex7", "altex8", "ex8", "ex9", "ex10")
mean.count_long <- pivot_longer(mean.count, cols = c(Group1, Group2,Group3), names_to = "Group", values_to = "MeanCount")
mean.count_long <- mean.count_long %>%
  mutate(Group = case_when(
    Group == "Group1" ~ "PBS",
    Group == "Group2" ~ "FNF",
    Group == "Group3" ~ "OA"
  ),
  Group = factor(Group, levels = c("PBS", "FNF","OA"))) # Set factor levels in the specific order


# Ensure the Exon column is ordered correctly
mean.count_long$Exon <- factor(mean.count_long$Exon,
                               levels = c("ex1", "ex2", "ex3", "ex4", "ex5", "ex6", "ex7", "altex8", "ex8", "ex9", "ex10"))
# Now plot with ggplot2
barplot_mean_counts <- ggplot(mean.count_long , aes(x = Exon, y = MeanCount, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.7) +
  labs(x = "Exons", y = "Exon Normalized Counts") +
  scale_fill_manual(values = c("PBS" = "grey90", "FNF" = "grey50","OA" = "grey10"),  # Custom colors
                    labels = c("PBS" = "PBS", "FNF" = "FN-f", 
                               "OA" = "OA")) +
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_markdown(size = 6, family = "Helvetica",
                                        margin = margin(r = -15)),
        text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.text.x = element_text(color = "black", size = 8, margin = margin(b = -1)),
        strip.background = element_blank(),
        #strip.text = element_blank(),
        strip.text.x.bottom = element_markdown(size = 8, margin = margin(t = 1)),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        legend.position = c(0.1,0.9),
        legend.title=element_blank(),
        legend.key.size = unit(0.25, 'cm'),
        legend.text = element_text(size = 5))


#------------------------------------------------------------------------------
#logplot
#-------------------------------------------------------------------------------
# CTL VS FNF
#-------------------------------------------------------------------------------
count_matrix_fnf <- exon_count_subset %>% dplyr::select(contains("CTL"), contains("FNF")) |> 
  as.data.frame()
rownames(count_matrix_fnf) <- paste0("SNRNP70:",exon_count_subset$exonid)

col_data_fnf <- data.frame(cbind(colnames(count_matrix_fnf),do.call(rbind,strsplit(colnames(count_matrix_fnf),"_"))))
colnames(col_data_fnf) <- c("id","donor","condition","rep","sex")

col_data_fnf$condition <- factor(col_data_fnf$condition, levels = c("CTL","FNF"))

original_order <- rownames(count_matrix_fnf)
design <- ~condition
ebsea.out_fnf <- EBSEA(count_matrix_fnf, col_data_fnf, design)
ebsea.out_fnf$Group <- factor(col_data_fnf$condition, levels = c("CTL","FNF"))

ebsea.out_fnf$ExonTable |> head()
ebsea.out_fnf$GeneTable |> head()


if (!identical(rownames(ebsea.out_fnf$NormCounts), original_order)) {
  ebsea.out_fnf$NormCounts <- ebsea.out_fnf$NormCounts[match(original_order, rownames(ebsea.out_fnf$NormCounts)), ]
}

if (!identical(rownames(ebsea.out_fnf$ExonTable), original_order)) {
  ebsea.out_fnf$ExonTable <- ebsea.out_fnf$ExonTable[match(original_order, rownames(ebsea.out_fnf$ExonTable)), ]
}
group_fnf <- as.factor(ebsea.out_fnf$Group)
exon.table_fnf <- ebsea.out_fnf$ExonTable
exon.table_fnf$GeneExon <- row.names(exon.table_fnf)
gene <- 'SNRNP70'
#Fetching Information
gene.exons_fnf <- exon.table_fnf[grep(gene, exon.table_fnf$GeneExon, fixed = TRUE), , drop = FALSE]
if(nrow(gene.exons_fnf) == 0) {
  stop('The gene provided is not in the list')
}

gene.info_fnf <- ebsea.out_fnf$GeneTable[grep(gene, ebsea.out_fnf$GeneTable$Gene, fixed = TRUE),
                                 , drop = FALSE]

# checking if the p-values have NA values. If NA is found it is converted to 1
gene.exons_fnf$pvalue[is.na(gene.exons_fnf$pvalue)] <- 1

# to make sure the order is this 



gene.exons_fnf$Exon <- c("ex1", "ex2", "ex3", "ex4", "ex5", "ex6", "ex7", "altex8", "ex8", "ex9", "ex10")
gene.exons_fnf$Exon <- factor(gene.exons_fnf$Exon,
                          levels = c("ex1", "ex2", "ex3", "ex4", "ex5", "ex6", "ex7", "altex8", "ex8", "ex9", "ex10"))
gene.exons_fnf$category <- ifelse(gene.exons_fnf$padj <= 0.001, '***',
                         ifelse(gene.exons_fnf$padj <= 0.01, '**',
                                ifelse(gene.exons_fnf$padj <= 0.05, '*', '')))
gene.exons_fnf <- gene.exons_fnf %>%
  mutate(status = case_when(
    log2FoldChange > 0 ~ "Upregulated",
    log2FoldChange < 0 ~ "Downregulated",
    TRUE ~ "Static"  # If log2FoldChange == 0 or any other cases you want to consider as Static
  ))



#-------------------------------------------------------------------------------
# CTL VS OA 
#-------------------------------------------------------------------------------
count_matrix_OA <- exon_count_subset %>% dplyr::select(contains("CTL"), contains("OA")) |> 
  as.data.frame()
rownames(count_matrix_OA) <- paste0("SNRNP70:",exon_count_subset$exonid)

col_data_oa <- data.frame(cbind(colnames(count_matrix_OA),do.call(rbind,strsplit(colnames(count_matrix_OA),"_"))))
colnames(col_data_oa) <- c("id","donor","condition","rep","sex")

col_data_oa$condition <- factor(col_data_oa$condition, levels = c("CTL","OA"))

original_order <- rownames(count_matrix_OA)
design <- ~condition
ebsea.out_oa <- EBSEA(count_matrix_OA, col_data_oa, design)
ebsea.out_oa$Group <- factor(col_data_oa$condition, levels = c("CTL","OA"))
ebsea.out_oa$ExonTable |> head()
ebsea.out_oa$GeneTable |> head()


if (!identical(rownames(ebsea.out_oa$NormCounts), original_order)) {
  ebsea.out_oa$NormCounts <- ebsea.out_oa$NormCounts[match(original_order, rownames(ebsea.out_oa$NormCounts)), ]
}

if (!identical(rownames(ebsea.out_oa$ExonTable), original_order)) {
  ebsea.out_oa$ExonTable <- ebsea.out_oa$ExonTable[match(original_order, rownames(ebsea.out_oa$ExonTable)), ]
}


group_oa <- as.factor(ebsea.out_oa$Group)
exon.table_oa <- ebsea.out_oa$ExonTable
exon.table_oa$GeneExon <- row.names(exon.table_oa)
gene <- 'SNRNP70'
# Fetching Information
gene.exons_oa <- exon.table_oa[grep(gene, exon.table_oa$GeneExon, fixed = TRUE), , drop = FALSE]
if(nrow(gene.exons_oa) == 0) {
  stop('The gene provided is not in the list')
}

gene.info_oa <- ebsea.out_oa$GeneTable[grep(gene, ebsea.out_oa$GeneTable$Gene, fixed = TRUE),
                                 , drop = FALSE]


gene.exons_oa$pvalue[is.na(gene.exons_oa$pvalue)] <- 1

# to make sure the order is this 

gene.exons_oa$Exon <- c("ex1", "ex2", "ex3", "ex4", "ex5", "ex6", "ex7", "altex8", "ex8", "ex9", "ex10")
gene.exons_oa$Exon <- factor(gene.exons_oa$Exon,
                          levels = c("ex1", "ex2", "ex3", "ex4", "ex5", "ex6", "ex7", "altex8", "ex8", "ex9", "ex10"))
gene.exons_oa$category <- ifelse(gene.exons_oa$padj <= 0.001, '***',
                         ifelse(gene.exons_oa$padj <= 0.01, '**',
                                ifelse(gene.exons_oa$padj <= 0.05, '*', '')))
gene.exons_oa <- gene.exons_oa %>%
  mutate(status = case_when(
    log2FoldChange > 0 ~ "Upregulated",
    log2FoldChange < 0 ~ "Downregulated",
    TRUE ~ "Static"  # If log2FoldChange == 0 or any other cases you want to consider as Static
  ))

gene.exons_oa <- gene.exons_oa %>% mutate(group="PBS vs. OA")
gene.exons_fnf <- gene.exons_fnf %>% mutate(group="PBS vs. FN-f")
gene_l2fc_all <- rbind(gene.exons_oa,gene.exons_fnf)

dodge_width=0.4
log2_plot <- ggplot(gene_l2fc_all, aes(x = Exon, y = log2FoldChange, fill = group)) +
  geom_col( position = position_dodge(width = dodge_width),width = 0.4) +  # Use geom_col here
  scale_fill_manual(values = c("PBS vs. FN-f" = "grey50", "PBS vs. OA" = "grey10"),
                    labels = c("PBS vs. FN-f" = "PBS vs. FN-f", "PBS vs. OA" = "PBS vs. OA")) +
  theme_classic() +
  scale_y_continuous(name = "Exons log~2~fold change", limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  coord_cartesian(clip = "off") +
  geom_hline(yintercept = 0, lty = 2, color = "grey25", linewidth = 0.25) +
  theme(strip.placement = "outside",
        axis.line = element_line(linewidth = 0.25),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_markdown(size = 6, family = "Helvetica",
                                        margin = margin(r = -15)),
        text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color = "black", size = 6),
        axis.text.x = element_text(color = "black", size = 8, margin = margin(b = -1)),
        strip.background = element_blank(),
        #strip.text = element_blank(),
        strip.text.x.bottom = element_markdown(size = 8, margin = margin(t = 1)),
        panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid = element_blank(),
        legend.position = c(0.1,0.9),
        legend.title=element_blank(),
        legend.key.size = unit(0.25, 'cm'),
        legend.text = element_text(size = 5),
        plot.title = element_text(hjust = 0.5, family = "Helvetica",
                                  size = 10, margin = margin(b=-3)))
#ggtitle("Exon-specific regulation SNRNP70 ")

#log2_plot_text <- log2_plot + geom_text(data = subset(gene_l2fc_all, category != ""),
#                                        aes(label = category, y = log2FoldChange-0.02 ),  # Adjust y offset as needed
#                                        position = position_dodge(width = dodge_width),
#                                        vjust = 0 ,color = "red", size = 2.5)

# Display the plot
library(patchwork)

combined_plot <- log2_plot / barplot_mean_counts

# Optionally, you can align the plots along the x-axis
combined_plot <- combined_plot & theme(
  plot.title.position = "plot",
  panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel background
  plot.background = element_rect(fill = "transparent", color = NA),  # Transparent plot background
  legend.background = element_rect(fill = "transparent", color = NA)  # Optional: Transparent legend background
)


# Print the combined plot
print(combined_plot)

ggsave(filename = "results_plots/Figure2_SNRNP70/snrnp70_combinedplot.pdf",
       plot =combined_plot, width = 5, height = 5, units = "in")
save(combined_plot, file = "results_plots/Figure2_SNRNP70/snrnp70_combinedplot_r43.rda")

