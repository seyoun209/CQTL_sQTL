# Differential Gene expression
## Author: Seyoun Byun
##Date:06.20.2024
#-------------------------------------------------------------------------------
setwd("/work/users/s/e/seyoun/CQTL_sQTL/")

## Load required libraries
library(data.table)
library(tximeta)
library(DESeq2)
library(dplyr)
library(plyranges)
library(liftOver)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


meta_cqtl <- fread("output/clu_fnf/meta_cqtl")
meta_cqtl_subset <- meta_cqtl[!meta_cqtl$Donor %in% config$samples_to_omit, ]
meta_ctl_oa <- meta_cqtl_subset %>% dplyr::filter(Condition %in% c('CTL','OA'))


## Add quant paths and names
meta_ctl_oa$files <- file.path("output", "quant", meta_ctl_oa$ID, "quant.sf")
colnames(meta_ctl_oa) <- gsub("ID", "names", colnames(meta_ctl_oa))
file.exists(meta_ctl_oa$files)

## Import data with tximeta & summarize to gene
se <- tximeta(meta_ctl_oa)
gse <- summarizeToGene(se)


#-------------------------------------------------------------------------------
#  Gene LEVEL
#-------------------------------------------------------------------------------
## Convert to factors (avoids a warning)
colData(gse)[] <- lapply(colData(gse), factor)

## Build DESeq object
dds_gene <- DESeqDataSet(gse, design = ~ Condition)

## Filter out lowly expressed genes
## (at least 10 counts in at least 2 samples)
keep <- rowSums(counts(dds_gene) >= 10) >= ceiling(nrow(colData(gse))*0.10)
dds_gene <- dds_gene[keep,]

## Fit model
dds_gene <- DESeq(dds_gene)
## Save dds
save(dds_gene, file = "output/quant/deGene_OA_expression_dds.rda")

load("output/quant/deGene_OA_expression_dds.rda")
#------------------------------------------------------------------------------
#Differential gene expression
de_genes_shrink <- lfcShrink(dds_gene,
                             coef = "Condition_OA_vs_CTL", format = "GRanges") |>
  plyranges::names_to_column("gene_id")

## outputs
# res <- results(dds_gene, name ="Condition_WD_vs_KD" )
# summary(res)
res <- lfcShrink(dds_gene, coef="Condition_OA_vs_CTL")
summary(res)
DESeq2::plotMA(res,ylim=c(-2,2),main = "Differential RNAseq CTL vs OA Analysis",
               ylab = "LFC",
               xlab = "mean of norm. counts")

## plot PCA
plotPCA(normTransform(dds_gene), intgroup = "Condition",ntop=23398 ) + ggplot2::theme(aspect.ratio = 1)
plotPCA(vsd, intgroup = c("dex", "cell"))
#plotPCA(normTransform(dds_gene), intgroup = "Condition",ntop=500) + ggplot2::theme(aspect.ratio = 1)

# Join results with gene info
de_genes_shrink <-
  inner_join(x = as.data.frame(de_genes_shrink),
             y = as.data.frame(rowData(gse)) %>%
               dplyr::select(c("gene_id", "tx_ids")),
             by = "gene_id") %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  keepStandardChromosomes(pruning.mode = "coarse") %>%
  as.data.frame()

write_csv(de_genes_shrink, file = "output/quant/de_genes_OA_results.csv")
de_genes_shrink <- read.table("/work/users/s/e/seyoun/CQTL_sQTL/output/quant/de_genes_OA_results.csv",header=T,sep=",") |> as.data.frame()
# Get significant genes
sig_OAdeGenes_pval01_l2fc15 <- de_genes_shrink |> 
  dplyr::filter(padj < 0.01 & abs(log2FoldChange) > 1.5) |> 
  write_csv("/work/users/s/e/seyoun/CQTL_sQTL/output/quant/sig_OAdeGenes_pval01_l2fc15.csv")
sig_OAdeGenes_pval01_l2fc1.2 <- de_genes_shrink |> 
  dplyr::filter(padj < 0.01 & abs(log2FoldChange) > 1.2) |> 
  write_csv("/work/users/s/e/seyoun/CQTL_sQTL/output/quant/sig_OAdeGenes_pval01_l2fc1.2.csv")


# Split into upregulated and downregulated
downsig_OAdeGenes_pval01_l2fc15 <- sig_OAdeGenes_pval01_l2fc15 |> 
  filter(log2FoldChange < 0) |> 
  write_csv("output/quant/downsig_OAdeGenes_pval01_l2fc15.csv")
upsig_OAdeGenes_pval01_l2fc15 <- sig_OAdeGenes_pval01_l2fc15 |> 
  filter(log2FoldChange > 0) |> 
  write_csv("output/quant/upsig_OAdeGenes_pval01_l2fc15.csv")

#-------------------------------------------------------------------------------

vsd_gene <- vst(dds_gene, blind=FALSE)
ntd_gene <- normTransform(dds_gene)
vsd_OAgene_update <- rownames_to_column(as.data.frame(assay(vsd_gene)),"ENSG")
ntd_OAgene_update <- rownames_to_column(as.data.frame(assay(ntd_gene)),"ENSG")

fwrite(vsd_OAgene_update,file="output/quant/normalized_vst_OAgene.txt",sep="\t",quote=F,row.names=F,col.names=T)
fwrite(ntd_OAgene_update,file="output/quant/normalized_log2_OAgene.txt",sep="\t",quote=F,row.names=F,col.names=T)

