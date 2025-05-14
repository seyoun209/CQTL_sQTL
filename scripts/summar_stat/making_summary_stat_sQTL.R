# Summary statistics getting ready (PBS and FNF)

setwd("/work/users/s/e/seyoun/CQTL_sQTL/")
library(data.table)
library(GenomicRanges)
suppressMessages(library(qvalue))
source("./scripts/sQTL_rscripts/utils.R")

norm_header <-  read.table("scripts/sQTL_rscripts/nominal_header.txt") 
#-------------------------------------------------------------------------------
#conditional  formatting
cond_header <- c("phe_id","phe_chr","phe_from",
                 "phe_to","phe_strd","n_var_in_cis","dist_phe_var","var_id",
                 "var_chr","var_from","var_to","rank",
                 "fwd_pval","fwd_r_squared","fwd_slope","fwd_best_hit","fwd_sig",
                 "bwd_pval","bwd_r_squared","bwd_slope","bwd_best_hit","bwd_sig")
conditional <- fread(paste0("output/01.qtltools_re/conditional_pbs/",i,"_pbs_condtional.txt"))
colnames(conditional) <- cond_header

pbs_QTL_cond_final <- readRDS("output/01.qtltools_re/pbs_cond_sig_beta_maf_added.rds")
fnf_QTL_cond_final <- readRDS("output/01.qtltools_re/fnf_cond_sig_beta_maf_added.rds")

pbs_cond_values <- pbs_QTL_cond_final %>%
  select(phe_id, var_id,rank) %>%
  distinct()

fnf_cond_values <- fnf_QTL_cond_final %>%
  select(phe_id, var_id,rank) %>%
  distinct()

#-------------------------------------------------------------------------------
#permuation formatting

perm_header <- read.table("scripts/sQTL_rscripts/header.txt")

perm <- fread(paste0("output/01.qtltools_re/perm_pbs/pc5/",i,".pbs.perm"))
colnames(perm) <- paste(perm_header)
perm_noNA <- perm |> filter(!is.na(var_id))

sig_perm_header <- cbind(perm_header, "perm_feature_qval","pval_nominal_threshold")

permQTL_fnf <- read.table("output/01.qtltools_re/01.significant/fnf_0.05_pc4.significant.txt") |> as.data.frame()
permQTL_pbs <- read.table("output/01.qtltools_re/01.significant/pbs_0.05_pc5.significant.txt") |> as.data.frame()
colnames(permQTL_pbs) <- sig_perm_header
colnames(permQTL_fnf) <- sig_perm_header

pbs_perm_values <- permQTL_pbs %>%
  select(phe_id, adj_beta_pval,perm_feature_qval, pval_nominal_threshold) %>%
  distinct()

fnf_perm_values <- permQTL_fnf %>%
  select(phe_id, adj_beta_pval, perm_feature_qval, pval_nominal_threshold) %>%
  distinct()

all_permQTL_fnf <- read.table("output/01.qtltools_re/01.significant/fnf_1_pc4.significant.txt") |> as.data.frame()
all_permQTL_pbs <- read.table("output/01.qtltools_re/01.significant/pbs_1_pc5.significant.txt") |> as.data.frame()
colnames(all_permQTL_fnf) <- sig_perm_header
colnames(all_permQTL_pbs) <- sig_perm_header


all_pbs_perm_values <- all_permQTL_pbs %>%
  select(phe_id, adj_beta_pval,perm_feature_qval) %>%
  distinct()

all_fnf_perm_values <- all_permQTL_fnf %>%
  select(phe_id, adj_beta_pval, perm_feature_qval) %>%
  distinct()




# Preparing for the MAF
maf_all <- fread("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/freq/cqtl.frq")
maf_subset <-  maf_all %>% dplyr::select(c("SNP","MAF","A1","A2")) %>% 
  dplyr::rename("var_id" =SNP, "Effect allele" =A1, "Other allele" = A2) 

# Function to read and process BED files
read_bed_file <- function(file_path) {
  fread(file_path, col.names = c("chr", "start", "var_id", "ref", "alt", "rsID"))
}

# Get list of BED files
bed_files <- list.files("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/rsid",pattern = "_rsid.bed",full.names = TRUE)
all_bed_data <- rbindlist(lapply(bed_files, read_bed_file))
all_bed_dataID <- all_bed_data |> dplyr::select(c("var_id","rsID"))

#merge rsIDs
maf_id_add <- left_join(maf_subset , all_bed_dataID,by="var_id")

#save(maf_id_add, file="/work/users/s/e/seyoun/CQTL_sQTL/output/geno/freq/maf_id_add_v2.rds")
load("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/freq/maf_id_add_v2.rds")


#Annotate Gene
exons <- read.delim(gzfile("/work/users/s/e/seyoun/CQTL_sQTL/tools/leafcutter/gencode_hg38_all_exons.txt.gz"), header=TRUE)
introns <- read.delim(gzfile("/work/users/s/e/seyoun/CQTL_sQTL/tools/leafcutter/gencode_hg38_all_introns.bed.gz"), 
                      header=FALSE,
                      col.names=c("chr","start","end","gene_name","gene_id","strand","transcript_id","number","type","annotations"))

fiveprime <- read.delim(gzfile("/work/users/s/e/seyoun/CQTL_sQTL/tools/leafcutter/gencode_hg38_fiveprime.bed.gz"), 
                        header=FALSE,
                        col.names=c("chr","start","end","gene_name","gene_id","strand","transcript_id","number","type","annotations"))

threeprime <- read.delim(gzfile("/work/users/s/e/seyoun/CQTL_sQTL/tools/leafcutter/gencode_hg38_threeprime.bed.gz"), 
                         header=FALSE,
                         col.names=c("chr","start","end","gene_name","gene_id","strand","transcript_id","number","type","annotations"))

all_annotations <- rbind(
  data.frame(introns[,c("chr","start","end","gene_id")], source="intron"),
  #data.frame(exons[,c("chr","start","end","gene_name")], source="exon"),
  data.frame(fiveprime[,c("chr","start","end","gene_id")], source="5prime"),
  data.frame(threeprime[,c("chr","start","end","gene_id")], source="3prime")
)

# Create GRanges object once for all annotations
annot_gr <- GRanges(all_annotations$chr, 
                    IRanges(all_annotations$start, all_annotations$end), 
                    gene_id=all_annotations$gene_id)

load("/work/users/s/e/seyoun/CQTL_sQTL/output/geno/freq/maf_id_add.rds")
chr <- paste0("chr",seq(1:22))

#-------------------------------------------------------------------------------
#PBS summary stat
# Create empty list to store data frames
all_chr_data <- list()

# Loop through chromosomes
for (i in chr) {
  nominal <- fread(paste0("output/01.qtltools_re/nominal_pbs/pc5/",i,".pbs.cis"))
  colnames(nominal) <- paste(norm_header)
  
  nominal_ensg <- annotate_genes(nominal, hg38_intron_sub_select_first, leafcutter_pheno_subset, txdb_genes)
  
  nominal_update <- nominal_ensg %>%
    left_join(maf_id_add, by="var_id") %>%
    left_join(all_pbs_perm_values, by="phe_id") %>%
    left_join(pbs_cond_values, by = c("phe_id", "var_id"))
  
  nominal_final <- nominal_update %>% select(
    Variant = var_id,
    Effect_allele = `Effect allele`,
    Other_allele = `Other allele`, 
    Feature = ensg,
    FeatureCoordinates = genomicLoc,
    `p-value` = nom_pval,
    beta = r_squared,
    se = slope_se,
    ID = phe_id,
    perm_adj_beta_pval = adj_beta_pval,
    perm_feature_qval = perm_feature_qval,
    conditional_signal = rank,
    EAF = MAF
  )
  
  # Add to list
  all_chr_data[[i]] <- nominal_final
}

# Combine all chromosomes
pbs_combined_data <- do.call(rbind, all_chr_data)

# Save combined data
fwrite(pbs_combined_data, "/work/users/s/e/seyoun/dbGap/sQTL/MSK/CHON_sQTL_PBS_summarystats.csv")

    
#-------------------------------------------------------------------------------
#FNF
    
# Create empty list to store data frames
all_chr_data_fnf <- list()

# Loop through chromosomes
for (i in chr) {
  nominal <- fread(paste0("output/01.qtltools_re/nominal_fnf/pc4/",i,".fnf.cis"))
  colnames(nominal) <- paste(norm_header)
  
  nominal_ensg <- annotate_genes(nominal, hg38_intron_sub_select_first, leafcutter_pheno_subset, txdb_genes)
  
  nominal_update <- nominal_ensg %>%
    left_join(maf_id_add, by="var_id") %>%
    left_join(all_fnf_perm_values, by="phe_id") %>%
    left_join(fnf_cond_values, by = c("phe_id", "var_id"))
  
  nominal_final <- nominal_update %>% select(
    Variant = var_id,
    Effect_allele = `Effect allele`,
    Other_allele = `Other allele`, 
    Feature = ensg,
    FeatureCoordinates = genomicLoc,
    `p-value` = nom_pval,
    beta = r_squared,
    se = slope_se,
    ID = phe_id,
    perm_adj_beta_pval = adj_beta_pval,
    perm_feature_qval = perm_feature_qval,
    conditional_signal = rank,
    EAF = MAF
  )
  
  # Add to list
  all_chr_data_fnf[[i]] <- nominal_final
}

# Combine all chromosomes
fnf_combined_data <- do.call(rbind, all_chr_data_fnf)

# Save combined data
fwrite(fnf_combined_data, "/work/users/s/e/seyoun/dbGap/sQTL/MSK/CHON_sQTL_FNF_summarystats.csv")



