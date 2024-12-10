# Summary statistics getting ready (PBS and FNF)

setwd("/work/users/s/e/seyoun/CQTL_sQTL/")
library(data.table)
library(GenomicRanges)

norm_header <-  read.table("scripts/sQTL_rscripts/nominal_header.txt") 
cond_header <- c("phe_id","phe_chr","phe_from",
                 "phe_to","phe_strd","n_var_in_cis","dist_phe_var","var_id",
                 "var_chr","var_from","var_to","rank",
                 "fwd_pval","fwd_r_squared","fwd_slope","fwd_best_hit","fwd_sig",
                 "bwd_pval","bwd_r_squared","bwd_slope","bwd_best_hit","bwd_sig")

perm_header <- read.table("scripts/sQTL_rscripts/header.txt")


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
for (i in chr) {
  nominal <- fread(paste0("output/01.qtltools_re/nominal_pbs/pc5/",i,".pbs.cis"))
  colnames(nominal) <- paste(norm_header)
  nom_gr <- GRanges(nominal$phe_chr, IRanges(nominal$phe_from, nominal$phe_to))
  overlaps <- findOverlaps(nom_gr, annot_gr)
  gene_mapping <- data.frame(
    idx = queryHits(overlaps),
    gene_id = annot_gr$gene_id[subjectHits(overlaps)]
  )
  gene_mapping <- gene_mapping[!duplicated(gene_mapping$idx), ]
  nominal$gene_id <- NA
  nominal$gene_id[gene_mapping$idx] <- gene_mapping$gene_id
  nominal$FeatureCoordinates <- paste0(nominal$phe_chr, ":", nominal$phe_from, "-", nominal$phe_to)
  
  nominal_update <- nominal %>% left_join(maf_id_add, by="var_id") %>%
    left_join(conditional,by=c("phe_id" = "phe_id", "var_id" = "var_id"))
  perm <- fread(paste0("output/01.qtltools_re/perm_pbs/pc5/",i,".pbs.perm"))
  colnames(perm) <- paste(perm_header)
  conditional <- fread(paste0("output/01.qtltools_re/conditional_pbs/",i,"_pbs_condtional.txt"))
  colnames(conditional) <- cond_header
  
}



# Add gene column to nominal data
nominal$gene_id <- NA
nominal$gene_id[gene_mapping$idx] <- gene_mapping$gene_id
pbs_chr22 <- fread("/work/users/s/e/seyoun/dbGap/nicole/chr22_filtered.txt")
"output/01.qtltools_re/nominal_pbs/pc5/chr1.pbs.cis"
