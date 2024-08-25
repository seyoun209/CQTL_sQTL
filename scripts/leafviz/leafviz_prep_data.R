library(dplyr)
library(purrr)
library(data.table)
library(leafviz)
library(stringr)
library(readr)
library(gdata)
library(sqtlviztools)
library(foreach)
options(echo=TRUE)



dir <- "/work/users/s/e/seyoun/CQTL_sQTL/output/"
VCF <- paste0(dir,"geno/fnf_geno/06.subset_sigSNps/sigSNPs_all_fnf_vcf.gz")
clusters_table <- paste0(dir,"clu_fnf/ctlvsfnf_perind_numers.counts.gz")
#permutation_res <- response_fnf_results
header <- read.table("scripts/sQTL_rscripts/nominal_header.txt") |> as.character()
permutation_full_re <- fread('output/01.qtltools_re/nominal_fnf/pc4_allchr.fnf.cis')
colnames(permutation_full_re) <- header

clusters <- read.table(clusters_table, header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)


# annotation - created by Leafcutter

annotation_code <- "/work/users/s/e/seyoun/CQTL_sQTL/tools/leafcutter/gencode_hg38"
exon_file <- paste0(annotation_code, "_all_exons.txt.gz")
all_introns <- paste0(annotation_code,"_all_introns.bed.gz" )
threeprime_file <- paste0( annotation_code,"_threeprime.bed.gz")
fiveprime_file <- paste0( annotation_code,"_fiveprime.bed.gz")

exons_table <- if (!is.null( exon_file )) {
  cat("Loading exons from",exon_file,"\n")
  as.data.frame(fread(cmd = paste("zless",exon_file)) )
} else {
  cat("No exon_file provided.\n")
  NULL
}


print("reading in results")


res <- response_fnf_results
genotypes <- unique(res$var_id)

######
# PREPARE VCF GENOTYPES
######

vcf <- fread(file =  VCF, header=TRUE,check.names = FALSE)
vcf[, (ncol(vcf) - 100):ncol(vcf)] <- lapply(vcf[, (ncol(vcf) - 100):ncol(vcf), with = FALSE], function(x) {
  ifelse(x == "0/0", 0,
         ifelse(x == "0/1", 1,
                ifelse(x == "1/1", 2, NA)))
})
vcf_meta <- cf_meta <- data.frame( SNP = vcf$ID, SNP_pos = vcf$POS, REF = vcf$REF, ALT = vcf$ALT, stringsAsFactors = FALSE)



##################
# PREPARE CLUSTERS
##################

# from significant associations
sigClusters <- str_split_fixed(res[,1], ":",4)[,4]
### add Nalls' GWAS SNP clusters and Yang's too!
#sigClusters <- c(sigClusters, gwas_clusters, yang_clusters)

introns <- leafviz::get_intron_meta(row.names(clusters) )
keepClusters <- match(introns$clu,sigClusters)

# remove non-significant (or non-GWAS SNP-associated) clusters
introns <- introns[ !is.na(keepClusters),]
clusters <- clusters[ !is.na(keepClusters),]

# rearrange sample columns in clusters so they match the VCF
samples <- names(vcf)[10:ncol(vcf)]
clusters <- clusters[, samples]

introns_to_plot <- get_intron_meta(row.names(clusters))

# thin down junctions in clusters - remove low contributing junctions

juncProp <- function(cluster){
  cluster$prop <- cluster$meanCount / sum(cluster$meanCount)
  return(cluster)
}

splitClusters <- introns_to_plot %>%
  mutate(
    clu = factor(.$clu, levels = unique(.$clu)),
    meanCount = rowMeans(clusters) ) %>%
  split( .$clu ) %>%
  purrr::map_df( juncProp ) %>%
  mutate( clu = as.character(.$clu))

# thinning out clusters - turn off for now - this needs to be worked on
#thinClusters <- FALSE

if( thinClusters == TRUE){
  introns_to_plot <- introns_to_plot[ splitClusters$prop >= 0.01,]
  clusters <- clusters[ splitClusters$prop >= 0.01,]
  introns <- introns[ splitClusters$prop >= 0.01,]
}
####################
# ANNOTATE JUNCTIONS
####################

intersect_introns <- function(introns){
  all.introns <- introns
  all.introns$start <- as.numeric(all.introns$start)
  all.introns$end <- as.numeric(all.introns$end)
  
  
  # for each splice site write out a bed file
  all.junctions <- dplyr::select(all.introns, chr, start, end, clusterID = clu)
  
  all.fiveprime <- data.frame( chr = all.introns$chr,
                               start = all.introns$start,
                               end = all.introns$start + 1,
                               clusterID = all.introns$clu)
  all.threeprime <- data.frame( chr = all.introns$chr,
                                start = all.introns$end,
                                end = all.introns$end + 1,
                                clusterID = all.introns$clu)
  all.file <- "all_junctions.bed"
  all.fiveprime.file <- "all_fiveprime.bed"
  all.threeprime.file <- "all_threeprime.bed"
  
  write.table( all.junctions, all.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t" )
  write.table( all.threeprime, all.threeprime.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t" )
  write.table( all.fiveprime, all.fiveprime.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t" )
  
  print( "BedTools intersect junctions with list of known splice sites")
  
  # first match junctions
  all.introns.cmd <- paste0("bedtools intersect -a ", all.file, " -b ", all_introns, " -wa -wb -loj -f 1" )
  all.introns_intersect <- data.table::fread("all_introns_intersect", header=FALSE)
  
  # intersect with bedtools to find the annotations of each splice site
  threeprime.cmd <- paste0( "bedtools intersect -a ", all.threeprime.file, " -b ",threeprime_file, " -wa -wb -loj -f 1" )
  threeprime_intersect <- data.table::fread("threeprime_intersect")
  
  fiveprime.cmd <- paste0( "bedtools intersect -a ", all.fiveprime.file, " -b ", fiveprime_file, " -wa -wb -loj -f 1" )
  fiveprime_intersect <- data.table::fread("fiveprime_intersect")
  
  # remove temporary files
  rm.cmd <- paste("rm ", all.file, all.fiveprime.file, all.threeprime.file)
  
  system(rm.cmd)
  
  return( list(threeprime_intersect, fiveprime_intersect,all.introns_intersect))
}

intersects <- list(threeprime_intersect, fiveprime_intersect,all.introns_intersect)

threeprime_intersect <- intersects[[1]]
fiveprime_intersect <- intersects[[2]]
all.introns_intersect <- intersects[[3]]

print("Annotating junctions")

uniqueClusters <- unique( introns$clu )

save.image("debug.RData")


annotatedClusters <- purrr::map_df( seq_along(uniqueClusters),
                                    ~annotate_single_cluster( introns,
                                                              clu = uniqueClusters[.],
                                                              cluIndex = .,
                                                              fiveprime=fiveprime_intersect,
                                                              threeprime=threeprime_intersect,
                                                              bothSS=all.introns_intersect
                                    )
)

annotatedClusters$gene[ is.na( annotatedClusters$gene) ] <- "."
annotatedClusters$ensemblID[ is.na( annotatedClusters$ensemblID) ] <- "."


#################
# PREPARE RESULTS - MOST SIGNIFICANT SNP x JUNCTION
#################

# Bind together metadata with original results
sigJunctions <- cbind( get_intron_meta( res[,"phe_id"]),
                       res[, c("var_id","rsID","FNF_p","var_chr","var_from")])

# sometimes there will be duplicates - remove!
sigJunctions <- dplyr::distinct(sigJunctions)

names(sigJunctions)[8] <- "bpval"
# present most significant junction for each SNP?
# or most significant SNP for each junction?
resultsByCluster <- dplyr::group_by(sigJunctions[order(sigJunctions$bpval),], clu) %>%
  dplyr::summarise( chr = first(chr),
                    start = min(start),
                    end = max(end),
                    var_id = first(var_id),
                    var_chr = first(var_chr),
                    var_from = first(var_from),
                    qval = first(bpval)) %>%
  dplyr::arrange(qval)

colnames(resultsByCluster) <- c("clu","chr","start","end","snp","snp_chr","pos","FDR")


####
## PREPARE FOR SHINY
####
annotation_code <- "/work/users/s/e/seyoun/CQTL_sQTL/tools/leafcutter/gencode_hg38"

resultsByCluster$gene <- annotatedClusters$gene[ match(resultsByCluster$clu, annotatedClusters$clusterID)]
resultsByCluster$SNP_pos <- paste0(resultsByCluster$snp_chr, ":", resultsByCluster$pos)

resultsByCluster$cluster_pos = paste0(resultsByCluster$chr,":", resultsByCluster$start,"-",resultsByCluster$end)

resultsToPlot <- as.data.frame( dplyr::select( resultsByCluster,
                                        SNP = snp,
                                        SNP_pos,
                                        gene = gene,
                                        cluster_pos,
                                        q = FDR 
) )

row.names(resultsToPlot) <- resultsByCluster$clu

resultsToPlot$q <- signif(resultsToPlot$q,  digits = 3)

perm_full <- permutation_full_re |> dplyr::select("phe_id","var_id","dist_phe_var","nom_pval","slope")
names(perm_full) <- c("clusterID", "SNP","X","FDR", "Beta" )

junctionsNeeded <- introns_to_plot %>%
  mutate( cluster = paste( chr, start, end, clu, sep = ":") ) %>%
  pull( cluster )



#perm_clean$RS_id <- snps$RS_id[ match( perm_clean$SNP, snps$dummy2)]
perm_clean <- perm_full %>%
  dplyr::filter( clusterID %in% junctionsNeeded ) %>%
  #mutate( RS_id = snps$RS_id[ match( SNP, snps$dummy2)] ) %>%
  dplyr::filter( !is.na(SNP)) %>%
  dplyr::mutate( SNP_ID = gsub(":", ".",SNP))

# sort out SNP IDs
perm_clean <- dplyr::select(perm_clean, clusterID, SNP2=SNP_ID, Beta, FDR,SNP)

perm_clean <- cbind( perm_clean,
                     get_intron_meta(perm_clean$clusterID),
                     get_snp_meta(perm_clean$SNP2)[,2:3])


perm_clean <- perm_clean |> dplyr::select(-c("SNP2")) |> dplyr::rename(snp_ID = "SNP")

# junction table - each junction with Beta, P value and annotation
junctionTable <- resultsToPlot %>%
  mutate( clu = row.names(resultsToPlot) ) %>%
  left_join(introns_to_plot, by = "clu" ) %>%
  dplyr::rename(snp_ID = SNP) %>%
  left_join( perm_clean,
             by = c("chr" = "snp_chr", "start", "end", "snp_ID", "clu", "middle")
  ) %>%
  mutate(coord = paste0( chr, ":", start, "-", end)) %>%
  left_join( annotatedClusters,
             by = c("clu" = "clusterID", "coord" )
  ) %>%
  dplyr::select(clu, coord, verdict, Beta, q = FDR) %>%
  mutate( Beta = signif(Beta, digits = 3),
          q = signif(q, digits = 3)) %>%
  mutate( Beta = ifelse(is.na(Beta), ".", Beta),
          q = ifelse(is.na(q), ".", q))

# why are some of the q values greater than 1?

#########
## SAVE OBJECTS
#########
# too big

#save.image("example_data.Rdata")
print("saving objects")


save( annotatedClusters, # every junction needed
      sigJunctions, # every junction x SNP interaction
      resultsToPlot, #significant clusters and the most significant SNP

      clusters, # junction counts for each sample
      vcf,# the genotypes of each sample
      vcf_meta, # the vcf metadata
      introns_to_plot, # all the intron positions

      exons_table, # the annotation
      junctionTable, # the junctions to display for each cluster

      annotation_code,
      code,
      file = "/work/users/s/e/seyoun/CQTL_sQTL/fnf_sQTL_results.Rdata"
)

load("/work/users/s/e/seyoun/CQTL_sQTL/sQTL_results.Rdata")


#-------------------------------------------------------------------------------
# Make the PBS sQTL results

dir <- "/work/users/s/e/seyoun/CQTL_sQTL/output/"
VCF <- paste0(dir,"geno/fnf_geno/06.subset_sigSNps/sigSNPs_all_fnf_vcf.gz")
clusters_table <- paste0(dir,"clu_fnf/ctlvsfnf_perind_numers.counts.gz")
#permutation_res <- response_fnf_results
header <- read.table("scripts/sQTL_rscripts/nominal_header.txt") |> as.character()
permutation_full_re <- fread('output/01.qtltools_re/nominal_pbs/pc5_allchr.pbs.cis')
colnames(permutation_full_re) <- header

clusters <- read.table(clusters_table, header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)

annotation_code <- "/work/users/s/e/seyoun/CQTL_sQTL/tools/leafcutter/gencode_hg38"
exon_file <- paste0(annotation_code, "_all_exons.txt.gz")
all_introns <- paste0(annotation_code,"_all_introns.bed.gz" )
threeprime_file <- paste0( annotation_code,"_threeprime.bed.gz")
fiveprime_file <- paste0( annotation_code,"_fiveprime.bed.gz")

exons_table <- if (!is.null( exon_file )) {
  cat("Loading exons from",exon_file,"\n")
  as.data.frame(fread(cmd = paste("zless",exon_file)) )
} else {
  cat("No exon_file provided.\n")
  NULL
}


print("reading in results")


res <- response_pbs_results
genotypes <- unique(res$var_id)


##################
# PREPARE CLUSTERS
##################

# from significant associations
sigClusters <- str_split_fixed(res[,1], ":",4)[,4]
### add Nalls' GWAS SNP clusters and Yang's too!
#sigClusters <- c(sigClusters, gwas_clusters, yang_clusters)

introns <- leafviz::get_intron_meta(row.names(clusters) )
keepClusters <- match(introns$clu,sigClusters)

# remove non-significant (or non-GWAS SNP-associated) clusters
introns <- introns[ !is.na(keepClusters),]
clusters <- clusters[ !is.na(keepClusters),]

# rearrange sample columns in clusters so they match the VCF
samples <- names(vcf)[10:ncol(vcf)]
clusters <- clusters[, samples]

introns_to_plot <- get_intron_meta(row.names(clusters))

# thin down junctions in clusters - remove low contributing junctions

juncProp <- function(cluster){
  cluster$prop <- cluster$meanCount / sum(cluster$meanCount)
  return(cluster)
}

splitClusters <- introns_to_plot %>%
  mutate(
    clu = factor(.$clu, levels = unique(.$clu)),
    meanCount = rowMeans(clusters) ) %>%
  split( .$clu ) %>%
  purrr::map_df( juncProp ) %>%
  mutate( clu = as.character(.$clu))

####################
# ANNOTATE JUNCTIONS
####################


intersects <- list(threeprime_intersect, fiveprime_intersect,all.introns_intersect)

threeprime_intersect <- intersects[[1]]
fiveprime_intersect <- intersects[[2]]
all.introns_intersect <- intersects[[3]]

print("Annotating junctions")

uniqueClusters <- unique( introns$clu )

save.image("debug.RData")


annotatedClusters <- purrr::map_df( seq_along(uniqueClusters),
                                    ~annotate_single_cluster( introns,
                                                              clu = uniqueClusters[.],
                                                              cluIndex = .,
                                                              fiveprime=fiveprime_intersect,
                                                              threeprime=threeprime_intersect,
                                                              bothSS=all.introns_intersect
                                    )
)

annotatedClusters$gene[ is.na( annotatedClusters$gene) ] <- "."
annotatedClusters$ensemblID[ is.na( annotatedClusters$ensemblID) ] <- "."

#################
# PREPARE RESULTS - MOST SIGNIFICANT SNP x JUNCTION
#################

# Bind together metadata with original results
sigJunctions <- cbind( get_intron_meta( res[,"phe_id"]),
                       res[, c("var_id","rsID","PBS_p","var_chr","var_from")])

# sometimes there will be duplicates - remove!
sigJunctions <- dplyr::distinct(sigJunctions)

names(sigJunctions)[8] <- "bpval"
# present most significant junction for each SNP?
# or most significant SNP for each junction?
resultsByCluster <- dplyr::group_by(sigJunctions[order(sigJunctions$bpval),], clu) %>%
  dplyr::summarise( chr = first(chr),
                    start = min(start),
                    end = max(end),
                    var_id = first(var_id),
                    var_chr = first(var_chr),
                    var_from = first(var_from),
                    qval = first(bpval)) %>%
  dplyr::arrange(qval)

colnames(resultsByCluster) <- c("clu","chr","start","end","snp","snp_chr","pos","FDR")


####
## PREPARE FOR SHINY
####
annotation_code <- "/work/users/s/e/seyoun/CQTL_sQTL/tools/leafcutter/gencode_hg38"

resultsByCluster$gene <- annotatedClusters$gene[ match(resultsByCluster$clu, annotatedClusters$clusterID)]
resultsByCluster$SNP_pos <- paste0(resultsByCluster$snp_chr, ":", resultsByCluster$pos)

resultsByCluster$cluster_pos = paste0(resultsByCluster$chr,":", resultsByCluster$start,"-",resultsByCluster$end)

resultsToPlot <- as.data.frame( dplyr::select( resultsByCluster,
                                               SNP = snp,
                                               SNP_pos,
                                               gene = gene,
                                               cluster_pos,
                                               q = FDR 
) )

row.names(resultsToPlot) <- resultsByCluster$clu

resultsToPlot$q <- signif(resultsToPlot$q,  digits = 3)

perm_full <- permutation_full_re |> dplyr::select("phe_id","var_id","dist_phe_var","nom_pval","slope")
names(perm_full) <- c("clusterID", "SNP","X","FDR", "Beta" )

junctionsNeeded <- introns_to_plot %>%
  mutate( cluster = paste( chr, start, end, clu, sep = ":") ) %>%
  pull( cluster )



#perm_clean$RS_id <- snps$RS_id[ match( perm_clean$SNP, snps$dummy2)]
perm_clean <- perm_full %>%
  dplyr::filter( clusterID %in% junctionsNeeded ) %>%
  #mutate( RS_id = snps$RS_id[ match( SNP, snps$dummy2)] ) %>%
  dplyr::filter( !is.na(SNP)) %>%
  dplyr::mutate( SNP_ID = gsub(":", ".",SNP))

# sort out SNP IDs
perm_clean <- dplyr::select(perm_clean, clusterID, SNP2=SNP_ID, Beta, FDR,SNP)

perm_clean <- cbind( perm_clean,
                     get_intron_meta(perm_clean$clusterID),
                     get_snp_meta(perm_clean$SNP2)[,2:3])


perm_clean <- perm_clean |> dplyr::select(-c("SNP2")) |> dplyr::rename(snp_ID = "SNP")

# junction table - each junction with Beta, P value and annotation
junctionTable <- resultsToPlot %>%
  mutate( clu = row.names(resultsToPlot) ) %>%
  left_join(introns_to_plot, by = "clu" ) %>%
  dplyr::rename(snp_ID = SNP) %>%
  left_join( perm_clean,
             by = c("chr" = "snp_chr", "start", "end", "snp_ID", "clu", "middle")
  ) %>%
  mutate(coord = paste0( chr, ":", start, "-", end)) %>%
  left_join( annotatedClusters,
             by = c("clu" = "clusterID", "coord" )
  ) %>%
  dplyr::select(clu, coord, verdict, Beta, q = FDR) %>%
  mutate( Beta = signif(Beta, digits = 3),
          q = signif(q, digits = 3)) %>%
  mutate( Beta = ifelse(is.na(Beta), ".", Beta),
          q = ifelse(is.na(q), ".", q))

# why are some of the q values greater than 1?

#########
## SAVE OBJECTS
#########
# too big

#save.image("example_data.Rdata")
print("saving objects")


save( annotatedClusters, # every junction needed
      sigJunctions, # every junction x SNP interaction
      resultsToPlot, #significant clusters and the most significant SNP
      
      clusters, # junction counts for each sample
      vcf,# the genotypes of each sample
      vcf_meta, # the vcf metadata
      introns_to_plot, # all the intron positions
      
      exons_table, # the annotation
      junctionTable, # the junctions to display for each cluster
      
      annotation_code,
      code,
      file = "/work/users/s/e/seyoun/CQTL_sQTL/pbs_sQTL_results.Rdata"
)

load("/work/users/s/e/seyoun/CQTL_sQTL/sQTL_results.Rdata")

