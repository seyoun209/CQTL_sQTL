setwd("/work/users/s/e/seyoun/CQTL_sQTL")
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(AnnotationDbi)
library(doParallel)
library(BiocParallel)
library(GeneStructureTools)
library(GenomicRanges)
library(stringr)
library(Gviz)
library(BSgenome)
library(Biostrings)
source("scripts/sQTL_rscripts/utils.R")

# load 
#txdb <- makeTxDbFromGFF(file="/work/users/s/e/seyoun/Ref/genome/gencode.v45.annotation.gtf")
#saveDb(x=txdb, file = "gencode.v45.annotation.TxDb")
txdb <- loadDb("../crispr/02.test_seq/gencode.v45.annotation.TxDb")
exon.re <- exonsBy(txdb,by="tx")
fi.cns <- c("TXCHROM","TXNAME","GENEID","TXSTART","TXEND","TXSTRAND","CDSID","CDSNAME","CDSSTART" ,"CDSEND" ,"EXONID","EXONNAME","EXONSTART","EXONEND")
txTable <- gsub(" ","",as.matrix(select(txdb,keys=names(exon.re),columns=fi.cns,keytype="TXID"))) |> as.data.frame()

FiveUTR <- fiveUTRsByTranscript(txdb)
ThreUTR <- threeUTRsByTranscript(txdb)
CDS.re <- cdsBy(txdb,by="tx")
intron.re <- intronsByTranscript(txdb)

un.CDS.re <- unlist(CDS.re)
un.fiv.re <- unlist(FiveUTR)
un.thr.re <- unlist(ThreUTR)
un.intron.re <- unlist(intron.re)


# MAPK8 (ENSG00000107643.17)
mapk8_txtable <- txTable |> dplyr::filter(GENEID == "ENSG00000107643.17")
mapk8_txtable_selectedTrans <- txTable |> dplyr::filter(TXNAME == "ENST00000374179.8")
txsubsetCDS <- mapk8_txtable_selectedTrans[!is.na(mapk8_txtable_selectedTrans$CDSID),]

response_pbs_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_pbs_resInv.rds")
response_fnf_results <- readRDS("output/01.qtltools_re/conditional_pbs/response_fnf_resInv.rds")



# First, Getting the transcript with the region of interest of mapk8
txAll_ran <- GRanges(
  seqnames = mapk8_txtable_selectedTrans$TXCHROM,
  ranges = IRanges(start = as.numeric(mapk8_txtable_selectedTrans$EXONSTART), end = as.numeric(mapk8_txtable_selectedTrans$EXONEND)),
  strand =mapk8_txtable_selectedTrans$TXSTRAND
)

txAllCDS.ran <- GRanges(
  seqnames = mapk8_txtable$TXCHROM,
  ranges = IRanges(start = as.numeric(mapk8_txtable$EXONSTART), end = as.numeric(mapk8_txtable$EXONEND)),
  strand =mapk8_txtable$TXSTRAND
)

# Most significant intron region MAPK8
mapk8_intron_sig <- response_fnf_results |> dplyr::filter(SYMBOL == "MAPK8")

#Whole cluster of intrested of the region MAPK8

mapk8_cluster_psi <- ratios_fnf[str_detect(rownames(ratios_fnf), mapk8_intron_sig$clusterID), ]
mapk8_cluster_region <- data.frame(str_split(rownames(mapk8_cluster_psi),":",simplify = TRUE))
colnames(mapk8_cluster_region) <- c("chr","start","end","clusterID")
mapk8_cluster_region <- mapk8_cluster_region |> mutate(strand = str_split(mapk8_cluster_region$clusterID,"_",simplify = TRUE)[,3])

intron_interest_ran <- GRanges(
  seqnames = mapk8_intron_sig$phe_chr,
  ranges = IRanges(start = as.numeric(mapk8_intron_sig$phe_from),end = as.numeric(mapk8_intron_sig$phe_to)),
  strand = mapk8_intron_sig$phe_strd
)

as_nt <- width(intron_interest_ran)

cluster_interest_ran <- GRanges(
  seqnames = mapk8_cluster_region$chr,
  ranges = IRanges(start = as.numeric(mapk8_cluster_region$start),end = as.numeric(mapk8_cluster_region$end)),
  strand = mapk8_cluster_region$strand
)



asAll.nm <- as.matrix(findOverlaps(txAllCDS.ran,intron_interest_ran))[1] # position of AS exon among CDS exons

altAll.tx.ran  <- txAllCDS.ran[-asAll.nm]
# AS.nm <- as.matrix(findOverlaps(tx.ran,as.exon.ran))[1] # position of AS exon among CDS exons
# alt.tx.ran <- tx.ran[-AS.nm] # exon ranges after removing AS exon
# 
wAllExon <- getSeq(BSgenome.Hsapiens.UCSC.hg38, txAll_ran)
woAllExon <- getSeq(BSgenome.Hsapiens.UCSC.hg38, altAll.tx.ran)
# wExon <- getSeq(BSgenome.Hsapiens.UCSC.hg38, tx.ran) #getseq
# woExon <- getSeq(BSgenome.Hsapiens.UCSC.hg38, alt.tx.ran) #without AS getseq
start_position <- regexpr("ATG", unlist(wAllExon))[1] #full form where the position starts
start_position0 <- regexpr("ATG", unlist(woAllExon))[1] # alternative from where the position starts
if(start_position0 <0 |start_position<0){next}
# start_position1 <- regexpr("ATG", unlist(wExon))[1]
# start_position2 <- regexpr("ATG", unlist(woExon))[1]

pro.seq <- as.character(translate(unlist(wAllExon)[start_position:length(unlist(wAllExon))])) # protein sequence withexon
if (!length(grep("[*]",pro.seq)))    pro.seq <- paste(pro.seq,"*",sep="")

alt.pro.seq  <- as.character(translate(unlist(woAllExon)[start_position0:length(unlist(woAllExon))]))
if (!length(grep("[*]",alt.pro.seq)))    alt.pro.seq <- paste(alt.pro.seq,"*",sep="")


pro.min <- min(which(strsplit(pro.seq,"")[[1]] == "*")) # position of stop codon in tx.seq
alt.min <- min(which(strsplit(alt.pro.seq,"")[[1]] == "*")) # position of stop codon in alt.tx.seq
ter.diff <- pro.min - alt.min # difference of alt stop codon from ref stop codon

wAllExon_position <- pro.min*3+start_position-3
woAllExon_position <- alt.min *3+start_position0-3

find_row_by_position <- function(position, exon) {
  cumulative_width <- 0
  
  for (i in 1:length(exon)) {
    cumulative_width <- cumulative_width + width(exon[i])
    
    if (cumulative_width >= position) {
      return(i)
    }
  }
  
  return(NA)
}

print_until_stop <- function(sequence) {
  stop_index <- regexpr("\\*", sequence)[1]  # Find the index of the first "*"
  
  if (stop_index == -1) {
    print(sequence)  # No "*", print the entire sequence
  } else {
    print(substr(sequence, 1, stop_index - 1))  # Print until the "*" character
  }
}


# Protein Sequence to check
