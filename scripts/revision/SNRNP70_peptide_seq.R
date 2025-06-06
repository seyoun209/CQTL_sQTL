# Finding the peptide sequence for SNRNP70

setwd("/work/users/s/e/seyoun/CQTL_sQTL/")

library(Biostrings)
library(AnnotationDbi)
library(BSgenome.Hsapiens.UCSC.hg38)
#-------------------------------------------------------------------------------
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


# SNRNP70 gene ID ENSG00000104852.15

snrnp70_txtable <- txTable |> dplyr::filter(GENEID == "ENSG00000104852.15")


# Transcript alternative exon 8 located it
snrnp70_txtable_alt8 <- txTable |> dplyr::filter(TXNAME == "ENST00000401730.5")
snrnp70_alt8_txsubsetCDS <- snrnp70_txtable_alt8[!is.na(snrnp70_txtable_alt8$CDSID),]

# First, Getting the transcript with the region of interest of SNRNP70 alternative 8 transcrip

snrnp70_alt8_All_ran <- GRanges(
  seqnames = snrnp70_txtable_alt8$TXCHROM,
  ranges = IRanges(start = as.numeric(snrnp70_txtable_alt8$EXONSTART), end = as.numeric(snrnp70_txtable_alt8$EXONEND)),
  strand =snrnp70_txtable_alt8$TXSTRAND
)

# This is the only CDS of the SNRNP70 alternative 8 Transcript
snrnp70_alt8_CDS_ran <- GRanges(
  seqnames = snrnp70_alt8_txsubsetCDS$TXCHROM,
  ranges = IRanges(start = as.numeric(snrnp70_alt8_txsubsetCDS$EXONSTART), end = as.numeric(snrnp70_alt8_txsubsetCDS$EXONEND)),
  strand =snrnp70_alt8_txsubsetCDS$TXSTRAND
)

# alternative exon 8 region SNRNP70

alt8_exon <- GRanges(
  seqnames = "chr19",
  ranges = IRanges(start = 49102114,end = 49103587),
  strand = "+"
)

as_exon_nt <- width(alt8_exon)

as_all.ma <- as.matrix(findOverlaps(snrnp70_alt8_All_ran,alt8_exon))[1]

snrnp70_alt8out_tx.ran  <- snrnp70_alt8_All_ran[-as_all.ma]

wAllExon <- getSeq(BSgenome.Hsapiens.UCSC.hg38, snrnp70_alt8_All_ran)
woAllExon <- getSeq(BSgenome.Hsapiens.UCSC.hg38, snrnp70_alt8out_tx.ran)

start_position <- regexpr("ATG", unlist(wAllExon))[1] #full form where the position starts
start_position0 <- regexpr("ATG", unlist(woAllExon))[1] # alternative from where the position starts

pro.seq <- as.character(translate(unlist(wAllExon)[start_position:length(unlist(wAllExon))])) # protein sequence withexon
alt.pro.seq  <- as.character(translate(unlist(woAllExon)[start_position0:length(unlist(woAllExon))]))

