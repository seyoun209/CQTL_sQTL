#annotate_single_cluster ####
annotate_single_cluster <- function(introns, clu, cluIndex, fiveprime, threeprime, bothSS){
  #print(clu)
  # for each intron in the cluster, check for coverage of both
  # output a vector of string descriptions
  cluster <- introns[introns$clu == clu, ]
  cluster$start <- as.integer(cluster$start)
  cluster$end <- as.integer(cluster$end)
  
  # Check if fiveprime, threeprime, and bothSS exist
  # Check if fiveprime, threeprime, and bothSS exist and if the required column exists
  fprimeClu <- if (!is.null(fiveprime) && "V4" %in% colnames(fiveprime)) {
    dplyr::filter(fiveprime, V4 == clu)
  } else {
    data.frame()
  }
  
  tprimeClu <- if (!is.null(threeprime) && "V4" %in% colnames(threeprime)) {
    dplyr::filter(threeprime, V4 == clu)
  } else {
    data.frame()
  }
  
  bothSSClu <- if (!is.null(bothSS) && "V4" %in% colnames(bothSS)) {
    dplyr::filter(bothSS, V4 == clu)
  } else {
    data.frame()
  }
  
  # for each intron in the cluster:
  #   create vector of overlapping splice sites, indexed by the row of the intersect
  if (nrow(cluster) == 1) {
    fprime <- list(fprimeClu)
    tprime <- list(tprimeClu)
    bothSS <- list(bothSSClu)
    cluster_genes <- unique(c(tprime[[1]]$V8, fprime[[1]]$V8, bothSS[[1]]$V8))
    cluster_ensemblIDs <- unique(c(tprime[[1]]$V9, fprime[[1]]$V9, bothSS[[1]]$V9))
  } else {
    fprime <- apply(cluster, MAR = 1, FUN = function(x) {
      chr <- which(names(cluster) == "chr")
      start <- which(names(cluster) == "start")
      if (!is.null(fprimeClu) && nrow(fprimeClu) > 0) {
        dplyr::filter(fprimeClu, V1 == x[chr] & V2 == as.numeric(x[start]))
      } else {
        data.frame()
      }
    })
    
    tprime <- apply(cluster, MAR = 1, FUN = function(x) {
      chr <- which(names(cluster) == "chr")
      end <- which(names(cluster) == "end")
      if (!is.null(tprimeClu) && nrow(tprimeClu) > 0) {
        dplyr::filter(tprimeClu, V1 == x[chr] & V2 == as.numeric(x[end]))
      } else {
        data.frame()
      }
    })
    
    bothSS <- apply(cluster, MAR = 1, FUN = function(x) {
      chr <- which(names(cluster) == "chr")
      start <- which(names(cluster) == "start")
      end <- which(names(cluster) == "end")
      if (!is.null(bothSSClu) && nrow(bothSSClu) > 0) {
        dplyr::filter(bothSSClu, V6 == as.numeric(x[start]) & V7 == as.numeric(x[end]))
      } else {
        data.frame()
      }
    })
    
    cluster_genes <- names(sort(table(do.call(what = rbind, c(tprime, fprime, bothSS))$V8), decreasing = TRUE))
    cluster_ensemblIDs <- names(sort(table(do.call(what = rbind, c(tprime, fprime, bothSS))$V9), decreasing = TRUE))
  }
  
  cluster_gene <- cluster_genes[cluster_genes != "."][1]
  if (length(cluster_gene) == 0) {
    cluster_gene <- "."
  }
  
  cluster_ensemblID <- cluster_ensemblIDs[cluster_ensemblIDs != "."][1]
  if (length(cluster_ensemblID) == 0) {
    cluster_ensemblID <- "."
  }
  
  verdict <- c()
  coord <- c()
  gene <- c()
  ensemblID <- c()
  
  for (intron in 1:nrow(cluster)) {
    coord[intron] <- paste0(cluster[intron,]$chr, ":", cluster[intron,]$start, "-", cluster[intron,]$end)
    
    gene[intron] <- cluster_gene
    ensemblID[intron] <- cluster_ensemblID
    
    verdict[intron] <- "error"
    if (all(tprime[[intron]]$V5 == ".") & all(fprime[[intron]]$V5 == ".")) {
      verdict[intron] <- "cryptic_unanchored"
    }
    if (all(tprime[[intron]]$V5 == ".") & all(fprime[[intron]]$V5 != ".")) {
      verdict[intron] <- "cryptic_threeprime"
    }
    if (all(tprime[[intron]]$V5 != ".") & all(fprime[[intron]]$V5 == ".")) {
      verdict[intron] <- "cryptic_fiveprime"
    }
    if (all(tprime[[intron]]$V5 != ".") & all(fprime[[intron]]$V5 != ".")) {
      if (all(bothSS[[intron]]$V5 != ".")) {
        verdict[intron] <- "annotated"
      } else {
        verdict[intron] <- "novel annotated pair"
      }
    }
  }
  
  if (cluIndex %% 500 == 0) {
    print(paste("processed", cluIndex, "clusters"))
  }
  
  return(
    data.frame(
      clusterID = clu,
      coord = coord,
      gene = gene,
      ensemblID = ensemblID,
      verdict = verdict,
      stringsAsFactors = FALSE
    )
  )
}

